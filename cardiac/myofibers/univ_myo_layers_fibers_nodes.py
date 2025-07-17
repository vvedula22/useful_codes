import os
import sys
import vtk
import time
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from datetime import timedelta

EPS = sys.float_info.epsilon
PI  = np.pi

DATASTR1 = 'Phi_EPI'
DATASTR2 = 'Phi_AB'

#----------------------------------------------------------------------
# Define fiber and sheet angles
# Fiber angle at endocardium
ALFA_END = (+40.0)*PI/180.0
# Fiber angle at epicardium
ALFA_EPI = (-50.0)*PI/180.0
# Sheet angle at endocardium
BETA_END = (-65.0)*PI/180.0
# Sheet angle at epicardium
BETA_EPI = (+25.0)*PI/180.0
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def isZero(a, b=None):

    absA = abs(a)
    if b == None:
        b = 0.0

    absA = abs(a)
    absB = abs(b)
    if absB > absA:
        absA, absB = absB, absA

    nrm  = max(absA, absB)
    flag = False
    if ((absA-absB)/nrm) < (10.0*EPS):
        flag = True

    return flag
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def norm(u, v=None):

    if v == None:
        v = u

    n = np.size(u)
    l2norm = 0.0
    for i in range(0, n):
        l2norm = l2norm + (u[i]*v[i])

    return l2norm
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def createVTKDataArray(dType, numComp, numTupls, dName):

    if dType == "double":
        D = vtk.vtkDoubleArray()
    elif dType == "int":
        D = vtk.vtkIntArray()

    D.SetNumberOfComponents(numComp)
    D.Allocate(numTupls)
    D.SetNumberOfTuples(numTupls)
    D.SetName(dName)

    return D
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def loadLaplaceSoln(fileName):

    print ("   Loading Laplace solution   <---   %s" % (fileName))
    vtuReader = vtk.vtkXMLUnstructuredGridReader()
    vtuReader.SetFileName(fileName)
    vtuReader.Update()

    vtuMesh = vtk.vtkUnstructuredGrid()

    print ("   Extracting solution and its gradients at points")

    gradFilter = vtk.vtkGradientFilter()
    gradFilter.SetInputConnection(vtuReader.GetOutputPort())

    gradFilter.SetInputScalars(vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS,\
        DATASTR1)
    gradFilter.SetResultArrayName(DATASTR1 + '_pgrad')
    gradFilter.Update()
    vtuMesh = gradFilter.GetOutput()
    nPhiEP  = vtuMesh.GetPointData().GetArray(DATASTR1)
    nGPhiEP = vtuMesh.GetPointData().GetArray(DATASTR1+'_pgrad')

    gradFilter.SetInputScalars(vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS,\
        DATASTR2)
    gradFilter.SetResultArrayName(DATASTR2 + '_pgrad')
    gradFilter.Update()
    vtuMesh = gradFilter.GetOutput()
    nPhiAB  = vtuMesh.GetPointData().GetArray(DATASTR2)
    nGPhiAB = vtuMesh.GetPointData().GetArray(DATASTR2+'_pgrad')

    # Clean unnecessary arrays
    vtuMesh.GetPointData().RemoveArray(DATASTR1)
    vtuMesh.GetPointData().RemoveArray(DATASTR1+'_grad')

    vtuMesh.GetPointData().RemoveArray(DATASTR2)
    vtuMesh.GetPointData().RemoveArray(DATASTR2+'_grad')

    return vtuMesh, nPhiEP, nPhiAB, nGPhiEP, nGPhiAB
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def loadLaplaceSolnCells(fileName):

    print ("   Loading Laplace solution at cells   <---   %s" % (fileName))
    vtuReader = vtk.vtkXMLUnstructuredGridReader()
    vtuReader.SetFileName(fileName)
    vtuReader.Update()

    vtuMesh = vtk.vtkUnstructuredGrid()

    print ("   Extracting solution at cells")

    pt2Cell = vtk.vtkPointDataToCellData()
    pt2Cell.SetInputConnection(vtuReader.GetOutputPort())
    pt2Cell.PassPointDataOn()

    gradFilter = vtk.vtkGradientFilter()
    gradFilter.SetInputConnection(pt2Cell.GetOutputPort())

    gradFilter.SetInputScalars(vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS,\
        DATASTR1)
    gradFilter.SetResultArrayName(DATASTR1 + '_grad')
    gradFilter.Update()
    vtuMesh = gradFilter.GetOutput()
    cPhiEP  = vtuMesh.GetCellData().GetArray(DATASTR1)
    cGPhiEP = vtuMesh.GetCellData().GetArray(DATASTR1+'_grad')

    gradFilter.SetInputScalars(vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS,\
        DATASTR2)
    gradFilter.SetResultArrayName(DATASTR2 + '_grad')
    gradFilter.Update()
    vtuMesh = gradFilter.GetOutput()
    cPhiAB  = vtuMesh.GetCellData().GetArray(DATASTR2)
    cGPhiAB = vtuMesh.GetCellData().GetArray(DATASTR2+'_grad')

    # Clean unnecessary arrays
    vtuMesh.GetPointData().RemoveArray(DATASTR1)
    vtuMesh.GetCellData().RemoveArray(DATASTR1)
    vtuMesh.GetCellData().RemoveArray(DATASTR1+'_grad')

    vtuMesh.GetPointData().RemoveArray(DATASTR2)
    vtuMesh.GetCellData().RemoveArray(DATASTR2)
    vtuMesh.GetCellData().RemoveArray(DATASTR2+'_grad')

    return cPhiEP, cPhiAB
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def setDomainID(vtuMesh, vtuPhiEP, vtuPhiAB):

    print ("   Setting domain IDs at cells")
    numCell = vtuMesh.GetNumberOfCells()
    dmnIDs  = createVTKDataArray("int", 1, numCell, "DOMAIN_ID")

    for iCell in range(0, numCell):
        phiEP = vtuPhiEP.GetTuple1(iCell)
        phiAB = vtuPhiAB.GetTuple1(iCell)

        dmnID = 2
        if phiEP >= EPI_MYO:
            dmnID = 3
        if (phiEP <= END_MYO):
            dmnID = 1
        if phiAB >= MID_BAS:
            dmnID = 0
        dmnIDs.SetTuple1(iCell, dmnID)

    return dmnIDs
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def getAlfaW(d):
    return ALFA_END*(1-d) + ALFA_EPI*d
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def getBetaW(d):
    return BETA_END*(1-d) + BETA_EPI*d
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def axis (u, v):

    n1 = np.sqrt(norm(u))
    if not isZero(n1):
        e1 = u / n1
    else:
        e1 = [1.0, 0.0, 0.0]
        print (" Zero gradient norm detected (e1)")

    e2 = v - np.dot(e1, v)*e1
    n2 = np.sqrt(norm(e2))
    if not isZero(n2):
        e2 = e2 / n2
    else:
        e2 = [0.0, 1.0, 0.0]
        print (" Zero gradient norm detected (e2)")

    e0 = np.cross(e1, e2)

    Q  = np.zeros((3,3))
    Q[:,0] = e0
    Q[:,1] = e1
    Q[:,2] = e2

    return Q
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def orient(Q, alpha, beta):

    ca = np.cos(alpha)
    sa = np.sin(alpha)
    cb = np.cos(beta)
    sb = np.sin(beta)

    Ra = np.array([ [ ca,  -sa,  0.0],
                    [ sa,   ca,  0.0],
                    [0.0,  0.0,  1.0]])

    Rb = np.array([ [1.0,  0.0,  0.0],
                    [0.0,   cb,  -sb],
                    [0.0,   sb,   cb]])

    Qt = np.matmul(Q, np.matmul(Ra, Rb))

    return Qt
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def getFiberDirections(vtuMesh, vtuPhiEP, vtuGPhiEP, vtuGPhiAB):

    numPoints = vtuMesh.GetNumberOfPoints()
    print ("   Computing fiber directions at points")
    F = createVTKDataArray("double", 3, numPoints, "FIB_DIR1")
    T = createVTKDataArray("double", 3, numPoints, "FIB_DIR2")
    S = createVTKDataArray("double", 3, numPoints, "FIB_DIR3")

    j = 1
    k = 1
    print ("      Progress "),
    sys.stdout.flush()
    for iPt in range(0, numPoints):
        phiEP  = vtuPhiEP.GetTuple1(iPt)
        gPhiEP = vtuGPhiEP.GetTuple3(iPt)
        gPhiAB = vtuGPhiAB.GetTuple3(iPt)

        alfaW = getAlfaW(phiEP)
        betaW = getBetaW(phiEP)

        Q_EPI = axis([gPhiAB[0], gPhiAB[1], gPhiAB[2]], \
            [ gPhiEP[0],  gPhiEP[1],  gPhiEP[2]])
        FST = orient(Q_EPI, alfaW, betaW)

        F.SetTuple3(iPt, FST[0,0], FST[1,0], FST[2,0])
        S.SetTuple3(iPt, FST[0,1], FST[1,1], FST[2,1])
        T.SetTuple3(iPt, FST[0,2], FST[1,2], FST[2,2])
        if iPt==j:
            print ("%d%%  " % ((k-1)*10)),
            sys.stdout.flush()
            k = k + 1
            j = int(float((k-1)*numPoints)/10.0)
    print ("[Done!]")

    return F, T, S
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# Main function
if __name__ == '__main__':

    t1 = time.perf_counter()
    print ("========================================================")
    fileName = "results_laplace/result_040.vtu"
    vtuMesh, phiEP, phiAB, gPhiEP, gPhiAB = loadLaplaceSoln(fileName)

    # Fiber directions use convention from Bayer et al.
    #    F: longitudinal fiber direction
    #    T: transverse or sheet orientation
    #    S: sheet-normal direction
    F, T, S = getFiberDirections(vtuMesh, phiEP, gPhiEP, gPhiAB)

    cphiEP, cphiAB = loadLaplaceSolnCells(fileName)
    dmnIDs = setDomainID(vtuMesh, cphiEP, cphiAB)

    print ("   Writing fibers and domains to VTK data structure")
    vtuMesh.GetPointData().AddArray(F)
    vtuMesh.GetPointData().AddArray(T)
    # vtuMesh.GetPointData().AddArray(S)

    vtuMesh.GetCellData().AddArray(dmnIDs)

    fileName = "fibers_domains.vtu"
    vtuWriter = vtk.vtkXMLUnstructuredGridWriter()
    vtuWriter.SetInputData(vtuMesh)
    vtuWriter.SetFileName(fileName)
    print ("   Writing to vtu file   --->   %s" % (fileName))
    vtuWriter.Write()
    dur = timedelta(seconds=time.perf_counter()-t1)
    print ("\n   Total time: ", dur)
    print ("========================================================")

#----------------------------------------------------------------------


