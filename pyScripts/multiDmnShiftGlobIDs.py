import os
import sys
import vtk
import glob
import shutil

# User defined inputs
MESH_COMPLETE_DOMAIN_1 = "pipe_demo-mesh-complete_domain-1"
MESH_COMPLETE_DOMAIN_2 = "pipe_demo-mesh-complete_domain-2"

NODEIDSTR = 'GlobalNodeID'
CELLIDSTR = 'GlobalElementID'

#----------------------------------------------------------------------
def loadVTU(fileName):

    print "   Loading vtu file   <---   %s" % (fileName)
    vtuReader = vtk.vtkXMLUnstructuredGridReader()
    vtuReader.SetFileName(fileName)
    vtuReader.Update()

    vtuMesh = vtk.vtkUnstructuredGrid()
    vtuMesh = vtuReader.GetOutput()

    return vtuMesh
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def loadVTP(fileName):

    print "   Loading vtp file   <---   %s" % (fileName)
    vtpReader = vtk.vtkXMLPolyDataReader()
    vtpReader.SetFileName(fileName)
    vtpReader.Update()

    vtpPoly = vtk.vtkPolyData()
    vtpPoly = vtpReader.GetOutput()

    return vtpPoly
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
def shiftMeshNodeIDs(vtuMesh, iShift):

    numPoints = vtuMesh.GetNumberOfPoints()
    globalNodeID = vtuMesh.GetPointData().GetArray(NODEIDSTR)
    for ipt in xrange(0, numPoints):
        tempID = globalNodeID.GetTuple1(ipt)
        tempID = tempID - iShift
        globalNodeID.SetTuple1(ipt, tempID)

    vtuMesh.GetPointData().RemoveArray(NODEIDSTR)
    vtuMesh.GetPointData().AddArray(globalNodeID)

    return vtuMesh
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def shiftMeshCellIDs(vtuMesh, iShift):

    numCells = vtuMesh.GetNumberOfCells()
    globalElemID = vtuMesh.GetCellData().GetArray(CELLIDSTR)
    for icel in xrange(0, numCells):
        tempID = globalElemID.GetTuple1(icel)
        tempID = tempID - iShift
        globalElemID.SetTuple1(icel, tempID)

    vtuMesh.GetCellData().RemoveArray(CELLIDSTR)
    vtuMesh.GetCellData().AddArray(globalElemID)

    return vtuMesh
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# Main function
if __name__ == '__main__':

    print "========================================================"
    path2Dmn1 = MESH_COMPLETE_DOMAIN_1
    path2Dmn2 = MESH_COMPLETE_DOMAIN_2

    # load mesh-complete vtu file for domain 1
    mshFile1 = "%s/mesh-complete.mesh.vtu" % (path2Dmn1)
    if not os.path.isfile(mshFile1):
        print("Error: failed to find mesh-complete.mesh.vtu")
        sys.exit()

    vtuMesh1 = loadVTU(mshFile1)
    numPoints1 = vtuMesh1.GetNumberOfPoints()
    numCells1 = vtuMesh1.GetNumberOfCells()

    print("   Total nodes in domain 1: %d" % (numPoints1))
    print("   Total cells in domain 1: %d" % (numCells1))

    # load mesh-complete vtu file for domain 2
    mshFile2 = "%s/mesh-complete.mesh.vtu" %(path2Dmn2)
    if not os.path.isfile(mshFile2):
        print("Error: failed to find mesh-complete.mesh.vtu")
        sys.exit()

    vtuMesh2 = loadVTU(mshFile2)
    numPoints2 = vtuMesh2.GetNumberOfPoints()
    numCells2 = vtuMesh2.GetNumberOfCells()

    print("   Total nodes in domain 2: %d" % (numPoints2))
    print("   Total cells in domain 2: %d" % (numCells2))

    # Shift global nodes and elements for domain mesh 2
    vtuMesh2 = shiftMeshNodeIDs(vtuMesh2, numPoints1)
    vtuMesh2 = shiftMeshCellIDs(vtuMesh2, numCells1)

    # Create a new mesh-complete folder for domain 2
    path2Dmn2_new = "%s_new" % os.path.abspath(path2Dmn2)
    if os.path.isdir(path2Dmn2_new):
        shutil.rmtree(path2Dmn2_new)
    os.mkdir(path2Dmn2_new)
    os.mkdir("%s/mesh-surfaces" % (path2Dmn2_new))

    # write a new mesh-complete vtu file for domain 2
    fileName = "%s/mesh-complete.mesh.vtu" % (path2Dmn2_new)
    vtuWriter = vtk.vtkXMLUnstructuredGridWriter()
    vtuWriter.SetInputData(vtuMesh2)
    vtuWriter.SetFileName(fileName)
    print "   Writing to vtu file   --->   %s" % (fileName)
    vtuWriter.Write()

    # process wall and exterior vtp files, if present
    os.chdir(path2Dmn2)
    faceFile = "%s/mesh-complete.exterior.vtp" %(path2Dmn2)
    for faceFile in glob.glob("*.vtp"):
        vtpPoly = loadVTP(faceFile)
        vtpPoly = shiftMeshNodeIDs(vtpPoly, numPoints1)
        vtpPoly = shiftMeshCellIDs(vtpPoly, numPoints1)
        fileName = "%s/%s" % (path2Dmn2_new, faceFile)
        vtpWriter = vtk.vtkXMLPolyDataWriter()
        vtpWriter.SetInputData(vtpPoly)
        vtpWriter.SetFileName(fileName)
        print "   Writing to vtp file   --->   %s" % (fileName)
        vtpWriter.Write()

    # now process all the faces file in mesh-surfaces
    if os.path.isdir("mesh-surfaces"):
        os.chdir("mesh-surfaces")
        for faceFile in glob.glob("*.vtp"):
            vtpPoly = loadVTP(faceFile)
            vtpPoly = shiftMeshNodeIDs(vtpPoly, numPoints1)
            vtpPoly = shiftMeshCellIDs(vtpPoly, numPoints1)
            fileName = "%s/mesh-surfaces/%s" % \
               (path2Dmn2_new, faceFile)
            vtpWriter = vtk.vtkXMLPolyDataWriter()
            vtpWriter.SetInputData(vtpPoly)
            vtpWriter.SetFileName(fileName)
            print "   Writing to vtp file   --->   %s" % (fileName)
            vtpWriter.Write()

    print "========================================================"

#----------------------------------------------------------------------

