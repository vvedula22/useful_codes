# Header
import numpy as np
import vtk

def getCV(fname1, fname2, dt):
    mshReader = vtk.vtkXMLUnstructuredGridReader()

    # Load dataset 1
    mshReader.SetFileName(fname1)
    mshReader.Update()
    msh = vtk.vtkUnstructuredGrid()
    msh = mshReader.GetOutput()
    msh_npts = msh.GetNumberOfPoints()
    msh_V1 = vtk.vtkDoubleArray()
    msh_V1 = msh.GetPointData().GetArray('EP_Action_potential')
    
    mshReader.SetFileName(fname2)
    mshReader.Update()
    msh_V2 = vtk.vtkDoubleArray()
    msh_V2 = msh.GetPointData().GetArray('EP_Action_potential')

    for ipt in range(0, msh_npts-1):
        v1 = msh_V1.GetTuple1(ipt)
        v2 = msh_V1.GetTuple1(ipt+1)
        if v1*v2<0:
            x1 = msh.GetPoint(ipt)
            x2 = msh.GetPoint(ipt+1)
            x0_1 = x1[1] + (x2[1]-x1[1])*v1/(v1-v2)
            break

    for ipt in range(0, msh_npts-1):
        v1 = msh_V2.GetTuple1(ipt)
        v2 = msh_V2.GetTuple1(ipt+1)
        if v1*v2<0:
            x1 = msh.GetPoint(ipt)
            x2 = msh.GetPoint(ipt+1)
            x0_2 = x1[1] + (x2[1]-x1[1])*v1/(v1-v2)
            break
    CV = (x0_2 - x0_1)/dt
    return CV
    

if __name__ == '__main__':
    h      = [0.10, 0.15, 0.20, 0.25, 0.3, \
              0.35, 0.40, 0.80, 1.00]
    dt     = 0.05
    ntime1 = 200
    ntime2 = 1800
    delta  = (ntime2-ntime1)*dt
    
    fname = "CV_dt%.2f.dat" %(dt)
    fid = open(fname,"a+")
    fid.write("Variables = h, CV\n")
    for i in range(0, len(h)):
        srcDir = "h%0.2f/1-procs_dt%.2f" %(h[i], dt)
        fname1 = "%s/result_%d.vtu" %(srcDir, ntime1)
        fname2 = "%s/result_%d.vtu" %(srcDir, ntime2)
        cv_cmps = 100.0*getCV(fname1, fname2, delta)
        fid.write("%.4f\t%.6f\n" %(h[i], cv_cmps))
    fid.close()

#=======================================================================
# EOF

