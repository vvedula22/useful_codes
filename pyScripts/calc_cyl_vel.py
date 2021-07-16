# Header
import numpy as np
import vtk

def calcAvgVel(FileName):
    mshReader = vtk.vtkXMLUnstructuredGridReader()
    mshReader.SetFileName(FileName)
    mshReader.Update()

    msh = vtk.vtkUnstructuredGrid()
    msh = mshReader.GetOutput()

    msh_npts = msh.GetNumberOfPoints()
    msh_vel  = vtk.vtkDoubleArray()
    msh_vel  = msh.GetPointData().GetArray('IB_Velocity')
    meanU = 0.0
    for ipt in range(0, msh_npts):
        v = msh_vel.GetTuple3(ipt)
        magV = np.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
        meanU = meanU + magV
    meanU = meanU / float(msh_npts)
    return meanU
    

if __name__ == '__main__':
    srcdir = "./24-procs/vtuFiles"
    nstart = 25
    nend   = 20000
    nfreq  = 25
    dt     = 0.0002

    nframe = int((nend - nstart)/nfreq) + 1
    for i in range(0, nframe):
        ntime = nstart + i*nfreq
        time  = float(ntime)*dt
        if (ntime < 100):
            fname = "%s/result_ib_%03d.vtu" %(srcdir, ntime)
        else:
            fname = "%s/result_ib_%d.vtu" %(srcdir, ntime)
        print("%.4f  %.6f" % (time, calcAvgVel(fname)))

#=======================================================================
# EOF

