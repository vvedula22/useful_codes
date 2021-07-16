# Header
import sympy as sp
import numpy as np
import vtk
import os

sp.init_printing(use_unicode = True)

class bcolors:
    HEADER  = '\033[95m'
    OKBLUE  = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL    = '\033[91m'
    ENDC    = '\033[0m'
    BOLD    = '\033[1m'
    ULINE   = '\033[4m'

def genFaceBC(fileName, fdir, nv) :
    x, y, z  = sp.symbols('x y z')    # coordinate axes
    t        = sp.symbols('t')        # time
    omega    = sp.symbols('omega')    # angular frequency (displacement)
    beta     = sp.symbols('beta')     # angular frequency (pressure)
    rho      = sp.symbols('rho')      # density of the structure
    mu       = sp.symbols('mu')       # shear modulus (Neo-Hookean model)

    u1   = t * t * sp.sin(omega*y) * sp.sin(omega*z)
    u2   = 0
    u3   = 0

    F    = sp.eye(3)
    F    = F + sp.Matrix([[sp.diff(u1, x), sp.diff(u1, y), sp.diff(u1, z)], \
                      [sp.diff(u2, x), sp.diff(u2, y), sp.diff(u2, z)], \
                      [sp.diff(u3, x), sp.diff(u3, y), sp.diff(u3, z)]])

    J    = F.det()

    C    = F.T * F
    Ci   = C**-1
    J23  = J**(-2/3)

    p    = t * t * sp.sin(beta*x) * sp.sin(beta*y) * sp.sin(beta*z)

    # p defined above is the negative of `p' defined in Holzapfel book
    Svol = (-p)*J*Ci

    Sb   = mu*sp.eye(3)

    TrS  = Sb[0,0]*C[0,0] + Sb[1,1]*C[1,1] + Sb[2,2]*C[2,2]

    Siso = J23 * (Sb - (TrS/3)*Ci)
    PK   = F*(Svol + Siso)

    hv = PK * nv
    # print(bcolors.HEADER + "="*80 + bcolors.ENDC)
    # print("")
    # print(bcolors.OKBLUE + "traction, h:" + bcolors.ENDC)
    # sp.pprint(nv)
    # print("")
    # sp.pprint(sp.simplify(hv))
    # print("")

    h_fn = sp.lambdify([t, omega, beta, mu, x, y, z], hv, "numpy")

    nt    = 52
    tr    = np.linspace(0, 0.51, nt)

    omega = 0.1*np.pi
    beta  = 0.2*np.pi

    mu    = 1.0

    vtpReader = vtk.vtkXMLPolyDataReader()
    vtpReader.SetFileName(fileName)
    vtpReader.Update()

    pdata = vtk.vtkPolyData()
    pdata = vtpReader.GetOutput()

    pdata_npts = pdata.GetNumberOfPoints()
    pdata_pts  = np.zeros((pdata_npts,3))
    pdata_nid  = np.zeros((pdata_npts,1),int)
    pdata_nid  = pdata.GetPointData().GetArray('GlobalNodeID')
    for ipt in range(0, pdata_npts):
        pdata_pts[ipt,:] = pdata.GetPoint(ipt)

    path, fname = os.path.split(fileName)
    fhdr, ext   = os.path.splitext(fname)

    fname = "%s/csv/bc_%s_nodeid.csv" % (fdir, fhdr)
    np.savetxt(fname, pdata_nid, delimiter=',', fmt='%d')

    print(bcolors.HEADER + "="*80 + bcolors.ENDC)
    print("")
    print(bcolors.OKBLUE + "Writing BC data for face %s from t=%.4f to t=%.4f" % \
        (fhdr, tr[0], tr[nt-1]) + bcolors.ENDC)
    print("")
    i = 0
    if (fhdr == 'X0') | (fhdr == 'X1') :
        for time in tr :
            hg = h_fn(time, omega, beta, mu, pdata_pts[:,0], pdata_pts[:,1], pdata_pts[:,2])
            i  = i + 1
            for j in range(0, np.size(hg,0)) :
                fname = "%s/csv/bc_%s_h%d_%02d.csv" % (fdir, fhdr, j+1, i)
                np.savetxt(fname, hg[j,0,:])
    elif (fhdr == 'Y0') | (fhdr == 'Y1') :
        for time in tr :
            hg = h_fn(time, omega, beta, mu, pdata_pts[:,0], pdata_pts[:,1], pdata_pts[:,2])
            i  = i + 1
            fname = "%s/csv/bc_%s_h1_%02d.csv" % (fdir, fhdr, i)
            np.savetxt(fname, hg[0,0])
            fname = "%s/csv/bc_%s_h2_%02d.csv" % (fdir, fhdr, i)
            np.savetxt(fname, hg[1,0])
    elif (fhdr == 'Z0') | (fhdr == 'Z1') :
        for time in tr :
            hg = h_fn(time, omega, beta, mu, pdata_pts[:,0], pdata_pts[:,1], pdata_pts[:,2])
            i  = i + 1
            fname = "%s/csv/bc_%s_h1_%02d.csv" % (fdir, fhdr, i)
            np.savetxt(fname, hg[0,0])
            fname = "%s/csv/bc_%s_h3_%02d.csv" % (fdir, fhdr, i)
            np.savetxt(fname, hg[2,0])



if __name__ == '__main__':
    N       = 4
    fhdr    =  '../N%03d' % (N)

    nv = sp.Matrix([-1, 0, 0])
    fname = '../../mesh/N%03d/mesh-surfaces/X0.vtp' %(N)
    genFaceBC(fname, fhdr, nv)

    nv = sp.Matrix([ 1, 0, 0])
    fname = '../../mesh/N%03d/mesh-surfaces/X1.vtp' %(N)
    genFaceBC(fname, fhdr, nv)

    nv = sp.Matrix([ 0,-1, 0])
    fname = '../../mesh/N%03d/mesh-surfaces/Y0.vtp' %(N)
    genFaceBC(fname, fhdr, nv)

    nv = sp.Matrix([ 0, 1, 0])
    fname = '../../mesh/N%03d/mesh-surfaces/Y1.vtp' %(N)
    genFaceBC(fname, fhdr, nv)

    nv = sp.Matrix([ 0, 0,-1])
    fname = '../../mesh/N%03d/mesh-surfaces/Z0.vtp' %(N)
    genFaceBC(fname, fhdr, nv)

    nv = sp.Matrix([ 0, 0, 1])
    fname = '../../mesh/N%03d/mesh-surfaces/Z1.vtp' %(N)
    genFaceBC(fname, fhdr, nv)