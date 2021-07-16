# Header
import sympy as sp
import numpy as np
import vtk
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

# Define symbols
x, y, z = sp.symbols('x y z')    # coordinate axes
t       = sp.symbols('t')        # time
omega   = sp.symbols('omega')    # angular frequency (displacement)
beta    = sp.symbols('beta')     # angular frequency (pressure)
rho     = sp.symbols('rho')      # density of the structure
mu      = sp.symbols('mu')       # shear modulus (Neo-Hookean model)

# Input parameters
nt      = 52                       # num. time steps
N       = 4
tmax    = 0.51
mshFile = '../../mesh/N%03d/mesh-complete.mesh.vtu' % (N)
fhdr    =  '../N%03d' % (N)

#=======================================================================
# Define displacement field
u1   = t * t * sp.sin(omega*y) * sp.sin(omega*z)
u2   = 0
u3   = 0
u    = sp.Matrix([u1, u2, u3])

print(bcolors.HEADER + "="*80 + bcolors.ENDC)
print("")
print(bcolors.OKBLUE + "Displacement field, u: " + bcolors.ENDC)
print("")
sp.pprint(u)
print("")

#=======================================================================
# Compute the deformation gradient tensor
F    = sp.eye(3)
F    = F + sp.Matrix([[sp.diff(u1, x), sp.diff(u1, y), sp.diff(u1, z)], \
                      [sp.diff(u2, x), sp.diff(u2, y), sp.diff(u2, z)], \
                      [sp.diff(u3, x), sp.diff(u3, y), sp.diff(u3, z)]])

# Jacobian
J    = F.det()

# Cauchy-Green deformation tensor and its inverse
C    = F.T * F
Ci   = C**-1
J23  = J**(-2/3)

#=======================================================================
# Now construct 1st & 2nd Piola-Kirchhoff Stress tensor
# Volumetric 2nd Piola-Kirchhoff stress tensor:

p    = t * t * sp.sin(beta*x) * sp.sin(beta*y) * sp.sin(beta*z)
print(bcolors.HEADER + "="*80 + bcolors.ENDC)
print("")
print(bcolors.OKBLUE + "Pressure, p: " + bcolors.ENDC)
print("")
sp.pprint(p)
print("")

# p defined above is the negative of `p' defined in Holzapfel book
Svol = (-p)*J*Ci

# Isochoric 2nd Piola-Kirhhoff stress tensor:
#    Neo-Hookean
Sb   = mu*sp.eye(3)

# TrS := Sb:C
TrS  = Sb[0,0]*C[0,0] + Sb[1,1]*C[1,1] + Sb[2,2]*C[2,2]

Siso = J23 * (Sb - (TrS/3)*Ci)
PK   = F*(Svol + Siso)

#=======================================================================
# Structure internal force, Div.(PK)
DivP = sp.diff(PK.col(0), x)
DivP = DivP.col_insert(1, sp.diff(PK.col(1),y))
DivP = DivP.col_insert(2, sp.diff(PK.col(2),z))
fs   = sp.Matrix([sum(DivP.row(0)), sum(DivP.row(1)), sum(DivP.row(2))]) 

#=======================================================================
# Inertial force, ddot(u)
fi   = sp.Matrix([sp.diff(u1, t, t), sp.diff(u2, t, t), sp.diff(u3, t, t)])

#=======================================================================
# Body force, fb
fb   = fi - (1/rho)*fs
fb_fn = sp.lambdify([t, omega, beta, rho, mu, x, y, z], fb, "numpy")

#=======================================================================
# Define time, domain and material parameters
# Set time parameters: t, T0
tr    = np.linspace(0, tmax, nt)

# Set parameters: omega, beta
omega = 0.1*np.pi
beta  = 0.2*np.pi

# Set material parameters:
# Density, g/cm3
rho   = 1.0

# Values chosen based on E = 3.0 dyn/cm2, nu = 0.5
# Isochoric: NeoHookean (mu, dyn/cm2)
mu    = 1.0

#=======================================================================
# Read VTU mesh file and load point coordinates
mshReader = vtk.vtkXMLUnstructuredGridReader()
mshReader.SetFileName(mshFile)
mshReader.Update()

msh = vtk.vtkUnstructuredGrid()
msh = mshReader.GetOutput()

msh_npts = msh.GetNumberOfPoints()
msh_pts  = np.zeros((msh_npts,3))
for ipt in range(0, msh_npts):
    msh_pts[ipt,:] = msh.GetPoint(ipt)

#=======================================================================
# Evaluate body force over the domain and write to file
bff = vtk.vtkDoubleArray()
bff.SetNumberOfComponents(3)
bff.Allocate(msh_npts)
bff.SetNumberOfTuples(msh_npts)
bff.SetName("FB")

print(bcolors.HEADER + "="*80 + bcolors.ENDC)
print("")
print(bcolors.OKBLUE + "Writing csv and vtu files from t=%.2f to t=%.2f" % \
    (tr[0], tr[nt-1]) + bcolors.ENDC)
print("")
i = 0
for time in tr:
    fbg = fb_fn(time, omega, beta, rho, mu, \
        msh_pts[:,0], msh_pts[:,1], msh_pts[:,2])
    i = i + 1
    for j in range(0, 3):
        fname = "%s/csv/fb%d_%02d.csv" % (fhdr, j+1, i)
        np.savetxt(fname, fbg[j,0,:], delimiter=',')
    for j in range(0, msh_npts):
        rtmp = fbg[:,0,j]
        bff.SetTuple3(j, rtmp[0], rtmp[1], rtmp[2])
    fname = "%s/vtu/fb_%02d.vtu" % (fhdr, i)
    msh.GetPointData().AddArray(bff)
    mshWrite = vtk.vtkXMLUnstructuredGridWriter()
    mshWrite.SetInputData(msh)
    mshWrite.SetFileName(fname)
    mshWrite.Write()
print(bcolors.HEADER + "="*80 + bcolors.ENDC)

#=======================================================================
# EOF

