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

# Define symbols
x, y, z = sp.symbols('x y z')    # coordinate axes
t       = sp.symbols('t')        # time
omega   = sp.symbols('omega')    # angular frequency (displacement)
beta    = sp.symbols('beta')     # angular frequency (pressure)
rho     = sp.symbols('rho')      # density of the structure
mu      = sp.symbols('mu')       # shear modulus (Neo-Hookean model)

# Input parameters
N       = 4
time    = 0.5
mshFile = '../../mesh/N%03d/mesh-complete.mesh.vtu' % (N)
fhdr    =  '../N%03d' % (N)

#=======================================================================
# Define displacement field
u1   = t * t * sp.sin(omega*y) * sp.sin(omega*z)
u2   = 0
u3   = 0

u  = sp.Matrix([u1])
u_fn = sp.lambdify([t, omega, x, y, z], u, "numpy")
print(bcolors.HEADER + "="*80 + bcolors.ENDC)
print("")
print(bcolors.OKBLUE + "Displacement field, u: " + bcolors.ENDC)
print("")
sp.pprint(u)
print("")

#=======================================================================
# Compute velocity field
v  = sp.Matrix([sp.diff(u1,t)])
v_fn = sp.lambdify([t, omega, x, y, z], v, "numpy")
print(bcolors.HEADER + "="*80 + bcolors.ENDC)
print("")
print(bcolors.OKBLUE + "Velocity field, v: " + bcolors.ENDC)
print("")
sp.pprint(v)
print("")

#=======================================================================
# Compute the deformation gradient tensor
F    = sp.eye(3)
F    = F + sp.Matrix([[sp.diff(u1, x), sp.diff(u1, y), sp.diff(u1, z)], \
                      [sp.diff(u2, x), sp.diff(u2, y), sp.diff(u2, z)], \
                      [sp.diff(u3, x), sp.diff(u3, y), sp.diff(u3, z)]])
#F_fn = sp.lambdify([t, omega, x, y, z], F, "numpy")
print(bcolors.HEADER + "="*80 + bcolors.ENDC)
print("")
print(bcolors.OKBLUE + "Def. gradient, F: " + bcolors.ENDC)
print("")
sp.pprint(F)
print("")

F12  = F[0,1]
F13  = F[0,2]
F12_fn = sp.lambdify([t, omega, x, y, z], F12, "numpy")
F13_fn = sp.lambdify([t, omega, x, y, z], F13, "numpy")

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
p_fn = sp.lambdify([t, beta, x, y, z], p, "numpy")
print(bcolors.HEADER + "="*80 + bcolors.ENDC)
print("")
print(bcolors.OKBLUE + "Pressure, p: " + bcolors.ENDC)
print("")
sp.pprint(p)
print("")

# Isochoric 2nd Piola-Kirhhoff stress tensor:
#    Neo-Hookean
Sb   = mu*sp.eye(3)

# TrS := Sb:C
TrS  = Sb[0,0]*C[0,0] + Sb[1,1]*C[1,1] + Sb[2,2]*C[2,2]

Siso = J23 * (Sb - (TrS/3)*Ci)
PK   = F*Siso

sigm = PK * F.T
sigm = sigm / J
s_fn = sp.lambdify([t, omega, mu, x, y, z], sigm, "numpy")

#=======================================================================
# Define time, domain and material parameters
# Set time parameters: t
ti    = time

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
# Read vtk mesh
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
# Evaluate data at the mesh point coordiantes
ug = u_fn(ti, omega, msh_pts[:,0], msh_pts[:,1], msh_pts[:,2])
vg = v_fn(ti, omega, msh_pts[:,0], msh_pts[:,1], msh_pts[:,2])
pg = p_fn(ti, beta, msh_pts[:,0], msh_pts[:,1], msh_pts[:,2])
sg = s_fn(ti, omega, mu, msh_pts[:,0], msh_pts[:,1], msh_pts[:,2])

F12g = F12_fn(ti, omega, msh_pts[:,0], msh_pts[:,1], msh_pts[:,2])
F13g = F13_fn(ti, omega, msh_pts[:,0], msh_pts[:,1], msh_pts[:,2])

#=======================================================================
# Copy data to vtk objects and write to file
vtkU = vtk.vtkDoubleArray()
vtkU.SetNumberOfComponents(3)
vtkU.Allocate(msh_npts)
vtkU.SetNumberOfTuples(msh_npts)
vtkU.SetName("ST_Displacement")

vtkV = vtk.vtkDoubleArray()
vtkV.SetNumberOfComponents(3)
vtkV.Allocate(msh_npts)
vtkV.SetNumberOfTuples(msh_npts)
vtkV.SetName("ST_Velocity")

vtkP = vtk.vtkDoubleArray()
vtkP.SetNumberOfComponents(1)
vtkP.Allocate(msh_npts)
vtkP.SetNumberOfTuples(msh_npts)
vtkP.SetName("ST_Pressure")

vtkS = vtk.vtkDoubleArray()
vtkS.SetNumberOfComponents(6)
vtkS.Allocate(msh_npts)
vtkS.SetNumberOfTuples(msh_npts)
vtkS.SetName("ST_Stress")

vtkF = vtk.vtkDoubleArray()
vtkF.SetNumberOfComponents(9)
vtkF.Allocate(msh_npts)
vtkF.SetNumberOfTuples(msh_npts)
vtkF.SetName("ST_Deformation_gradient")

print(bcolors.HEADER + "="*80 + bcolors.ENDC)
print("")
print(bcolors.OKBLUE + "Writing vtu files at t=%.4f" % \
    (time) + bcolors.ENDC)
print("")

sg11 = sg[0,0]
sg22 = sg[1,1]
sg33 = sg[2,2]
sg12 = sg[0,1]
sg13 = sg[0,2]

for j in range(0, msh_npts):
    vtkP.SetTuple1(j, pg[j])
    vtkU.SetTuple3(j, ug[0,0,j], 0.0, 0.0)
    vtkV.SetTuple3(j, vg[0,0,j], 0.0, 0.0)
    vtkS.SetTuple6(j, sg11[j], sg22[j], sg33[j], \
                      sg12[j], sg13[j], 0.0)
    vtkF.SetTuple9(j, 1.0, F12g[j], F13g[j], \
                      0.0, 1.0, 0.0,\
                      0.0, 0.0, 1.0)

fname = "%s/exact_soln_t%.4f.vtu" % (fhdr, ti)
msh.GetPointData().AddArray(vtkP)
msh.GetPointData().AddArray(vtkU)
msh.GetPointData().AddArray(vtkV)
msh.GetPointData().AddArray(vtkS)
msh.GetPointData().AddArray(vtkF)
mshWrite = vtk.vtkXMLUnstructuredGridWriter()
mshWrite.SetInputData(msh)
mshWrite.SetFileName(fname)
mshWrite.Write()

print(bcolors.HEADER + "="*80 + bcolors.ENDC)

#=======================================================================
# EOF

