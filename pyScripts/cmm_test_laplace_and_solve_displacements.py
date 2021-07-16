import os
import sys
import numpy as np
import vtk
import time
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve, cg
from scipy import sparse

MESH_COMPLETE_PATH = '../short_cyl_lumen-mesh-complete/'
INFLATION_PRESSURE = 1.081279E5
DEFORMABLE_E = 10000000.0
DEFORMABLE_T = 0.016
POISSON_RATIO = 0.49
KCONS = 5.0/6.0
PS_ITER = 100

# ==============================================================================

def computeDisplacementNorm(disp_sol):

  norm_val = 0.0
  
  for i in xrange(0, len(disp_sol)):
    temp_val = np.sqrt(disp_sol[i][0]*disp_sol[i][0] + \
                       disp_sol[i][1]*disp_sol[i][1] + \
                       disp_sol[i][2]*disp_sol[i][2])
    norm_val = norm_val + temp_val
  
  return norm_val

# ==============================================================================

def loadMeshAndWalls():

  vtu_reader = vtk.vtkXMLUnstructuredGridReader()
  vtu_reader.SetFileName(MESH_COMPLETE_PATH + 'mesh-complete.mesh.vtu')
  vtu_reader.Update()
  vtu_mesh = vtk.vtkUnstructuredGrid()
  vtu_mesh = vtu_reader.GetOutput()
  
  vtp_reader = vtk.vtkXMLPolyDataReader()
  vtp_reader.SetFileName(MESH_COMPLETE_PATH + 'walls_combined.vtp')
  vtp_reader.Update()
  walls_combined_vtp = vtk.vtkPolyData()
  walls_combined_vtp = vtp_reader.GetOutput()
  
  return vtu_mesh, walls_combined_vtp

# ==============================================================================

def assignDirichletBC(mesh_vtu, bc_face_vtp, value, dirichlet_BC):

  numMeshPoints = mesh_vtu.GetNumberOfPoints()
  wallID = mesh_vtu.GetPointData().GetArray('GlobalNodeID')
  numFacePoints = bc_face_vtp.GetNumberOfPoints()
  faceID = bc_face_vtp.GetPointData().GetArray('GlobalNodeID')
  
  num_matches = 0
  for i_face in xrange(0, numFacePoints):
    iid = faceID.GetTuple1(i_face)
    for i_wall in xrange(0, numMeshPoints):
      iid_wall = wallID.GetTuple1(i_wall)
      if(iid == iid_wall):
        dirichlet_BC[i_wall] = value
        num_matches = num_matches + 1
        break

# ==============================================================================

def laplace_stiffness_local_tet(p0, p1, p2, p3):

  # First, compute the derivatives of the spatial coordinates wrt. parent coordinates
  # Using the shape functions in the parent element. For this example, we use
  # the following shape functions:
  #  1. N0 = r
  #  2. N1 = s
  #  3. N2 = t
  #  4. N3 = 1 - r - s - t
  dxdr = p0[0] - p3[0]
  dxds = p1[0] - p3[0]
  dxdt = p2[0] - p3[0]
  
  dydr = p0[1] - p3[1]
  dyds = p1[1] - p3[1]
  dydt = p2[1] - p3[1]
  
  dzdr = p0[2] - p3[2]
  dzds = p1[2] - p3[2]
  dzdt = p2[2] - p3[2]
  
  # We then put these into a 3x3 jacobian matrix and take the inverse to solve for the
  # inverse jacobians that we need to construct the element stiffness matrix.
  # The resulting expressions for the derivatives are as follows:
  drdx = dzdt*dyds - dzds*dydt
  dsdx = dzdr*dydt - dzdt*dydr
  dtdx = dzds*dydr - dyds*dzdr
  
  drdy = dxdt*dzds - dzdt*dxds
  dsdy = dzdt*dxdr - dzdr*dxdt
  dtdy = dxds*dzdr - dzds*dxdr
  
  drdz = dydt*dxds - dxdt*dyds
  dsdz = dxdt*dydr - dydt*dxdr
  dtdz = dyds*dxdr - dxds*dydr
  
  # Solve for the determinant of the 3x3 jacobian matrix as part of the inverse
  detj = dxdr*drdx + dxds*dsdx + dxdt*dtdx
  
  # Manually input the entries of the element stiffness matrix
  kab = np.zeros((4, 4))
  
  # We divide each entry by 6.0 because this is the volume of the parent
  # tetrahedra and the weights of the Gauss quadrature points have to sum
  # up to the volume of the parent element. We also divide by the determinant
  # of the jacobian. We get 2 determinants in the denominator from the change
  # variables of the derivatives, and 1 determinant in the numerator from
  # the integral change of variables.
  if(abs(detj) > 0.00000001):
  
    dN1dx = drdx
    dN1dy = drdy
    dN1dz = drdz
    
    dN2dx = dsdx
    dN2dy = dsdy
    dN2dz = dsdz
    
    dN3dx = dtdx
    dN3dy = dtdy
    dN3dz = dtdz
    
    dN4dx = -drdx - dsdx - dtdx
    dN4dy = -drdy - dsdy - dtdy
    dN4dz = -drdz - dsdz - dtdz
  
    kab[0][0] = (dN1dx*dN1dx + dN1dy*dN1dy + dN1dz*dN1dz)/(6.0 * detj)
    kab[0][1] = (dN1dx*dN2dx + dN1dy*dN2dy + dN1dz*dN2dz)/(6.0 * detj)
    kab[0][2] = (dN1dx*dN3dx + dN1dy*dN3dy + dN1dz*dN3dz)/(6.0 * detj)
    kab[0][3] = (dN1dx*dN4dx + dN1dy*dN4dy + dN1dz*dN4dz)/(6.0 * detj)
    
    kab[1][0] = kab[0][1]
    kab[1][1] = (dN2dx*dN2dx + dN2dy*dN2dy + dN2dz*dN2dz)/(6.0 * detj)
    kab[1][2] = (dN2dx*dN3dx + dN2dy*dN3dy + dN2dz*dN3dz)/(6.0 * detj)
    kab[1][3] = (dN2dx*dN4dx + dN2dy*dN4dy + dN2dz*dN4dz)/(6.0 * detj)
    
    kab[2][0] = kab[0][2]
    kab[2][1] = kab[1][2]
    kab[2][2] = (dN3dx*dN3dx + dN3dy*dN3dy + dN3dz*dN3dz)/(6.0 * detj)
    kab[2][3] = (dN3dx*dN4dx + dN3dy*dN4dy + dN3dz*dN4dz)/(6.0 * detj)
    
    kab[3][0] = kab[0][3]
    kab[3][1] = kab[1][3]
    kab[3][2] = kab[2][3]
    kab[3][3] = (dN4dx*dN4dx + dN4dy*dN4dy + dN4dz*dN4dz)/(6.0 * detj)
  else:
    tet_vol = vtk.vtkTetra().ComputeVolume(p0, p1, p2, p3)
    print('Singular jacobian detected with volume: ' + str(tet_vol) + ' and detj: ' + str(detj))
  
  return kab

# ==============================================================================

def cmm_stiffness_local(p0, p1, p2, local_stresses):

  # Create containers for the element stiffness and element forcing vector
  kab = np.zeros((9, 9))
  fa = np.zeros(9)
  
  # Create the rotation matrix by getting a set of coordinates in the element
  # i.e., two coordinate directions will be in the plane of the element and
  # one coordinate will be perpendicularly out
  e0_local = (p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2])
  temp_mag = np.sqrt(e0_local[0]*e0_local[0] + e0_local[1]*e0_local[1] \
                     + e0_local[2]*e0_local[2])
  e0_local = (e0_local[0]/temp_mag, e0_local[1]/temp_mag, e0_local[2]/temp_mag)
  
  etemp_local = (p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2])
  temp_mag = np.sqrt(etemp_local[0]*etemp_local[0] + \
                     etemp_local[1]*etemp_local[1] + \
                     etemp_local[2]*etemp_local[2])
  etemp_local = (etemp_local[0]/temp_mag, etemp_local[1]/temp_mag, \
                 etemp_local[2]/temp_mag)
                 
  e2_local = np.cross(e0_local, etemp_local)
  temp_mag = np.sqrt(e2_local[0]*e2_local[0] + e2_local[1]*e2_local[1] \
                     + e2_local[2]*e2_local[2])
  e2_local = (e2_local[0]/temp_mag, e2_local[1]/temp_mag, e2_local[2]/temp_mag)
  
  e1_local = np.cross(e2_local, e0_local)
  temp_mag = np.sqrt(e1_local[0]*e1_local[0] + e1_local[1]*e1_local[1] \
                     + e1_local[2]*e1_local[2])
  e1_local = (e1_local[0]/temp_mag, e1_local[1]/temp_mag, e1_local[2]/temp_mag)
  
  # Making sure all the local bases are orthogonal
  assert(abs(np.dot(e0_local, e1_local)) < 1e-6)
  assert(abs(np.dot(e0_local, e2_local)) < 1e-6)
  assert(abs(np.dot(e1_local, e2_local)) < 1e-6)
  
  Theta = np.zeros((3, 3))
  Theta[0][0] = e0_local[0]
  Theta[0][1] = e0_local[1]
  Theta[0][2] = e0_local[2]
  Theta[1][0] = e1_local[0]
  Theta[1][1] = e1_local[1]
  Theta[1][2] = e1_local[2]
  Theta[2][0] = e2_local[0]
  Theta[2][1] = e2_local[1]
  Theta[2][2] = e2_local[2]
  
  # Compute the coordinates of the three vertices in the local system using the
  # rotation matrix Theta
  p0_local = np.dot(Theta, p0)
  p1_local = np.dot(Theta, p1)
  p2_local = np.dot(Theta, p2)
  
  # Compute the determinant of the Jacobian matrix and the area of the elements.
  # Here, we have taken the following to be the shape functions:
  #  1. N0 = r
  #  2. N1 = s
  #  3. N2 = 1 - r - s
  detJ = (p0_local[0]-p2_local[0])*(p1_local[1]-p2_local[1]) - \
         (p0_local[1]-p2_local[1])*(p1_local[0]-p2_local[0])
  area = 0.5*detJ
  
  # These are the derivatives of the parent coordinates wrt local coordinates
  dr_dxl = (1.0/detJ)*(p1_local[1]-p2_local[1])
  ds_dxl = (1.0/detJ)*(p2_local[1]-p0_local[1])
  dr_dyl = (1.0/detJ)*(p2_local[0]-p1_local[0])
  ds_dyl = (1.0/detJ)*(p0_local[0]-p2_local[0])
  
  # Now, write the shape function derivatives in terms of these
  dN0_dxl = dr_dxl
  dN0_dyl = dr_dyl
  dN1_dxl = ds_dxl
  dN1_dyl = ds_dyl
  dN2_dxl = -dr_dxl - ds_dxl
  dN2_dyl = -dr_dyl - ds_dyl
  
  # Now, we assume a strain vector like the following:
  # e_l = (du_dxl, dv_dyl, du_dyl+dv_dxl, dw_dxl, dw_dyl)
  # We wish to compute the B matrix that relates the strain vector to the
  # displacement vector u_l = (u0_l, v0_l, w0_l, u1_l, v1_l, w1_l, u2_l, v2_l, w2_l)
  # i.e. => e_l = B*u_l
  B_mat = np.zeros((5, 9))
  B_mat[0][0] = dN0_dxl
  B_mat[0][3] = dN1_dxl
  B_mat[0][6] = dN2_dxl
  B_mat[1][1] = dN0_dyl
  B_mat[1][4] = dN1_dyl
  B_mat[1][7] = dN2_dyl
  B_mat[2][0] = dN0_dyl
  B_mat[2][1] = dN0_dxl
  B_mat[2][3] = dN1_dyl
  B_mat[2][4] = dN1_dxl
  B_mat[2][6] = dN2_dyl
  B_mat[2][7] = dN2_dxl
  B_mat[3][2] = dN0_dxl
  B_mat[3][5] = dN1_dxl
  B_mat[3][8] = dN2_dxl
  B_mat[4][2] = dN0_dyl
  B_mat[4][5] = dN1_dyl
  B_mat[4][8] = dN2_dyl
  
  # Form the D matrix of material coefficients
  D_mat = np.zeros((5, 5))
  D_coef = DEFORMABLE_E/(1.0 - POISSON_RATIO*POISSON_RATIO)
  D_mat[0][0] = D_coef
  D_mat[0][1] = POISSON_RATIO*D_coef
  D_mat[1][0] = POISSON_RATIO*D_coef
  D_mat[1][1] = D_coef
  D_mat[2][2] = 0.5*D_coef*(1.0 - POISSON_RATIO)
  D_mat[3][3] = 0.5*D_coef*KCONS*(1.0 - POISSON_RATIO)
  D_mat[4][4] = 0.5*D_coef*KCONS*(1.0 - POISSON_RATIO)
  
  # Compute the element stiffness matrix as 0.5 * (thickness)* B^T * D * B * detJ
  # The 0.5 factor comes from the Gauss quadrature rules which says that the
  # Gauss weights for the parent domain triangle must sum up to the area
  # of the triangle, which is 0.5. The detJ comes from the change in coordinates
  # for the stiffness integral
  B_transpose = np.transpose(B_mat)
  k_local = np.dot(B_transpose, D_mat)
  k_local = np.dot(k_local, B_mat)
  k_coef = area*DEFORMABLE_T
  for i in xrange(0, 9):
    for j in xrange(0, 9):
      k_local[i][j] = k_coef*k_local[i][j]
  
  # Now we need to rotate this k_local back into the global coordinate frame
  # Form the 9x9 rotation tensor
  rot_mat = np.zeros((9,9))
  for i in xrange(0, 3):
    for j in xrange(0, 3):
      rot_mat[i][j] = Theta[i][j]
      rot_mat[i+3][j+3] = Theta[i][j]
      rot_mat[i+6][j+6] = Theta[i][j]
  
  rot_mat_t = np.transpose(rot_mat)
  Theta_t = np.transpose(Theta)
      
  k_local = np.dot(rot_mat_t, k_local)
  kab = np.dot(k_local, rot_mat)
  
  # Add in prestress terms
  assert(len(local_stresses) == 6)
  
  # Form the 3x3 Cauchy stress tensor associated with the pre-stress from
  # the local_stresses input vector
  # These stress components are: {o_xx, o_yy, o_zz, tau_xy, tau_xz, tau_yz}
  pre_stress_global = np.zeros((3,3))
  pre_stress_global[0][0] = local_stresses[0]
  pre_stress_global[1][1] = local_stresses[1]
  pre_stress_global[2][2] = local_stresses[2]
  pre_stress_global[0][1] = local_stresses[3]
  pre_stress_global[1][0] = local_stresses[3]
  pre_stress_global[0][2] = local_stresses[4]
  pre_stress_global[2][0] = local_stresses[4]
  pre_stress_global[1][2] = local_stresses[5]
  pre_stress_global[2][1] = local_stresses[5]
  
  # Rotate this global pre-stress tensor into the local coordinates. Not sure if the order of these operations is correct?
  pre_stress_local = np.dot(Theta, pre_stress_global)
  pre_stress_local = np.dot(pre_stress_local, Theta_t)
  
  # Assert that this local pre-stress tensor is symmetric still
  assert(abs(pre_stress_local[0][1] - pre_stress_local[1][0]) < 1e-3)
  assert(abs(pre_stress_local[0][2] - pre_stress_local[2][0]) < 1e-3)
  assert(abs(pre_stress_local[2][1] - pre_stress_local[1][2]) < 1e-3)
  
  # Apply the components of this local pre-stress vector to a local residual
  local_res = np.zeros(9)
  local_res[0] = dN0_dxl*pre_stress_local[0][0] + dN0_dyl*pre_stress_local[0][1]
  local_res[1] = dN0_dxl*pre_stress_local[0][1] + dN0_dyl*pre_stress_local[1][1]
  local_res[2] = dN0_dxl*pre_stress_local[0][2] + dN0_dyl*pre_stress_local[2][1]
  local_res[3] = dN1_dxl*pre_stress_local[0][0] + dN1_dyl*pre_stress_local[0][1]
  local_res[4] = dN1_dxl*pre_stress_local[0][1] + dN1_dyl*pre_stress_local[1][1]
  local_res[5] = dN1_dxl*pre_stress_local[0][2] + dN1_dyl*pre_stress_local[2][1]
  local_res[6] = dN2_dxl*pre_stress_local[0][0] + dN2_dyl*pre_stress_local[0][1]
  local_res[7] = dN2_dxl*pre_stress_local[0][1] + dN2_dyl*pre_stress_local[1][1]
  local_res[8] = dN2_dxl*pre_stress_local[0][2] + dN2_dyl*pre_stress_local[2][1]
  
  # Rotate these local residuals back to the global coordinate system
  global_res = np.dot(rot_mat_t, local_res)
  
  # Add contributions from these residuals to the element forcing vector, fa
  for i in xrange(0, 9):
    fa[i] = fa[i] - area*DEFORMABLE_T*global_res[i]
  
  # Now form the right hand side forcing vector as a normal force of the deformable
  # pressure
  for ip in xrange(0, 3):
    for jp in xrange(0, 3):
      fa[3*ip+jp] = fa[3*ip+jp] + area*INFLATION_PRESSURE*e2_local[jp]/3.0
  
  return kab, fa

# ==============================================================================

def performLaplaceSolve(mesh_vtu, dirichlet_BC):

  numPoints = mesh_vtu.GetNumberOfPoints()
  numCells = mesh_vtu.GetNumberOfCells()

  # Define an ID array to store the global equation number of every unknown
  ID = np.zeros(numMeshPoints)
  num_unknowns = 0
  for i_wall in xrange(0, numMeshPoints):
    if(dirichlet_BC[i_wall] < 0.0):
      ID[i_wall] = num_unknowns
      num_unknowns = num_unknowns + 1
    else:
      ID[i_wall] = -1.0
      
  print('Number of mesh unknowns: ' + str(num_unknowns))
  
  # Define the LM array to store the global ID number for all element nodes
  LM = np.zeros((numCells, 4))
  for i_cell in xrange(0, numCells):
    temp_cell = mesh_vtu.GetCell(i_cell)
    pts_cell = temp_cell.GetPointIds()
    for i_tet in xrange(0, 4):
      iid = pts_cell.GetId(i_tet)
      LM[i_cell][i_tet] = ID[iid]
      
  # Loop over all the elements and add their element stiffness matrix into the
  # global stiffness matrix
  globalK_i = []
  globalK_j = []
  globalK_v = []
  globalF = np.zeros(num_unknowns)
  for i_cell in xrange(0, numCells):
    
    temp_cell = mesh_vtu.GetCell(i_cell)
    pts_cell = temp_cell.GetPointIds()
    # Get the nodal coordinates
    p0 = temp_cell.GetPoints().GetPoint(0)
    p1 = temp_cell.GetPoints().GetPoint(1)
    p2 = temp_cell.GetPoints().GetPoint(2)
    p3 = temp_cell.GetPoints().GetPoint(3)
    
    kab = laplace_stiffness_local_tet(p0, p1, p2, p3)
    
    # Put this element matrix contribution into the global stiffness matrix and
    # global right hand side vector
    for a in xrange(0, 4):
      for b in xrange(0, 4):
        id_a = int(LM[i_cell][a])
        id_b = int(LM[i_cell][b])
        if(id_a >= 0 and id_b >= 0):
          globalK_i.append(id_a)
          globalK_j.append(id_b)
          globalK_v.append(kab[a][b])
    
        if(id_a >= 0 and id_b == -1):
          iid = pts_cell.GetId(b)
          assert(dirichlet_BC[iid] > 0.0)
          globalF[id_a] = globalF[id_a] - dirichlet_BC[iid]*kab[a][b]
        
  print('Equations assembled! Making sparse matrix...')
  
  sparseK = sparse.coo_matrix((globalK_v, (globalK_i, globalK_j)),shape=(num_unknowns, num_unknowns)).tocsr()
  
  print('Sparse matrix made! Solving sparse system...')
  
  laplace_sol = spsolve(sparseK, globalF)
  
  print('Sparse solve complete!')
  
  # Save the laplace solution together with the dirichlet BC's and return this
  # combined array
  out_sol = np.zeros(numPoints)
  
  for i_pt in xrange(0, numPoints):
    iid = int(ID[i_pt])
    if(iid >= 0):
      out_sol[i_pt] = laplace_sol[iid]
    else:
      out_sol[i_pt] = dirichlet_BC[i_pt]
  
  return out_sol     
          
# ==============================================================================

def performCMMInflation(wall_vtp, dirichlet_BC):

  numPts = wall_vtp.GetNumberOfPoints()
  numCells = wall_vtp.GetNumberOfCells()
  
  # Define an ID array to store the global equation number of every unknown
  ID = np.zeros((numPts, 3))
  num_unknowns = 0
  for i_wall in xrange(0, numPts):
    if(dirichlet_BC[i_wall] < 0.0):
      ID[i_wall][0] = num_unknowns
      ID[i_wall][1] = num_unknowns+1
      ID[i_wall][2] = num_unknowns+2
      num_unknowns = num_unknowns+3
    else:
      ID[i_wall][0] = -1.0
      ID[i_wall][1] = -1.0
      ID[i_wall][2] = -1.0
      
  print('Number of wall unknowns: ' + str(num_unknowns))
  
  # Define the LM array to store the global ID number for all element nodes
  LM = np.zeros((numCells, 3, 3)) # ElementId, local_node_number, dof_number
  for i_cell in xrange(0, numCells):
    temp_cell = wall_vtp.GetCell(i_cell)
    pts_cell = temp_cell.GetPointIds()
    for i_tri in xrange(0, 3):
      iid = pts_cell.GetId(i_tri)
      LM[i_cell][i_tri][0] = ID[iid][0]
      LM[i_cell][i_tri][1] = ID[iid][1]
      LM[i_cell][i_tri][2] = ID[iid][2]
      
  # Loop over all the elements and add their element stiffness matrix into the
  # global stiffness matrix
  globalK_i = []
  globalK_j = []
  globalK_v = []
  globalF = np.zeros(num_unknowns)
  
  for i_cell in xrange(0, numCells):
  
    temp_cell = wall_vtp.GetCell(i_cell)
    pts_cell = temp_cell.GetPointIds()
    # Get the nodal coordinates
    p0 = temp_cell.GetPoints().GetPoint(0)
    p1 = temp_cell.GetPoints().GetPoint(1)
    p2 = temp_cell.GetPoints().GetPoint(2)
    
    local_stresses = np.zeros(6)
    
    # Compute the element stiffness matrix and element forcing vector
    kab, fa = cmm_stiffness_local(p0, p1, p2, local_stresses)
    
    # Put this element matrix contribution into the global stiffness 
    # matrix and global right hand side vector
    for a in xrange(0, 3):
      for b in xrange(0, 3):
        for a_dof in xrange(0, 3):
          for b_dof in range(0, 3):
            id_a = int(LM[i_cell][a][a_dof])
            id_b = int(LM[i_cell][b][b_dof])
            if(id_a >= 0 and id_b >= 0):
              globalK_i.append(id_a)
              globalK_j.append(id_b)
              globalK_v.append(kab[3*a+a_dof][3*b+b_dof])
  
    for af in xrange(0, 3):
      for af_dof in xrange(0, 3):
        id_af = int(LM[i_cell][af][af_dof])
        if(id_af >= 0):
          globalF[id_af] = globalF[id_af] + fa[3*af+af_dof]
  
  print('CMM Equations assembled! Making sparse matrix...')
  
  sparseK = sparse.coo_matrix((globalK_v, (globalK_i, globalK_j)),shape=(num_unknowns, num_unknowns)).tocsr()
  
  print('CMM sparse matrix made! Solving sparse system...')
  
  disp_sol = spsolve(sparseK, globalF)
  
  print('Sparse solve complete!')
  
  out_sol = np.zeros((numPts, 3))
  
  for i_pt in xrange(0, numPts):
    iid_1 = int(ID[i_pt][0])
    if(iid_1 >= 0):
      out_sol[i_pt][0] = disp_sol[iid_1]
    else:
      out_sol[i_pt][0] = dirichlet_BC[i_pt]
      
    iid_2 = int(ID[i_pt][1])
    if(iid_2 >= 0):
      out_sol[i_pt][1] = disp_sol[iid_2]
    else:
      out_sol[i_pt][1] = dirichlet_BC[i_pt]
      
    iid_3 = int(ID[i_pt][2])
    if(iid_1 >= 0):
      out_sol[i_pt][2] = disp_sol[iid_3]
    else:
      out_sol[i_pt][2] = dirichlet_BC[i_pt]
  
  return out_sol

# ==============================================================================

def computeCMMStresses(wall_vtp, disp_sol):

  print('Computing CMM wall stresses!')
  
  numPts = wall_vtp.GetNumberOfPoints()
  numCells = wall_vtp.GetNumberOfCells()
  assert(len(disp_sol) == numPts)
  assert(len(disp_sol[0]) == 3)
  
  scaling_areas = np.zeros(numPts)
  out_stresses = np.zeros((numPts, 6))
  
  for i_cell in xrange(0, numCells):
  
    temp_cell = wall_vtp.GetCell(i_cell)
    pts_cell = temp_cell.GetPointIds()
    # Get the nodal coordinates
    p0 = temp_cell.GetPoints().GetPoint(0)
    p1 = temp_cell.GetPoints().GetPoint(1)
    p2 = temp_cell.GetPoints().GetPoint(2)
    
    # Compute area of the cell
    tri_area = temp_cell.TriangleArea(p0, p1, p2)
    
    id0 = pts_cell.GetId(0)
    id1 = pts_cell.GetId(1)
    id2 = pts_cell.GetId(2)
    
    # Get the nodal displacements in global coordinates
    d0 = disp_sol[id0]
    d1 = disp_sol[id1]
    d2 = disp_sol[id2]
    
    cell_stress = computeCellStress(p0, p1, p2, d0, d1, d2)
    
    scaling_areas[id0] = scaling_areas[id0] + tri_area
    scaling_areas[id1] = scaling_areas[id1] + tri_area
    scaling_areas[id2] = scaling_areas[id2] + tri_area
    
    for i in xrange(0, 6):
      out_stresses[id0][i] = out_stresses[id0][i] + cell_stress[i]*tri_area
      out_stresses[id1][i] = out_stresses[id1][i] + cell_stress[i]*tri_area
      out_stresses[id2][i] = out_stresses[id2][i] + cell_stress[i]*tri_area
      
  # Scale the outstresses by the scaling area
  for i_pt in xrange(0, numPts):
    for i in xrange(0, 6):
      out_stresses[i_pt][i] = out_stresses[i_pt][i] / scaling_areas[i_pt]
      
  print('Done computing CMM wall stresses!')
  
  return out_stresses

# ==============================================================================

def computeCellStress(p0, p1, p2, d0, d1, d2):

  # Create the rotation matrix by getting a set of coordinates in the element
  # i.e., two coordinate directions will be in the plane of the element and
  # one coordinate will be perpendicularly out
  e0_local = (p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2])
  temp_mag = np.sqrt(e0_local[0]*e0_local[0] + e0_local[1]*e0_local[1] \
                     + e0_local[2]*e0_local[2])
  e0_local = (e0_local[0]/temp_mag, e0_local[1]/temp_mag, e0_local[2]/temp_mag)
  
  etemp_local = (p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2])
  temp_mag = np.sqrt(etemp_local[0]*etemp_local[0] + \
                     etemp_local[1]*etemp_local[1] + \
                     etemp_local[2]*etemp_local[2])
  etemp_local = (etemp_local[0]/temp_mag, etemp_local[1]/temp_mag, \
                 etemp_local[2]/temp_mag)
                 
  e2_local = np.cross(e0_local, etemp_local)
  temp_mag = np.sqrt(e2_local[0]*e2_local[0] + e2_local[1]*e2_local[1] \
                     + e2_local[2]*e2_local[2])
  e2_local = (e2_local[0]/temp_mag, e2_local[1]/temp_mag, e2_local[2]/temp_mag)
  
  e1_local = np.cross(e2_local, e0_local)
  temp_mag = np.sqrt(e1_local[0]*e1_local[0] + e1_local[1]*e1_local[1] \
                     + e1_local[2]*e1_local[2])
  e1_local = (e1_local[0]/temp_mag, e1_local[1]/temp_mag, e1_local[2]/temp_mag)
  
  # Making sure all the local bases are orthogonal
  assert(abs(np.dot(e0_local, e1_local)) < 1e-6)
  assert(abs(np.dot(e0_local, e2_local)) < 1e-6)
  assert(abs(np.dot(e1_local, e2_local)) < 1e-6)
  
  Theta = np.zeros((3, 3))
  Theta[0][0] = e0_local[0]
  Theta[0][1] = e0_local[1]
  Theta[0][2] = e0_local[2]
  Theta[1][0] = e1_local[0]
  Theta[1][1] = e1_local[1]
  Theta[1][2] = e1_local[2]
  Theta[2][0] = e2_local[0]
  Theta[2][1] = e2_local[1]
  Theta[2][2] = e2_local[2]
  Theta_t = np.transpose(Theta)
  
  # Compute the coordinates of the three vertices in the local system using the
  # rotation matrix Theta
  p0_local = np.dot(Theta, p0)
  p1_local = np.dot(Theta, p1)
  p2_local = np.dot(Theta, p2)
  
  # Compute the displacements in the local system using the rotation matrix Theta
  d0_local = np.dot(Theta, d0)
  d1_local = np.dot(Theta, d1)
  d2_local = np.dot(Theta, d2)
  
  # Compute the determinant of the Jacobian matrix and the area of the elements.
  # Here, we have taken the following to be the shape functions:
  #  1. N0 = r
  #  2. N1 = s
  #  3. N2 = 1 - r - s
  detJ = (p0_local[0]-p2_local[0])*(p1_local[1]-p2_local[1]) - \
         (p0_local[1]-p2_local[1])*(p1_local[0]-p2_local[0])
  area = 0.5*detJ
  
  # Making sure we computed the correct area and that our transformation was
  # area preserving
  vtk_area = vtk.vtkTriangle().TriangleArea(p0, p1, p2)
  vtk_area2 = vtk.vtkTriangle().TriangleArea(p0_local, p1_local, p2_local)
  assert(abs(area - vtk_area)/vtk_area < 1e-5)
  assert(abs(area - vtk_area2)/vtk_area2 < 1e-5)
  
  # These are the derivatives of the parent coordinates wrt local coordinates
  dr_dxl = (1.0/detJ)*(p1_local[1]-p2_local[1])
  ds_dxl = (1.0/detJ)*(p2_local[1]-p0_local[1])
  dr_dyl = (1.0/detJ)*(p2_local[0]-p1_local[0])
  ds_dyl = (1.0/detJ)*(p0_local[0]-p2_local[0])
  
  # Now, write the shape function derivatives in terms of these
  dN0_dxl = dr_dxl
  dN0_dyl = dr_dyl
  dN1_dxl = ds_dxl
  dN1_dyl = ds_dyl
  dN2_dxl = -dr_dxl - ds_dxl
  dN2_dyl = -dr_dyl - ds_dyl
  
  # Form the D matrix of material constants
  D_mat = np.zeros((5, 5))
  D_coef = DEFORMABLE_E/(1.0 - POISSON_RATIO*POISSON_RATIO)
  D_mat[0][0] = D_coef
  D_mat[0][1] = POISSON_RATIO*D_coef
  D_mat[1][0] = POISSON_RATIO*D_coef
  D_mat[1][1] = D_coef
  D_mat[2][2] = 0.5*D_coef*(1.0 - POISSON_RATIO)
  D_mat[3][3] = 0.5*D_coef*KCONS*(1.0 - POISSON_RATIO)
  D_mat[4][4] = 0.5*D_coef*KCONS*(1.0 - POISSON_RATIO)
  
  # Compute the strains in reduced vector form: {du/dx, dv/dy, du/dy + dv/dx,
  # dw/dx, dw/dy} all in the local coordinates
  strain_vec = np.zeros(5)
  strain_vec[0] = dN0_dxl*d0_local[0] + dN1_dxl*d1_local[0] + dN2_dxl*d2_local[0] #du_dx
  strain_vec[1] = dN0_dyl*d0_local[1] + dN1_dyl*d1_local[1] + dN2_dyl*d2_local[1] #dv_dy
  strain_vec[2] = (dN0_dyl*d0_local[0] + dN0_dxl*d0_local[1]) + \
           (dN1_dyl*d1_local[0] + dN1_dxl*d1_local[1]) + \
           (dN2_dyl*d2_local[0] + dN2_dxl*d2_local[1]) # tau_xy
  strain_vec[3] = dN0_dxl*d0_local[2] + dN1_dxl*d1_local[2] + dN2_dxl*d2_local[2] # dw_dx
  strain_vec[4] = dN0_dyl*d0_local[2] + dN1_dyl*d1_local[2] + dN2_dyl*d2_local[2] #dw_dy
  
  stress_local = np.dot(D_mat, strain_vec)
  
  # Form the 3x3 stress tensor from this collapsed notation
  stress_tensor_local = np.zeros((3,3))
  stress_tensor_local[0][0] = stress_local[0]
  stress_tensor_local[0][1] = stress_local[2]
  stress_tensor_local[0][2] = stress_local[3]
  stress_tensor_local[1][0] = stress_local[2]
  stress_tensor_local[1][1] = stress_local[1]
  stress_tensor_local[1][2] = stress_local[4]
  stress_tensor_local[2][0] = stress_local[3]
  stress_tensor_local[2][1] = stress_local[4]
  stress_tensor_local[2][2] = 0.0
  
  # Rotate this local stress tensor back into the global coordinates
  stress_tensor_global = np.dot(Theta_t, stress_tensor_local)
  stress_tensor_global = np.dot(stress_tensor_global, Theta)
  
  # Assert that this global stress tensor is symmetric
  assert(abs(stress_tensor_global[0][1] - stress_tensor_global[1][0]) < 1e-8)
  assert(abs(stress_tensor_global[0][2] - stress_tensor_global[2][0]) < 1e-8)
  assert(abs(stress_tensor_global[2][1] - stress_tensor_global[1][2]) < 1e-8)
  
  # These stress components are: {o_xx, o_yy, o_zz, tau_xy, tau_xz, tau_yz}
  stress_global = np.zeros(6)
  stress_global[0] = stress_tensor_global[0][0]
  stress_global[1] = stress_tensor_global[1][1]
  stress_global[2] = stress_tensor_global[2][2]
  stress_global[3] = stress_tensor_global[0][1]
  stress_global[4] = stress_tensor_global[2][0]
  stress_global[5] = stress_tensor_global[2][1]
  
  return stress_global

# ==============================================================================

def performCMMPrestress(wall_vtp, dirichlet_BC):

  numPts = wall_vtp.GetNumberOfPoints()
  numCells = wall_vtp.GetNumberOfCells()
  pre_stresses = np.zeros((numPts, 6))
  
  # Define an ID array to store the global equation number of every unknown
  ID = np.zeros((numPts, 3))
  num_unknowns = 0
  for i_wall in xrange(0, numPts):
    if(dirichlet_BC[i_wall] < 0.0):
      ID[i_wall][0] = num_unknowns
      ID[i_wall][1] = num_unknowns+1
      ID[i_wall][2] = num_unknowns+2
      num_unknowns = num_unknowns+3
    else:
      ID[i_wall][0] = -1.0
      ID[i_wall][1] = -1.0
      ID[i_wall][2] = -1.0
      
  print('Number of wall unknowns: ' + str(num_unknowns))
  
  # Define the LM array to store the global ID number for all element nodes
  LM = np.zeros((numCells, 3, 3)) # ElementId, local_node_number, dof_number
  for i_cell in xrange(0, numCells):
    temp_cell = wall_vtp.GetCell(i_cell)
    pts_cell = temp_cell.GetPointIds()
    for i_tri in xrange(0, 3):
      iid = pts_cell.GetId(i_tri)
      LM[i_cell][i_tri][0] = ID[iid][0]
      LM[i_cell][i_tri][1] = ID[iid][1]
      LM[i_cell][i_tri][2] = ID[iid][2]
      
  # Initialize displacement norm and counter
  disp_norm = 99999.0
  count = 0
  
  disp_norm_tracker = []
  
  while(disp_norm > 0.0001 and count < PS_ITER):
  
    # Loop over all the elements and add their element stiffness matrix into the
    # global stiffness matrix
    globalK_i = []
    globalK_j = []
    globalK_v = []
    globalF = np.zeros(num_unknowns)
    
    for i_cell in xrange(0, numCells):
  
      temp_cell = wall_vtp.GetCell(i_cell)
      pts_cell = temp_cell.GetPointIds()
      # Get the nodal coordinates
      p0 = temp_cell.GetPoints().GetPoint(0)
      p1 = temp_cell.GetPoints().GetPoint(1)
      p2 = temp_cell.GetPoints().GetPoint(2)
    
      # Load in the pre-stresses from the nodes
      local_stresses = np.zeros(6)
      for i_pt in xrange(0, 3):
        iid = pts_cell.GetId(i_pt)
        for i in xrange(0, 6):
          local_stresses[i] = local_stresses[i] + pre_stresses[iid][i]/3.0
    
      # Compute the element stiffness matrix and element forcing vector
      kab, fa = cmm_stiffness_local(p0, p1, p2, local_stresses)
    
      # Put this element matrix contribution into the global stiffness 
      # matrix and global right hand side vector
      for a in xrange(0, 3):
        for b in xrange(0, 3):
          for a_dof in xrange(0, 3):
            for b_dof in range(0, 3):
              id_a = int(LM[i_cell][a][a_dof])
              id_b = int(LM[i_cell][b][b_dof])
              if(id_a >= 0 and id_b >= 0):
                globalK_i.append(id_a)
                globalK_j.append(id_b)
                globalK_v.append(kab[3*a+a_dof][3*b+b_dof])
  
      for af in xrange(0, 3):
        for af_dof in xrange(0, 3):
          id_af = int(LM[i_cell][af][af_dof])
          if(id_af >= 0):
            globalF[id_af] = globalF[id_af] + fa[3*af+af_dof]
  
    print('CMM Equations assembled! Making sparse matrix...')
  
    sparseK = sparse.coo_matrix((globalK_v, (globalK_i, globalK_j)),shape=(num_unknowns, num_unknowns)).tocsr()
  
    print('CMM sparse matrix made! Solving sparse system...')
  
    disp_sol = spsolve(sparseK, globalF)
    
    out_sol = np.zeros((numPts, 3))
  
    for i_pt in xrange(0, numPts):
      iid_1 = int(ID[i_pt][0])
      if(iid_1 >= 0):
        out_sol[i_pt][0] = disp_sol[iid_1]
      else:
        out_sol[i_pt][0] = dirichlet_BC[i_pt]
      
      iid_2 = int(ID[i_pt][1])
      if(iid_2 >= 0):
        out_sol[i_pt][1] = disp_sol[iid_2]
      else:
        out_sol[i_pt][1] = dirichlet_BC[i_pt]
      
      iid_3 = int(ID[i_pt][2])
      if(iid_3 >= 0):
        out_sol[i_pt][2] = disp_sol[iid_3]
      else:
        out_sol[i_pt][2] = dirichlet_BC[i_pt]
  
    print('Sparse solve complete!')
  
    # Check displacement norm
    disp_norm = computeDisplacementNorm(out_sol)
    disp_norm_tracker.append(disp_norm)
    count = count + 1
  
    print('Displacement norm = ' + str(disp_norm))
  
    temp_stresses = computeCMMStresses(wall_vtp, out_sol)
  
    for i_pt in xrange(0, len(temp_stresses)):
      for i_comp in xrange(0, len(temp_stresses[i_pt])):
        pre_stresses[i_pt][i_comp] = pre_stresses[i_pt][i_comp] + temp_stresses[i_pt][i_comp]
        
  return out_sol, pre_stresses, disp_norm_tracker

# ==============================================================================

def writeOutArray(mesh_vtu, array_name, out_solution):

  numPoints = mesh_vtu.GetNumberOfPoints()
  assert(len(out_solution) == numPoints)
  
  out_array = vtk.vtkDoubleArray()
  out_array.SetNumberOfComponents(1)
  out_array.Allocate(numPoints, 10000)
  out_array.SetNumberOfTuples(numPoints)
  out_array.SetName(array_name)
  
  for i_sol in xrange(0, numPoints):
    out_array.SetTuple1(i_sol, out_solution[i_sol])
  
  mesh_vtu.GetPointData().AddArray(out_array)

# ==============================================================================

def writeOutArray3(vtk_model, array_name, out_solution):

  numPoints = vtk_model.GetNumberOfPoints()
  assert(len(out_solution) == numPoints)
  assert(len(out_solution[0]) == 3)
  
  out_array = vtk.vtkDoubleArray()
  out_array.SetNumberOfComponents(3)
  out_array.Allocate(numPoints, 10000)
  out_array.SetNumberOfTuples(numPoints)
  out_array.SetName(array_name)
  
  for i_sol in xrange(0, numPoints):
    out_array.SetTuple3(i_sol, out_solution[i_sol][0], out_solution[i_sol][1], out_solution[i_sol][2])
    
  vtk_model.GetPointData().AddArray(out_array)
  
# ==============================================================================

def writeOutArray6(vtk_model, array_name, out_solution):

  numPoints = vtk_model.GetNumberOfPoints()
  assert(len(out_solution) == numPoints)
  assert(len(out_solution[0]) == 6)
  
  out_array = vtk.vtkDoubleArray()
  out_array.SetNumberOfComponents(6)
  out_array.Allocate(numPoints, 10000)
  out_array.SetNumberOfTuples(numPoints)
  out_array.SetName(array_name)
  
  for i_sol in xrange(0, numPoints):
    out_array.SetTuple6(i_sol, out_solution[i_sol][0], out_solution[i_sol][1], \
      out_solution[i_sol][2], out_solution[i_sol][3], \
      out_solution[i_sol][4], out_solution[i_sol][5])
  
  vtk_model.GetPointData().AddArray(out_array)

# ==============================================================================

def writeDisplacementToVtu(vtu_model, vtp_model, out_name, out_solution):

  numPts_vtu = vtu_model.GetNumberOfPoints()
  numPts_vtp = vtp_model.GetNumberOfPoints()
  
  assert(len(out_solution) == numPts_vtp)
  assert(len(out_solution[0]) == 3)
  
  out_array = vtk.vtkDoubleArray()
  out_array.SetNumberOfComponents(3)
  out_array.Allocate(numPts_vtu, 10000)
  out_array.SetNumberOfTuples(numPts_vtu)
  out_array.SetName(out_name)
  
  # Initialize out_array to zero
  for i_set in xrange(0, numPts_vtu):
    out_array.SetTuple3(i_set, 0.0, 0.0, 0.0)
  
  nodeID_vtu = vtu_model.GetPointData().GetArray('GlobalNodeID')
  nodeID_vtp = vtp_model.GetPointData().GetArray('GlobalNodeID')
  
  for i_vtp in xrange(0, numPts_vtp):
    temp_ID_vtp = nodeID_vtp.GetTuple1(i_vtp)
    for i_vtu in xrange(0, numPts_vtu):
      temp_ID_vtu = nodeID_vtu.GetTuple1(i_vtu)
      if(temp_ID_vtp == temp_ID_vtu):
        temp_disp = out_solution[i_vtp]
        out_array.SetTuple3(i_vtu, temp_disp[0], temp_disp[1], temp_disp[2])
        break
  
  vtu_model.GetPointData().AddArray(out_array)

# ==============================================================================

def writePrestressToVtu(vtu_model, vtp_model, out_name, out_solution):

  numPts_vtu = vtu_model.GetNumberOfPoints()
  numPts_vtp = vtp_model.GetNumberOfPoints()
  
  assert(len(out_solution) == numPts_vtp)
  assert(len(out_solution[0]) == 6)
  
  out_array = vtk.vtkDoubleArray()
  out_array.SetNumberOfComponents(6)
  out_array.Allocate(numPts_vtu, 10000)
  out_array.SetNumberOfTuples(numPts_vtu)
  out_array.SetName(out_name)
  
  # Initialize out_array to zero
  for i_set in xrange(0, numPts_vtu):
    out_array.SetTuple6(i_set, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
  
  nodeID_vtu = vtu_model.GetPointData().GetArray('GlobalNodeID')
  nodeID_vtp = vtp_model.GetPointData().GetArray('GlobalNodeID')
  
  for i_vtp in xrange(0, numPts_vtp):
    temp_ID_vtp = nodeID_vtp.GetTuple1(i_vtp)
    for i_vtu in xrange(0, numPts_vtu):
      temp_ID_vtu = nodeID_vtu.GetTuple1(i_vtu)
      if(temp_ID_vtp == temp_ID_vtu):
        temp_ps0 = out_solution[i_vtp]
        out_array.SetTuple6(i_vtu, temp_ps0[0], temp_ps0[1], temp_ps0[2], \
                       temp_ps0[3], temp_ps0[4], temp_ps0[5])
        break
  
  vtu_model.GetPointData().AddArray(out_array)

# ==============================================================================

if __name__ == '__main__':

  t_start = time.time()

  mesh_vtu, walls_combined_vtp = loadMeshAndWalls()
  numMeshPoints = mesh_vtu.GetNumberOfPoints()
  numMeshCells = mesh_vtu.GetNumberOfCells()
  numWallPoints = walls_combined_vtp.GetNumberOfPoints()
  numWallCells = walls_combined_vtp.GetNumberOfCells()
  
  print('Mesh numPoints: ' + str(numMeshPoints))
  print('Mesh numCells: ' + str(numMeshCells))
  print('Wall numPoints: ' + str(numWallPoints))
  print('Wall numCells: ' + str(numWallCells) + '\n')
  
  # Define numpy arrays for the E and thickness distributions
  dirichlet_BC = np.zeros(numMeshPoints)
  
  for i_sol in xrange(0, numMeshPoints):
    dirichlet_BC[i_sol] = -1.0
    
  # Test applying dirichlet BC on E
  vtp_reader = vtk.vtkXMLPolyDataReader()
  vtp_reader.SetFileName(MESH_COMPLETE_PATH + 'mesh-surfaces/cap_inlet.vtp')
  vtp_reader.Update()
  inflow_vtp = vtk.vtkPolyData()
  inflow_vtp = vtp_reader.GetOutput()
  
  assignDirichletBC(mesh_vtu, inflow_vtp, 500.0, dirichlet_BC)
  
  vtp_reader = vtk.vtkXMLPolyDataReader()
  vtp_reader.SetFileName(MESH_COMPLETE_PATH + 'mesh-surfaces/cap_outlet.vtp')
  vtp_reader.Update()
  outflow_vtp = vtk.vtkPolyData()
  outflow_vtp = vtp_reader.GetOutput()
  
  assignDirichletBC(mesh_vtu, outflow_vtp, 300.0, dirichlet_BC)
  
  # Solve the Laplace equation to solve for the E distribution
  #laplace_sol = performLaplaceSolve(mesh_vtu, dirichlet_BC)
  
  # Save the laplace solution to a vtk array
  #writeOutArray(mesh_vtu, "test_sol", laplace_sol)
  
  # Write out the mesh_vtu to file for checking
  #vtu_writer = vtk.vtkXMLUnstructuredGridWriter()
  #vtu_writer.SetInputData(mesh_vtu)
  #vtu_writer.SetFileName("test_out.vtu")
  #vtu_writer.Write()
  
  # Now write the routines to compute the initial displacement for a CMM formulation
  dirichlet_BC = np.zeros(numWallPoints)
  
  for i_sol in xrange(0, numWallPoints):
    dirichlet_BC[i_sol] = -1.0
    
  assignDirichletBC(walls_combined_vtp, inflow_vtp, 0.0, dirichlet_BC)
  assignDirichletBC(walls_combined_vtp, outflow_vtp, 0.0, dirichlet_BC)
  
  # Solve for the initial displacements
  inflate_sol = performCMMInflation(walls_combined_vtp, dirichlet_BC)
  
  writeOutArray3(walls_combined_vtp, "Displacement", inflate_sol)
  writeDisplacementToVtu(mesh_vtu, walls_combined_vtp, "Displacement", inflate_sol)
  
  # Write out the solution to .vtu file
  vtu_writer = vtk.vtkXMLUnstructuredGridWriter()
  vtu_writer.SetInputData(mesh_vtu)
  vtu_writer.SetFileName("CMM_inflate.vtu")
  vtu_writer.Write()
  
  # Compute the CMM stresses
  #disp_sol, stress_sol, disp_norm_vec = performCMMPrestress(walls_combined_vtp, dirichlet_BC)
  
  # Write out the displacement solution
  #writeOutArray3(walls_combined_vtp, "PS_Displacement", disp_sol)
  
  # Write out pre-stress solution
  #writeOutArray6(walls_combined_vtp, "PS_Prestress", stress_sol)
  
  # Write out the solution to file to investigate
  #vtp_writer = vtk.vtkXMLPolyDataWriter()
  #vtp_writer.SetInputData(walls_combined_vtp)
  #vtp_writer.SetFileName("test_prestress.vtp")
  #vtp_writer.Write()
  
  # Write out the displacement norm history
  #norm_file = open('displacement_norms.dat', 'w')
  
  #for i in xrange(0, len(disp_norm_vec)):
    #write_string = str(disp_norm_vec[i]) + '\n'
    #norm_file.write(write_string)
  
  #norm_file.close()
  
  #writeDisplacementToVtu(mesh_vtu, walls_combined_vtp, "PS_Displacement", disp_sol)
  #writePrestressToVtu(mesh_vtu, walls_combined_vtp, "PS_Prestress", stress_sol)
  
  # Write out the solution to .vtu file
  #vtu_writer = vtk.vtkXMLUnstructuredGridWriter()
  #vtu_writer.SetInputData(mesh_vtu)
  #vtu_writer.SetFileName("CMM_prestress.vtu")
  #vtu_writer.Write()
  
  time_end = time.time()
  time_elapsed = time_end - t_start
  print('Total time: ' + str(time_elapsed))
  
# ==============================================================================



















































