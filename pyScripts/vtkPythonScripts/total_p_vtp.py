#!/usr/bin/env python

# --- VTK-PYTHON SCRIPT FOR READING VTP OF MESH WITH FLOW & PRESSURE ARRAYS 
# --- AND COMPUTING INSTANTANEOUS TOTAL PRESSURE FLUX AT A TARGET SURFACE (VTP) -- AV 03/05/18

# TOTAL PRESSURE FLUX IN A CLOSED SET OF SURFACES IS EQUIVALENT TO ENERGY DISSIPATION
# For formula, refer to the appendix in :
# "A New Multiparameter Approach to Computational Simulation for Fontan Assessment and Redesign"
# Marsden, A., Reddy, V.M., Shadden, S.C., Chan, F.P., Taylor, C.A., Feinstein, J.A.
# Congenital Heart Disese, 2010, Pgs. 104-17 

# We assume the inputs are in consistent units to compute total pressures

import sys
import numpy as np

def area_calc(p0,p1,p2) :

        # -- AREA = 0.5*|AB|*|AC|*sin(theta)
        # Theta -- angle b/w AB and AC

        # A : p0, B : p1 , C : p2 

        import numpy as np

        q0 = np.array(p0)
        q1 = np.array(p1)
        q2 = np.array(p2)

        dq1 = q1-q0
        dq2 = q2-q0

        costh = np.dot(dq1,dq2)/((np.linalg.norm(dq1))*(np.linalg.norm(dq2)))
        sinth = np.sqrt(1.0 - costh*costh)

        return 0.5*(np.linalg.norm(dq1))*(np.linalg.norm(dq2))*sinth



# Compute centroid value using nodal values
def centval_calc(p0,p1,p2,v0,v1,v2) :

	# Determine centroid coordinates
	x0 = (p0[0]+p1[0]+p2[0])/3.0
	x1 = (p0[1]+p1[1]+p2[1])/3.0
	x2 = (p0[2]+p1[2]+p2[2])/3.0

	x = np.array([x0,x1,x2])

	# Barycentric interpolation 

	a0 = area_calc(x,p1,p2)
	a1 = area_calc(p0,x,p2)
	a2 = area_calc(p0,p1,x)

	# Return weighted sum, based on deviations
	return (v0*a0 + v1*a1 + v2*a2)/(a0+a1+a2)


# Main total pressure flux computation routine

def ptot_flux(FlowFileName,SurfFileName,velname,pressurename):

	import vtk as vtk
	
	rho = 1.06
	
	# Read your flow data into another polydata variable for manipulation
	datareader=vtk.vtkXMLPolyDataReader()
	datareader.SetFileName(FlowFileName)
	datareader.Update()
	model=vtk.vtkPolyData()
	model=datareader.GetOutput()
	

	# Read your surf-of-interest data into another polydata variable for manipulation
	datareader_surf=vtk.vtkXMLPolyDataReader()
	datareader_surf.SetFileName(SurfFileName)
	datareader_surf.Update()
	surf=vtk.vtkPolyData()
	surf=datareader_surf.GetOutput()

	numPts=model.GetNumberOfPoints()
	numCells=model.GetNumberOfCells()
	surf_numPts=surf.GetNumberOfPoints()
	surf_numCells=surf.GetNumberOfCells()

	
	# Read your fav variables from the model vtp 
	temp_pressure=vtk.vtkDoubleArray()
	temp_pressure=model.GetPointData().GetArray(pressurename)
	temp_vel=vtk.vtkDoubleArray()
	temp_vel=model.GetPointData().GetArray(velname)
	model_id = model.GetPointData().GetArray('GlobalNodeID')
	surf_id = surf.GetPointData().GetArray('GlobalNodeID')
  
	# print model_id

 	# Find the normal vector to this face
 	normalGenerator = vtk.vtkPolyDataNormals()
 	normalGenerator.SetInput(surf)
 	normalGenerator.ComputeCellNormalsOff()
 	normalGenerator.ComputePointNormalsOn()
 	normalGenerator.Update()
 	normals_test = normalGenerator.GetOutput()
 	normals_arr = normals_test.GetPointData().GetArray('Normals')

	# --- LOOP THROUGH SURFACE POINTS FOR MODEL AND SURF OF INTEREST
	
	ptot_flux = np.zeros(surf_numPts)

	for ipts_surf in xrange(0,surf_numPts):
		temp_sid = surf_id.GetTuple1(ipts_surf)
		for ipts_model in xrange(0,numPts):
			if temp_sid == model_id.GetTuple1(ipts_model):
			 	pt_normal = np.asarray(normals_arr.GetTuple3(ipts_surf))
				pt_p = temp_pressure.GetTuple1(ipts_model)
				pt_vel = np.asarray(temp_vel.GetTuple3(ipts_model))
				# print pt_p, pt_vel, pt_normal
				ptot_flux[ipts_surf] = (pt_p + 0.5*rho*sum(pt_vel*pt_vel))*np.dot(pt_vel,pt_normal)
				break

	# --- END LOOP	

	totalflux = 0.0

	# INTEGRATE PTOT_FLUX POINT DATA INTO A SINGLE VARIABLE - CELL LEVEL LOOP
	for icell in xrange(0,surf_numCells):
		temp_cell = surf.GetCell(icell)
		pts_cell_id = temp_cell.GetPointIds()
		vtkpt = temp_cell.GetPoints()

		p0 = vtkpt.GetPoint(0)
		p1 = vtkpt.GetPoint(1)
		p2 = vtkpt.GetPoint(2)		

		v = []

		for kpt in range(pts_cell_id.GetNumberOfIds()) :
			kid = pts_cell_id.GetId(kpt)
			v.append(ptot_flux[kid])

		# print v, temp_cell.TriangleArea(p0,p1,p2)


		totalflux = totalflux + (centval_calc(p0,p1,p2,v[0],v[1],v[2]))*area_calc(p0,p1,p2)

	return totalflux


if __name__ == '__main__':
    # Map command line arguments to function arguments.
    print ptot_flux(*sys.argv[1:])
