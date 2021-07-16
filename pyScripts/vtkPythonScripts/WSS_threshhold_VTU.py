#!/usr/bin/env python

# --- VTK UNSTRUCTURED GRID MODIFICATION FOR CALCULATION OF LOW-WSS REGIONS -- BASED OFF A SCRIPT BY JUSTIN TRAN 
# --- 05/09/2017 AV

# This script computes the (raw) area of the regions below a certain wall shear stress threshhold
# A simple binary threshhold is used in this script

# It needs three inputs: input VTK (VTU) filename, cut-off threshhold for WSS, WSS array name 

import sys


def bound_box(p0,p1,p2) :
	# --- Defines a geometrical domain of interest (here: a box)
	# If true: point lies within domain of interest ; else false

	bdy_low = -0.5
	bdy_upp =  0.8

	bdz_low =  0.8
	bdz_upp =  5.0

	y0 = (bdy_low <= p0[1] <= bdy_upp)
	y1 = (bdy_low <= p1[1] <= bdy_upp)
	y2 = (bdy_low <= p2[1] <= bdy_upp)

	z0 = (bdz_low <= p0[2] <= bdz_upp)
	z1 = (bdz_low <= p1[2] <= bdz_upp)
	z2 = (bdz_low <= p2[2] <= bdz_upp)


	if( y0 and y1 and y2 and z0 and z1 and z2) :
		return True
	else :
		return False


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


def boundary_extract(model) : 


	numPts=model.GetNumberOfPoints()
	numCells=model.GetNumberOfCells()
	
	# Note : A VTK Cell is equivalent ot an element


	# ====================================================
	#	Find boundary faces 
	# ====================================================	

	face_list = []

	# --- Loop over number of cells 
	for i_cell in xrange(0,numCells):

		temp_cell = model.GetCell(i_cell)
		pts_cell = temp_cell.GetPointIds()

		# Get sorted list of IDs in the element

		local_IDs =  [int(pts_cell.GetId(0)), int(pts_cell.GetId(1)), int(pts_cell.GetId(2)), int(pts_cell.GetId(3))]
		local_IDs = sorted(local_IDs)

		# Add all the faces of the element to the face list 

		face_list.append((local_IDs[0],local_IDs[1],local_IDs[2]))
		face_list.append((local_IDs[0],local_IDs[1],local_IDs[3]))		
		face_list.append((local_IDs[0],local_IDs[2],local_IDs[3]))
		face_list.append((local_IDs[1],local_IDs[2],local_IDs[3]))		


	# Sort face list
	face_list = sorted(face_list)


	bd_face = []
	counter = 0

	# Remove faces that are shared between elements (internal faces) and filter out boundary faces

	# Note that this loop because each face is shared by atmost two elements 

	while(counter < len(face_list) -1):

		face1 = face_list[counter]
		face2 = face_list[counter+1]


		if ( (face1[0]==face2[0]) and (face1[1]==face2[1]) and (face1[2]==face2[2]) ) : 

			counter = counter + 1
		else :

			bd_face.append((face1[0],face1[1],face1[2]))


		counter = counter + 1 

	return bd_face



def wss_thresh(FileName,thresh,wssname) :

	import vtk as vtk
	import numpy as np
	
	th = float(thresh)

	# ---  CREATE A UnstructuredGrid READER. ASSOCIATE FILE NAME HERE 
	datareader=vtk.vtkXMLUnstructuredGridReader()
	datareader.SetFileName(FileName)
	datareader.Update()
	
	# Read your data into another UnstructuredGrid variable for manipulation
	model=vtk.vtkUnstructuredGrid()
	model=datareader.GetOutput()

	
	# Read your fav variable from the vtu 
	temp=vtk.vtkDoubleArray()
	temp=model.GetPointData().GetArray(wssname)

	# # Read coordinate array
	# coord=model.GetPoints()

	# print coord

	# Extract boundary from vtu
	bd_face = boundary_extract(model)

	nfaces = len(bd_face)


	
	totalarea = 0;
	lowwssarea = 0;

	# --- LOOP THROUGH CELLS 
	
	for iface in xrange(0,nfaces):

		face = bd_face[iface]

		id0 = face[0]
		id1 = face[1]
		id2 = face[2]

		p0 = model.GetPoint(id0)
		p1 = model.GetPoint(id1)
		p2 = model.GetPoint(id2)

		wss0 = np.array(temp.GetTuple3(id0))
		wss1 = np.array(temp.GetTuple3(id1))
		wss2 = np.array(temp.GetTuple3(id2))

		wss_avg_face = np.linalg.norm(wss0)+ np.linalg.norm(wss1) + np.linalg.norm(wss2)

		wss_avg_face = wss_avg_face/3.0

		flag_wss = 0

		if (wss_avg_face >= th):
			flag_wss = 1

		if(bound_box(p0,p1,p2)):
			totalarea = totalarea + area_calc(p0,p1,p2)

		if(flag_wss == 0 and bound_box(p0,p1,p2)):
			lowwssarea = lowwssarea + area_calc(p0,p1,p2)


	# print lowwssarea
	print totalarea

	return (lowwssarea)




if __name__ == '__main__':
    # Map command line arguments to function arguments.
    print wss_thresh(*sys.argv[1:])