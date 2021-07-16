#!/usr/bin/env python

# --- VTU MODIFICATION FOR CALCULATION OF VISCOUS DISSIPATION FIELD
# --- LEAVE ALL YOUR HATE AT THE LOGIN NODE - AV
# --- 09/15/2017

# This scripts takes the name of a file as input and reads it as a 3D unstructured grid
# Based on the velocity values (named as NS_Velocity), it computes 
# the velocity-gradient tensor and the viscous dissipation field inside the domain

import vtk as vtk
import numpy as np
import os
import sys


def vis_diss(FileName) :

    # ---  CREATE A UnstructuredGrid READER. ASSOCIATE FILE NAME HERE 
    datareader=vtk.vtkUnstructuredGridReader()
    datareader.SetFileName(FileName)
    datareader.Update()

    # Read your data into another UnstructuredGrid variable for manipulation
    model=vtk.vtkUnstructuredGrid()
    model=datareader.GetOutput()


    # Read in velocity field from the vtu 
    temp=vtk.vtkDoubleArray()
    temp=model.GetPointData().GetArray("NS_Velocity")


    # -- Get total number of cells and faces in the model
    ncell = model.GetNumberOfCells()
    npts = model.GetNumberOfPoints()

    # -- Create sum-of-gradients array and a counter array, for computing average gradient at each point
    # Format: gradsum(id,k) = du_j/dx_i where k = 3*i+j
    count = np.zeros(npts)
    gradsum = np.zeros((npts,9))


    # -- MASTER LOOP 
    for icell in xrange(0,ncell):

        currcell = model.GetCell(icell)
        pts_cell = currcell.GetPointIds()

        # -- Increment count at the respective point 
        for i in xrange(0,4):
            count[pts_cell.GetId(i)] = count[pts_cell.GetId(i)] + 1 

        # -- Extract points & velocities for current cell

        p0 = model.GetPoint(pts_cell.GetId(0))
        p1 = model.GetPoint(pts_cell.GetId(1))
        p2 = model.GetPoint(pts_cell.GetId(2))
        p3 = model.GetPoint(pts_cell.GetId(3))


        v0 = temp.GetTuple3(pts_cell.GetId(0))
        v1 = temp.GetTuple3(pts_cell.GetId(1))
        v2 = temp.GetTuple3(pts_cell.GetId(2))
        v3 = temp.GetTuple3(pts_cell.GetId(3))

        # -- Compute simplex gradient using linear solves for each component
        # -- For this, 1. first compute centroid 

        ctr_x = 0.25*(p0[0]+p1[0]+p2[0]+p3[0])
        ctr_y = 0.25*(p0[1]+p1[1]+p2[1]+p3[1])
        ctr_z = 0.25*(p0[2]+p1[2]+p2[2]+p3[2])

        # -- 2. Compute deviations from centroid dx_i,dy_i,dz_i

        dx0 = p0[0]-ctr_x 
        dx1 = p1[0]-ctr_x 
        dx2 = p2[0]-ctr_x 
        dx3 = p3[0]-ctr_x 

        dy0 = p0[1]-ctr_y 
        dy1 = p1[1]-ctr_y 
        dy2 = p2[1]-ctr_y 
        dy3 = p3[1]-ctr_y 

        dz0 = p0[2]-ctr_z 
        dz1 = p1[2]-ctr_z 
        dz2 = p2[2]-ctr_z 
        dz3 = p3[2]-ctr_z 

        # --- 3. Prepare 4x4 matrix for computing gradient

        D = np.matrix([[1,dx0,dy0,dz0],[1,dx1,dy1,dz1],[1,dx2,dy2,dz2],[1,dx3,dy3,dz3]])

        # --- 4. Do linear solve for matrix gradient and centroid velocity values

        # Format: du = [u_ctr,du/dx,du/dy,du/dz] and so on for v and w 
        du = np.linalg.solve(D,[v0[0],v1[0],v2[0],v3[0]])
        dv = np.linalg.solve(D,[v0[1],v1[1],v2[1],v3[1]])
        dw = np.linalg.solve(D,[v0[2],v1[2],v2[2],v3[2]])

        # --- Now use these to update gradsum 
        for i in xrange(0,4):
            gradsum[pts_cell.GetId(i),0] = gradsum[pts_cell.GetId(i),0] + du[1] 
            gradsum[pts_cell.GetId(i),1] = gradsum[pts_cell.GetId(i),1] + du[2] 
            gradsum[pts_cell.GetId(i),2] = gradsum[pts_cell.GetId(i),2] + du[3]
            gradsum[pts_cell.GetId(i),3] = gradsum[pts_cell.GetId(i),3] + dv[1] 
            gradsum[pts_cell.GetId(i),4] = gradsum[pts_cell.GetId(i),4] + dv[2] 
            gradsum[pts_cell.GetId(i),5] = gradsum[pts_cell.GetId(i),5] + dv[3]
            gradsum[pts_cell.GetId(i),6] = gradsum[pts_cell.GetId(i),6] + dw[1] 
            gradsum[pts_cell.GetId(i),7] = gradsum[pts_cell.GetId(i),7] + dw[2] 
            gradsum[pts_cell.GetId(i),8] = gradsum[pts_cell.GetId(i),8] + dw[3]

    # -- END MASTER LOOP

    # -- COMPUTE AVERAGE GRADIENTS -- the following computes a row-by-row division 
    gradavg = gradsum/count[:,None]


    # -- Compute viscous dissipation field -- assume cgs units for using viscosity
    mu = 0.04
    phi = vtk.vtkDoubleArray()
    phi.SetNumberOfComponents(1)
    phi.Allocate(npts,128)
    phi.SetNumberOfTuples(npts)
    phi.SetName("phi")

    for i in xrange(0,npts):

        phitmp = mu*(   2*gradavg[i,0]*gradavg[i,0] \
            +           2*gradavg[i,4]*gradavg[i,4] \
            +           2*gradavg[i,8]*gradavg[i,8] \
            +           pow(gradavg[i,1]+gradavg[i,3],2) \
            +           pow(gradavg[i,2]+gradavg[i,6],2) \
            +           pow(gradavg[i,5]+gradavg[i,7],2) )

        phi.SetTuple1(i,phitmp)

    # -- Store viscous dissipation function into new file 

    FileNew = "Updated_"+'.vtu'.join(FileName.rsplit('.vtk',1))
    model.GetPointData().AddArray(phi)

    modelwrite = vtk.vtkXMLUnstructuredGridWriter()
    # Depending on version, the following command may be SetInput or SetInputData
    modelwrite.SetInput(model)
    modelwrite.SetFileName(FileNew)
    modelwrite.Write()


if __name__ == '__main__':
    # Map command line arguments to function arguments.
    vis_diss(*sys.argv[1:])



