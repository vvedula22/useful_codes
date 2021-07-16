import os
import sys
import numpy as np
import glob
import vtk

result_dir = '../24-procs/vtuFiles/'
result_prefix = 'result_'
velocity_name = 'Velocity'

# Get all the result files and sort them
result_list = glob.glob(result_dir + result_prefix + '*')
result_list = sorted(result_list)

result_flows_1 = []
result_flows_2 = []
result_flows_3 = []

for i_result in result_list:

  print(i_result)
  result_reader = vtk.vtkXMLUnstructuredGridReader()
  result_reader.SetFileName(i_result)
  result_reader.Update()

  result_model = vtk.vtkUnstructuredGrid()
  result_model = result_reader.GetOutput()

  # Make a slice near the inlet of the model
  cut_origin = (0.0, 0.0, 0.25)
  cut_direction = (0.0, 0.0, 1.0)

  cutPlane = vtk.vtkPlane()
  cutPlane.SetOrigin(cut_origin)
  cutPlane.SetNormal(cut_direction)
  cutter = vtk.vtkCutter()
  cutter.SetCutFunction(cutPlane)
  cutter.SetInputData(result_model)
  cutter.Update()

  slice_model = vtk.vtkPolyData()
  slice_model = cutter.GetOutput()
  slice_model_numPts = slice_model.GetNumberOfPoints()
  slice_model_numCells = slice_model.GetNumberOfCells()
  slice_model_velocities = slice_model.GetPointData().GetArray(velocity_name)

  # Dot product the velocity with the normal vector to and integrate to compute flow
  total_area = 0.0
  total_flow = 0.0
  for i_cell in xrange(0, slice_model_numCells):

    temp_cell = slice_model.GetCell(i_cell)
    pts_cell = temp_cell.GetPointIds()
    p0 = temp_cell.GetPoints().GetPoint(0)
    p1 = temp_cell.GetPoints().GetPoint(1)
    p2 = temp_cell.GetPoints().GetPoint(2)

    local_area = vtk.vtkTriangle().TriangleArea(p0, p1, p2)
    total_area = total_area + local_area

    local_flow = 0.0
    for ipt in xrange(0, pts_cell.GetNumberOfIds()):
      iid = pts_cell.GetId(ipt)
      nodal_velocity = slice_model_velocities.GetTuple3(iid)
      normal_velocity = -1.0*nodal_velocity[2]
      local_flow = local_flow + normal_velocity

    local_flow = local_flow * local_area / 3.0
    total_flow = total_flow + local_flow

  result_flows_1.append(total_flow)

  # Slice the vessel in the middle
  cut_origin = (0.0, 0.0, 15.0)
  cut_direction = (0.0, 0.0, 1.0)

  cutPlane = vtk.vtkPlane()
  cutPlane.SetOrigin(cut_origin)
  cutPlane.SetNormal(cut_direction)
  cutter = vtk.vtkCutter()
  cutter.SetCutFunction(cutPlane)
  cutter.SetInputData(result_model)
  cutter.Update()

  slice_model = vtk.vtkPolyData()
  slice_model = cutter.GetOutput()
  slice_model_numPts = slice_model.GetNumberOfPoints()
  slice_model_numCells = slice_model.GetNumberOfCells()
  slice_model_velocities = slice_model.GetPointData().GetArray(velocity_name)

  # Search for the maximum velocity
  total_area = 0.0
  total_flow = 0.0
  for i_cell in xrange(0, slice_model_numCells):

    temp_cell = slice_model.GetCell(i_cell)
    pts_cell = temp_cell.GetPointIds()
    p0 = temp_cell.GetPoints().GetPoint(0)
    p1 = temp_cell.GetPoints().GetPoint(1)
    p2 = temp_cell.GetPoints().GetPoint(2)

    local_area = vtk.vtkTriangle().TriangleArea(p0, p1, p2)
    total_area = total_area + local_area

    local_flow = 0.0
    for ipt in xrange(0, pts_cell.GetNumberOfIds()):
      iid = pts_cell.GetId(ipt)
      nodal_velocity = slice_model_velocities.GetTuple3(iid)
      normal_velocity = 1.0*nodal_velocity[2]
      local_flow = local_flow + normal_velocity

    local_flow = local_flow * local_area / 3.0
    total_flow = total_flow + local_flow

  result_flows_2.append(total_flow)

  # Slice the vessel at the end
  cut_origin = (0.0, 0.0, 29.75)
  cut_direction = (0.0, 0.0, 1.0)

  cutPlane = vtk.vtkPlane()
  cutPlane.SetOrigin(cut_origin)
  cutPlane.SetNormal(cut_direction)
  cutter = vtk.vtkCutter()
  cutter.SetCutFunction(cutPlane)
  cutter.SetInputData(result_model)
  cutter.Update()

  slice_model = vtk.vtkPolyData()
  slice_model = cutter.GetOutput()
  slice_model_numPts = slice_model.GetNumberOfPoints()
  slice_model_numCells = slice_model.GetNumberOfCells()
  slice_model_velocities = slice_model.GetPointData().GetArray(velocity_name)

  # Search for the maximum velocity
  total_area = 0.0
  total_flow = 0.0
  for i_cell in xrange(0, slice_model_numCells):

    temp_cell = slice_model.GetCell(i_cell)
    pts_cell = temp_cell.GetPointIds()
    p0 = temp_cell.GetPoints().GetPoint(0)
    p1 = temp_cell.GetPoints().GetPoint(1)
    p2 = temp_cell.GetPoints().GetPoint(2)

    local_area = vtk.vtkTriangle().TriangleArea(p0, p1, p2)
    total_area = total_area + local_area

    local_flow = 0.0
    for ipt in xrange(0, pts_cell.GetNumberOfIds()):
      iid = pts_cell.GetId(ipt)
      nodal_velocity = slice_model_velocities.GetTuple3(iid)
      normal_velocity = 1.0*nodal_velocity[2]
      local_flow = local_flow + normal_velocity

    local_flow = local_flow * local_area / 3.0
    total_flow = total_flow + local_flow

  result_flows_3.append(total_flow)

# Write out to file

centerline_output_file = open('cylinder_flows.dat', 'w')

for i in xrange(0, len(result_flows_3)):
  write_string = str(result_flows_1[i]) +  ' ' + str(result_flows_2[i]) +  ' ' + str(result_flows_3[i]) +  '\n'
  centerline_output_file.write(write_string)

centerline_output_file.close()

