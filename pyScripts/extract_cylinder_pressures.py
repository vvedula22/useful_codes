import os
import sys
import numpy as np
import glob
import vtk

result_dir = '../24-procs/vtuFiles/'
result_prefix = 'result_'
pressure_name = 'Pressure'

# Get all the result files and sort them
result_list = glob.glob(result_dir + result_prefix + '*')
result_list = sorted(result_list)

pressure_output_1 = []
pressure_output_2 = []
pressure_output_3 = []

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
  slice_model_pressure = slice_model.GetPointData().GetArray(pressure_name)

  pressure_list = []
  for i in xrange(0, slice_model_numPts):
    temp_pressure = slice_model_pressure.GetTuple1(i)
    pressure_list.append(temp_pressure)

  avg_pressure = np.mean(pressure_list)
  pressure_output_1.append(avg_pressure)

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
  slice_model_pressure = slice_model.GetPointData().GetArray(pressure_name)

  pressure_list = []
  for i in xrange(0, slice_model_numPts):
    temp_pressure = slice_model_pressure.GetTuple1(i)
    pressure_list.append(temp_pressure)

  avg_pressure = np.mean(pressure_list)
  pressure_output_2.append(avg_pressure)

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
  slice_model_pressure = slice_model.GetPointData().GetArray(pressure_name)

  pressure_list = []
  for i in xrange(0, slice_model_numPts):
    temp_pressure = slice_model_pressure.GetTuple1(i)
    pressure_list.append(temp_pressure)

  avg_pressure = np.mean(pressure_list)
  pressure_output_3.append(avg_pressure)

# Write out to file

pressure_output_file = open('cylinder_pressures.dat', 'w')

for i in xrange(0, len(pressure_output_3)):
  write_string = str(pressure_output_1[i]) + ' ' + str(pressure_output_2[i]) + ' ' + str(pressure_output_3[i]) + '\n'
  pressure_output_file.write(write_string)

pressure_output_file.close()

