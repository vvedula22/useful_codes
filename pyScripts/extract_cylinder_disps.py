import os
import sys
import numpy as np
import glob
import vtk

result_dir = '../24-procs/vtuFiles/'
result_prefix = 'result_'
displacement_name = 'Displacement'
init_displacement_name = 'Initial_displacement'

# Get all the result files and sort them
result_list = glob.glob(result_dir + result_prefix + '*')
result_list = sorted(result_list)

output_displacements = []
output_net_displacements = []

for i_result in result_list:

  print(i_result)
  result_reader = vtk.vtkXMLUnstructuredGridReader()
  result_reader.SetFileName(i_result)
  result_reader.Update()

  result_model = vtk.vtkUnstructuredGrid()
  result_model = result_reader.GetOutput()

  # Extract out the surface to get the shell
  surface_extractor = vtk.vtkDataSetSurfaceFilter()
  surface_extractor.SetInputData(result_model)
  surface_extractor.Update()

  surface_model = vtk.vtkPolyData()
  surface_model = surface_extractor.GetOutput()

  # Slice the vessel in the middle
  cut_origin = (0.0, 0.0, 15.0)
  cut_direction = (0.0, 0.0, 1.0)

  cutPlane = vtk.vtkPlane()
  cutPlane.SetOrigin(cut_origin)
  cutPlane.SetNormal(cut_direction)
  cutter = vtk.vtkCutter()
  cutter.SetCutFunction(cutPlane)
  cutter.SetInputData(surface_model)
  cutter.Update()

  cut_model = vtk.vtkPolyData()
  cut_model = cutter.GetOutput()
  cut_model_numPts = cut_model.GetNumberOfPoints()
  cut_model_displacements  = cut_model.GetPointData().GetArray(displacement_name)
  cut_model_displacements0 = cut_model.GetPointData().GetArray(init_displacement_name)

  # Extract out the displacements and save them into a list
  cut_disp = []
  cut_net_disp = []
  for i_node in xrange(0, cut_model_numPts):
    temp_disp = cut_model_displacements.GetTuple3(i_node)
    temp_disp0 = cut_model_displacements0.GetTuple3(i_node)
    temp_disp_mag = np.sqrt(temp_disp[0]*temp_disp[0] + temp_disp[1]*temp_disp[1] + temp_disp[2]*temp_disp[2])
    temp_net_disp_mag = np.sqrt((temp_disp[0]-temp_disp0[0])*(temp_disp[0]-temp_disp0[0]) \
        + (temp_disp[1]-temp_disp0[1])*(temp_disp[1]-temp_disp0[1]) \
        + (temp_disp[2]-temp_disp0[2])*(temp_disp[2]-temp_disp0[2]))
    cut_disp.append(temp_disp_mag)
    cut_net_disp.append(temp_net_disp_mag)

  # Get the average displacement and save to output vector
  mean_disp = np.mean(cut_disp)
  mean_net_disp = np.mean(cut_net_disp)
  output_displacements.append(mean_disp)
  output_net_displacements.append(mean_net_disp)

# Write out to file
output_file = open('cylinder_displacements.dat', 'w')

for i in xrange(0, len(output_displacements)):
  write_string = str(output_displacements[i]) + ' ' + str(output_net_displacements[i]) + '\n'
  output_file.write(write_string)

output_file.close()




























