
import vtk


FileName  =  '/home/ngrandeg/Documents/Toronto_KDProject/KDR01/SimResults/Rigid/all_results.vtu'

#Read vtu file
datareader = vtk.vtkXMLUnstructuredGridReader()
datareader.SetFileName(FileName)
datareader.Update()
data = datareader.GetOutput()

array_name = 'velocity_05000'

#Apply gradient filter
gradient = vtk.vtkGradientFilter()
gradient.SetInputConnection(datareader.GetOutputPort())
gradient.SetInputScalars(vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, array_name)
gradient.SetResultArrayName(array_name + '_gradient')
gradient.Update()


#Write vtu file
out = vtk.vtkXMLUnstructuredGridWriter()
out.SetInputData(gradient.GetOutput())
out.SetFileName('all_results_ext.vtu')
out.Write()