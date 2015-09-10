#/usr/bin/env python

import vtk
import math
from icqsol.shapes.icqShapeManager import ShapeManager

# Create structured grid and save grid in VTK file
nLon, nLat = 16, 8
nLon1, nLat1 = nLon + 1, nLat + 1
structuredGrid = vtk.vtkStructuredGrid()
points = vtk.vtkPoints()
points.SetNumberOfPoints(nLon1 * nLat1)

indx = 0
for i in range(nLon1):
    for j in range(nLat1):
        lon = i * 2. * math.pi / float(nLon)
        lat = -math.pi/2.0 + j * math.pi / float(nLat)
        x = math.cos(lat) * math.cos(lon)
        y = math.cos(lat) * math.sin(lon)
        z = math.sin(lat)
        points.SetPoint(indx, (x, y, z))
        indx += 1

structuredGrid.SetDimensions(1, nLat1, nLon1)
structuredGrid.SetPoints(points)

writer = vtk.vtkStructuredGridWriter()
if vtk.VTK_MAJOR_VERSION >= 6:
    writer.SetInputData(structuredGrid)
else:
    writer.SetInput(structuredGrid)
writer.SetFileName('sphereStructGrid.vtk')
writer.SetFileTypeToASCII()
writer.Write()

# Try reading the grid
shp_manager = ShapeManager(file_format='vtk', 
                           vtk_dataset_type='STRUCTURED_GRID')
shp = shp_manager.loadAsShape('sphereStructGrid.vtk')

