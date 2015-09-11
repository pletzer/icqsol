#/usr/bin/env python

import vtk
import math
from icqsol.shapes.icqShapeManager import ShapeManager
from icqsol import util

# Create unstructured grid and save grid in VTK file
nLon, nLat = 16, 8
nLon1, nLat1 = nLon + 1, nLat + 1

points = vtk.vtkPoints()
points.SetNumberOfPoints(nLon * nLat1)
indx = 0
for j in range(nLat1):
    for i in range(nLon):
        lon = i * 2. * math.pi / float(nLon)
        lat = -math.pi/2.0 + j*math.pi/float(nLat)
        x = math.cos(lat) * math.cos(lon)
        y = math.cos(lat) * math.sin(lon)
        z = math.sin(lat)
        points.SetPoint(indx, (x, y, z))
        indx += 1

unstructuredGrid = vtk.vtkUnstructuredGrid()
unstructuredGrid.SetPoints(points)

quad = vtk.vtkQuad()
ptIds = vtk.vtkIdList()
ptIds.SetNumberOfIds(4)
for j0 in range(nLat):
    j1 = j0 + 1
    for i0 in range(nLon):
        i1 = (i0 + 1) % nLon
        ptIds.SetId(0, j0*nLon + i0)
        ptIds.SetId(1, j0*nLon + i1)
        ptIds.SetId(2, j1*nLon + i1)
        ptIds.SetId(3, j1*nLon + i0)
        unstructuredGrid.InsertNextCell(quad.GetCellType(), ptIds)

writer = vtk.vtkUnstructuredGridWriter()
if vtk.VTK_MAJOR_VERSION >= 6:
    writer.SetInputData(unstructuredGrid)
else:
    writer.SetInput(unstructuredGrid)
writer.SetFileName('sphereUnstructGrid.vtk')
writer.SetFileTypeToASCII()
writer.Write()

# Try reading the grid
shp_manager = ShapeManager(file_format=util.VTK_FORMAT,  vtk_dataset_type=util.UNSTRUCTURED_GRID)
shp = shp_manager.loadAsShape('sphereUnstructGrid.vtk')
