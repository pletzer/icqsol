#!/usr/bin/env python

"""
@brief An abstract class for constructing shapes
@author pletzer@psu.edu
"""

# standard python modules
import re

# extensions
import vtk
import numpy

class BaseShape:

  """
  Base class for shapes
  """

  def __init__(self):
    """
    Constructor
    """

    self.points = []
    self.surfaceMesh = {}
    self.volumeMesh = []

  def evaluate(self, points):
    """
    Evaluate whether points are inside/outside
    @return > 0 for points inside, 0 for points outside
    """
    raise NotImplementedError

  def computeSurfaceMeshes(self, maxTriArea):
    """
    Compute the surface meshes
    @param maxTriArea maximum triangle area
    """
    raise NotImplementedError

  def computeVOlumeMesh(self, maxTetVol):
    """
    Compute the volume mesh 
    @param maxTetVol maximum tetrahedron volume
    """
    raise NotImplementedError

  def getPoints(self):
    """
    Get the points
    """
    return self.points

  def getSurfaceMesh(self, name):
    """
    Get surface mesh
    @param name of the of the surface
    @return node connectivity array
    """
    return self.surfaceMesh

  def getVolumeMesh(self):
    """
    Get the volume mesh
    @return node connectivity array
    """
    return self.volumeMesh

  def save(self, plyFilename):
    """
    Save data in PLY file 
    @param plyFilename PLY file
    """

    xyzArr = vtk.vtkDoubleArray()
    points = vtk.vtkPoints()
    cellArr = vtk.vtkCellArray()
    surfPolyData = vtk.vtkPolyData()

    numPoints = self.points.shape[0]

    xyzArr.SetNumberOfTuples(numPoints)
    xyzArr.SetNumberOfComponents(3)
    xyzArr.SetVoidArray(self.points, numPoints*3, 1)
    points.SetData(xyzArr)

    numCells = self.surfaceMesh.shape[0]

    cell = vtk.vtkTriangle()
    cellArr.Allocate(numCells, 1)
    for i in range(numCells):
      cell.GetPointIds().SetId(0, self.surfaceMesh[i, 0])
      cell.GetPointIds().SetId(1, self.surfaceMesh[i, 1])
      cell.GetPointIds().SetId(2, self.surfaceMesh[i, 2])
      cellArr.InsertNextCell(cell.GetCellType(), cell.GetPointIds())

    surfPolyData.SetPoints(points)
    surfPolyData.SetPolys(cellArr)

    writer = vtk.vtkPLYWriter()
    writer.SetFileName(plyFilename)
    if vtk.VTK_MAJOR_VERSION >= 6:
      writer.SetInputData(surfPolyData)
    else:
      writer.SetInput(surfPolyData)
    writer.SetInput(surfPolyData)
    writer.Write()







