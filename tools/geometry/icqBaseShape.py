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
    self.surfaceMeshes = []

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

  def getPoints(self):
    """
    Get the points
    """
    return self.points

  def getSurfaceMesh(self, faceId):
    """
    Get surface mesh
    @param faceId Id of the surface
    @return node connectivity array
    """
    return self.surfaceMeshes[faceId]

  def save(self, filename):
    """
    Save data in file 
    @param filename file name
    """

    # point array
    xyzArr = vtk.vtkDoubleArray()

    # cell connectivity array 
    cellConnectivity = vtk.vtkIdTypeArray()

    # vertices
    points = vtk.vtkPoints()

    # array of cell Ids
    cellIds = vtk.vtkIdTypeArray()

    # cell connectivity
    cellArr = vtk.vtkCellArray()

    # points and cell connectivity
    polyData = vtk.vtkPolyData()

    numPoints = self.points.shape[0]

    xyzArr.SetNumberOfTuples(numPoints)
    xyzArr.SetNumberOfComponents(3) # 3D
    xyzArr.SetVoidArray(self.points, numPoints*3, 1)

    # set the vertices
    points.SetData(xyzArr)

    # compute the total number of cells
    numCells = 0
    for face in self.surfaceMeshes:
      numCells += face.shape[0]

    # copy all the surface meshes into a single connectivity array
    allCells = numpy.zeros( (numCells, 4), numpy.int )
    allCells[:, 0] = 3 # triangles
    numFaces = len(self.surfaceMeshes)
    cellCount = 0
    for i in range(numFaces):
      sMesh = self.surfaceMeshes[i]
      nCells = sMesh.shape[0]
      iBeg, iEnd = cellCount, cellCount + nCells
      allCells[iBeg:iEnd, 1:] = sMesh[:, :]
      cellCount += nCells

    cellIds.SetVoidArray(allCells, cellCount*4, 1)

    # build the cell connectivity object
    cellArr.SetNumberOfCells(cellCount)
    cellArr.SetCells(cellCount, cellIds)

    # build the polydata object
    polyData.SetPoints(points)
    polyData.SetPolys(cellArr)

    # write to file
    writer = None
    if filename.lower().find('.ply') > 0:
      # PLY format
      writer = vtk.vtkPLYWriter()
    else:
      # will default to VTK
      writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(filename)
    if vtk.VTK_MAJOR_VERSION >= 6:
      writer.SetInputData(polyData)
    else:
      writer.SetInput(polyData)
    writer.Write()
    writer.Update()







