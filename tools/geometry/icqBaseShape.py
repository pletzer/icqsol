#!/usr/bin/env python

"""
@brief A base class for constructing shapes
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

    self.points = numpy.array([])
    self.surfaceMeshes = []
    self.surfaceNormals = []

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

  def getSurfaceNormals(self, faceId):
    """
    Get the surface normal vectors
    @param faceId Id of the surface
    @return array of cell area normals
    @note the normals are not normalized
    """
    return self.surfaceNormals[faceId]

  def computeSurfaceNormals(self):
    """
    Compute the normals for each surface triangle
    @note the vectors are NOT normalized, the magnitudes are the 
    cell areas
    """

    if self.surfaceNormals:
      # nothing to do
      return 

    for face in self.surfaceMeshes:

      # index 0 is cell index
      # index 1 is vertex index
      # index 2 is axis index
      verts = self.points[face]

      verts0 = verts[:, 0, :]
      verts1 = verts[:, 1, :] - verts0
      verts2 = verts[:, 2, :] - verts0
      # verts1 and verts2 are the vector distances from the triangle
      # base (point 0) to points 1 and 2

      normals = numpy.zeros( (verts.shape[0], 3), numpy.float64 )
      # cross product
      normals[:, 0] = verts1[:, 1]*verts2[:, 2] - verts1[:, 2]*verts2[:, 1]
      normals[:, 1] = verts1[:, 2]*verts2[:, 0] - verts1[:, 0]*verts2[:, 2]
      normals[:, 2] = verts1[:, 0]*verts2[:, 1] - verts1[:, 1]*verts2[:, 0]

      self.surfaceNormals.append(normals)

  def save(self, filename):
    """
    Save data in file 
    @param filename file name, either a PLY or VTK file (suffix will determine
                    which file format will be used)
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

  def load(self, filename):
    """
    Load object from file
    @param filename file name, either a PLY or VTK file (suffix will determine
                    which file format will be used)
    """

    reader = None
    if filename.lower().find('ply') > 0:
      # PLY format
      reader = vtk.vtkPLYReader()
    else:
      # will default to VTK
      reader = vtk.vtkPolyDataReader()

    reader.SetFileName(filename)

    pdata = reader.GetOutput()
    pts = pdata.GetPoints()
    cells = pdata.GetPolys()

    numPoints = pts.GetNumberOfPoints()
    numCells = cells.GetNumberOfCells()

    self.points = numpy.zeros( (numPoints, 3), numpy.float64 )
    allTriangles = numpy.zeros( (numCells, 3), numpy.int )

    for i in range(numPoints):
      self.points[i, :] = pts.GetPoint(i)

    cells.InitTraversal()
    ptIds = vtk.vtkIdList()
    for i in range(numCells):
      cells.GetNextCell(ptIds)
      allTriangles[i, :] = ptIds.GetId(0), ptIds.GetId(1), ptIds.GetId(2)

    self.surfaceMeshes = [allTriangles]





