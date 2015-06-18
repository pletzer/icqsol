#!/usr/bin/env python

"""
@brief A class for writing data to a PLY compliant file
@author pletzer@psu.edu
"""

import os
import time
import numpy

class PLYWriter:

  def __init__(self, filename):

    self.fileHandle = file(filename, 'w')

    self.vertices = []
    self.triangles = []

  def setVertices(self, verts):
    """
    Set vertex data
    @param verts array of x, y, z vertices
    """
    self.vertices = verts

  def setTriangles(self, triangs):
    """
    Set triangle data
    @param triangs vertex indices (triangles)
    """
    self.triangles = triangs

  def write(self):
    """
    Write data to file
    """

    pointData = None
    if self.vertices.dtype == numpy.float64:
      pointData = vtk.vtkDoubleArray()
    else:
      pointData = vtk.vtkFloatArray()
    numPoints, numDims = self.vertices.shape
    pointData.SetNumberOfTuples(numPoints)
    pointData.SetNumberOfComponents(numDims)
    pointData.SetVoidArray(self.vertices, numPoints*numDims, 1)

    points = vtk.vtkPoints()
    points.SetData(pointData)

    pd = vtkPolyDataMapper()


################################################################################

def test():

  numPoints = 3
  verts = numpy.zeros((numPoints,), dtype=[('x', 'f4'), 
                                           ('y', 'f4'), 
                                           ('z', 'f4'),
                                           ])
  verts['x'] = [0., 0., 1.]
  verts['y'] = [0., 1., 0.]
  verts['z'] = [1., 0., 0.]

  numTriangs = 1
  triangs = numpy.zeros((numTriangs,), dtype=[('ia', 'i4'),
                                              ('ib', 'i4'),
                                              ('ic', 'i4'), 
                                              ('r', 'uint8'),
                                              ('g', 'uint8'),
                                              ('b', 'uint8')
                                             ])
  triangs['ia'][0] = 0
  triangs['ib'][0] = 1
  triangs['ic'][0] = 2

  pw = PLYWriter('test.ply')
  pw.setVertices(verts)
  pw.setTriangles(triangs)
  pw.write()

if __name__ == '__main__':
  test()
