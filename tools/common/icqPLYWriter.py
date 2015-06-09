#!/usr/bin/env python

"""
@brief A class for writing data to a PLY compliant file
@author pletzer@psu.edu
"""

import os
import time
import numpy

class PLYWriter:

  dt2Type = {
      numpy.dtype('float32'): 'float',
      numpy.dtype('float64'): 'double',
      numpy.dtype('int32'): 'int',
      numpy.dtype('uint8'): 'uchar',
  }

  dt2Fmt = {
      numpy.dtype('float32'): '%g',
      numpy.dtype('float64'): '%g',
      numpy.dtype('int32'): '%d',
      numpy.dtype('uint8'): '%d',
  }

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
    Write data
    """
    self._writeHeader()

    # write vertices
    numpy.savetxt(self.fileHandle, self.vertices, fmt='%g')

    # write triangles
    numpy.savetxt(self.fileHandle, self.triangles, fmt='%d')

  def _writeHeader(self):
    """
    Write header
    """

    # global metadata
    print >> self.fileHandle, 'ply'
    print >> self.fileHandle, 'format ascii 1.0'
    userName = os.environ.get('USER', '')
    print >> self.fileHandle, 'comment author: {0}'.format(userName)
    print >> self.fileHandle, 'comment date: {0}'.format(time.asctime())

    # vertices
    numVerts = len(self.vertices)
    print >> self.fileHandle, 'element vertex {0}'.format(numVerts)
    for p in 'x', 'y', 'z':
      print >> self.fileHandle, 'property float {0}'.format(p)

    # triangles
    numTriangs = len(self.triangles)
    print >> self.fileHandle, 'element face {0}'.format(numTriangs)
    print >> self.fileHandle, 'property list uchar int vertex_index'

    # close
    print >> self.fileHandle, 'end_header'    

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
