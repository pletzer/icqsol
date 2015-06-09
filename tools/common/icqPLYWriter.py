#!/usr/bin/env python

"""
@brief A class for writing data to a PLY compliant file
@author pletzer@psu.edu
"""

import os
import time

class PLYWriter:

	def __init__(self, filename):

		self.fileHandle = file(filename, 'w')

    self.vertices = []
    self.triangles = []

  def setVertices(self, verts):
    """
    Set vertex data
    @param verts array of x, y, z vertices and corresponding propeerties as a 
                 1-d array of structs. Example:
                 numpy.zeros((n,), dtype = [('x', 'f4'), 
                                            ('y', 'f4'), 
                                            ('z', 'f4'), 
                                            ('color', '|S10')])
    """
    self.vertices = verts

  def setTriangles(self, triangs):
    """
    Set triangle data
    @param triangs vertex indices (triangles) and properties as a 1-d array 
           of structs. Example: 
           numpy.zeros((n,), dtype = [('ia', 'i4'), ('ib', 'i4'), ('ic', 'i4'), 
                                      ('color', '|S10')])
    """
    self.triangles = triangs

	def write(self):
    """
    Write data
    """
		self._writeHeader()
    self._writeVertices()
    self._writeSurfaceTriangle()

  def _writeHeader(self):
    """
    Write header
    """

    dt2Type = {
      numpy.dtype('float32'): 'float',
      numpy.dtype('float64'): 'double',
      numpy.dtype('int32'): 'int',
      numpy.dtype('uint8'): 'uchar',
    }

    # global metadata
    print >> self.fileHandle, 'ply'
    print >> self.fileHandle, 'format ascii 1.0'
    userName = os.environ.get('USER', '')
    print >> self.fileHandle, 'comment author: {0}'.format(userName)
    print >> self.fileHandle, 'comment date: {0}'.format(time.asctime())

    # vertices
    numVerts = len(self.vertices)
    print >> self.fileHandle, 'element vertex {0}'.format(numVerts)
    for fld in self.vertices.dtype.fields.items():
      print >> self.fileHandle, 
               'property {0} {1}'.format(fld[0], dt2Type[fld[1][0]])

    # surface triangles
    numTriangs = len(self.triangles)
    print >> self.fileHandle, 'element face {0}'.format(numTriangs)
    print >> self.fileHandle, 'property list uchar int vertex_index'
    for fld in self.triangles.dtype.fields.items():
      name = fld[0]
      # skip since this is part of the above list
      if name == 'ia' or name == 'ib' or name == 'ic':
        continue
      typ = dt2Type[fld[1][0]]
      print >> self.fileHandle, 'property {0} {1}'.format(name, typ)

    # close
    print >> self.fileHandle, 'end_header'

  def _writeVertices(self):

    pass

  def _writeTriangles(self):

    pass

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
  triangs['ib'][1] = 1
  triangs['ic'][2] = 2

  pw = PLYWriter('test.ply')
  pw.setVertices(verts)
  pw.setTriangles(triangs)
  pw.write()

if __name__ == '__main__':
  test()
