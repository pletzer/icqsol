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
    self.vertexProps = []
    self.surfaceTriangles = []
    self.surfaceTriangleProps = []

  def setVertices(self, verts, properties = []):
    """
    Set vertex data
    @param verts array of x, y, z vertices
    """
    pass

  def setTriangles(self, triangles, properties = []):
    pass


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

    # global metadata
    print >> self.fileHandle, 'ply'
    print >> self.fileHandle, 'format ascii 1.0'
    print >> self.fileHandle, 'comment author: {0}'.format(os.environ['USER'])
    print >> self.fileHandle, 'comment date: {0}'.format(time.asctime())

    # vertices
    numVerts = self.vertices.shape[0]
    print >> self.fileHandle, 'element vertex {0}'.format(numVerts)
    print >> self.fileHandle, 'property float x'
    print >> self.fileHandle, 'property float y'
    print >> self.fileHandle, 'property float z'
    for i in range(len(self.vertexProps)):
      pname = self.vertProps[i]['name']
      ptype = self.vertProps[i]['type']
      print >> self.fileHandle, 'property {0} {1}'.format(ptype, pname)

    # surface triangles
    numSurfaceTriangles = self.surfaceTriangles.shape[0]
    print >> self.fileHandle, 'element face {0}'.format(numSurfaceTriangles)
    print >> self.fileHandle, 'property list uchar int vertex_index'
    for i in range(len(self.surfaceTriangleProps)):
      pame = self.surfaceTraingleProps[i]['name']
      ptype = self.surfaceTraingleProps[i]['type']
      print >> self.fileHandle, 'property {0} {1}'.format(ptype, pname)

    # close
    print >> self.fileHandle, 'end_header'

  def _writeVertices(self):

    pass

