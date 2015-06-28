#!/usr/bin/env python

import vtk
import numpy

from icqShape import Shape
from icqUniformGrid import uniformGrid2D, uniformIndexGrid2D

class Cylinder(Shape):

  def __init__(self, radius, origin, length, n_theta = 32, n_rho=10, n_z=2):
    """
    Constructor
    @param radius radius
    @param origin center of low end disk
    @param length length of the cylinder in the z direction
    @param n_theta number of theta cells
    @param n_rho number of radial cells
    @param n_z number of z cells
    """

    Shape.__init__(self)
    
    self.pointList = []
    self.connectivityList = []

    # create outer cylinder
    # u ~ theta, v ~ z
    uu, zz = uniformGrid2D([0., 0.], [2*numpy.pi, length], (n_theta, n_z))
    xx = radius*numpy.cos(uu)
    yy = radius*numpy.sin(uu)
    xx += origin[0]
    yy += origin[1]
    zz += origin[2]
    numPoints = len(xx.flat)
    pointsOuter = numpy.zeros( (numPoints, 3), numpy.float32 )
    pointsOuter[:, 0] = xx.flat
    pointsOuter[:, 1] = yy.flat
    pointsOuter[:, 2] = zz.flat
    
    # cell connectivity
    numQuadCells = n_theta * n_z
    numCells = 2 * numQuadCells
    iiTheta, jjZ = uniformIndexGrid2D((0, 0), (n_theta, n_z))
    bigII = iiTheta*(n_z + 1) + jjZ
    
    connectivityOuter = numpy.zeros( (numCells, 4), numpy.int )
    
    # lower triangles
    connectivityOuter[:numQuadCells, 1] = numpy.ravel(bigII[:-1, :-1])
    connectivityOuter[:numQuadCells, 2] = numpy.ravel(bigII[1:, :-1])
    connectivityOuter[:numQuadCells, 3] = numpy.ravel(bigII[1:, 1:])
    # upper triangles
    connectivityOuter[numQuadCells:, 1] = numpy.ravel(bigII[:-1, :-1])
    connectivityOuter[numQuadCells:, 2] = numpy.ravel(bigII[1:, 1:])
    connectivityOuter[numQuadCells:, 3] = numpy.ravel(bigII[:-1, 1:])
    connectivityOuter[:, 0] = 3 # triangles
    
    # make sure to save the point list and connectivity list
    self.pointList.append(pointsOuter)
    self.connectivityList.append(connectivityOuter)
    
    # set the points and connectivity array
    self.pointArray.SetNumberOfTuples(numPoints)
    self.pointArray.SetVoidArray(pointsOuter, numPoints*3, 1)
    
    self.cellInds.SetNumberOfTuples(numCells)
    self.cellInds.SetVoidArray(connectivityOuter, numCells*4, 1)
    
    self.cells.SetCells(numCells, self.cellInds)

################################################################################
def test():

  cyl = Cylinder(radius=1.0, origin=(0., 0., 0.), length=1., n_theta=4, n_z=1)
  cyl.save('cyl.vtk', file_format='vtk', file_type='ascii')
  cyl.show()

if __name__ == '__main__':
  test()

