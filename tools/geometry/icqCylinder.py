#!/usr/bin/env python

import vtk
import numpy

from icqShape import Shape
from icqUniformGrid import uniformGrid2D, uniformIndexGrid2D

class Cylinder(Shape):

  def __init__(self, radius, origin, length,
               n_rho=10, n_theta = 32, n_z=2):
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

    # create the outer cylinder
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
    
    # create the low end disk
    # u ~ rho, v ~ theta
    uu, vv = uniformGrid2D([0., 0.], [radius, 2*numpy.pi], (n_rho, n_theta))
    xx = uu*numpy.cos(vv)
    yy = uu*numpy.sin(vv)
    zz = numpy.zeros( xx.shape, xx.dtype )
    xx += origin[0]
    yy += origin[1]
    zz += origin[2]
    numPoints = len(xx.flat)
    pointsLo = numpy.zeros( (numPoints, 3), numpy.float32 )
    pointsLo[:, 0] = xx.flat
    pointsLo[:, 1] = yy.flat
    pointsLo[:, 2] = zz.flat
    
    # cell connectivity
    numQuadCells = n_rho * n_theta
    numCells = 2 * numQuadCells
    iiRho, jjTheta = uniformIndexGrid2D((0, 0), (n_rho, n_theta))
    bigII = iiRho*(n_theta + 1) + jjTheta
    
    connectivityLo = numpy.zeros( (numCells, 4), numpy.int )
    
    # lower triangles
    connectivityLo[:numQuadCells, 1] = numpy.ravel(bigII[:-1, :-1])
    connectivityLo[:numQuadCells, 2] = numpy.ravel(bigII[1:, :-1])
    connectivityLo[:numQuadCells, 3] = numpy.ravel(bigII[1:, 1:])
    # upper triangles
    connectivityLo[numQuadCells:, 1] = numpy.ravel(bigII[:-1, :-1])
    connectivityLo[numQuadCells:, 2] = numpy.ravel(bigII[1:, 1:])
    connectivityLo[numQuadCells:, 3] = numpy.ravel(bigII[:-1, 1:])
    connectivityLo[:, 0] = 3 # triangles
    
    # make sure to save the point list and connectivity list
    self.pointList.append(pointsLo)
    self.connectivityList.append(connectivityLo)

    # concatenate all the meshes
    totNumPoints = 0
    totNumCells = 0
    numMeshes = len(self.pointList)
    for i in range(numMeshes):
      self.connectivityList[i][:, 1:] += totNumPoints
      numPoints = self.pointList[i].shape[0]
      numCells = self.connectivityList[i].shape[0]
      totNumPoints += numPoints
      totNumCells += numCells
    
    allPoints = numpy.zeros( (totNumPoints, 3), numpy.float32 )
    allConnectivity = numpy.zeros( (totNumCells, 4), numpy.int )
    
    totNumPoints = 0
    totNumCells = 0
    for i in range(numMeshes):
      numPoints = self.pointList[i].shape[0]
      numCells = self.connectivityList[i].shape[0]
      allPoints[totNumPoints: totNumPoints + numPoints, :] = \
            self.pointList[i]
      allConnectivity[totNumCells: totNumCells + numCells, :] = \
            self.connectivityList[i]
      totNumPoints += numPoints
      totNumCells += numCells

    # set the points and connectivity array
    self.pointArray.SetNumberOfTuples(totNumPoints)
    self.pointArray.SetVoidArray(allPoints, totNumPoints*3, 1)
    
    self.cellInds.SetNumberOfTuples(totNumCells)
    self.cellInds.SetVoidArray(allConnectivity, totNumCells*4, 1)
    
    self.cells.SetCells(totNumCells, self.cellInds)

################################################################################
def test():

  cyl = Cylinder(radius=1.0, origin=(0., 0., 0.), length=1.,
                 n_rho=2, n_theta=4, n_z=1)
  cyl.save('cyl.vtk', file_format='vtk', file_type='ascii')
  #cyl.show()

if __name__ == '__main__':
  test()

