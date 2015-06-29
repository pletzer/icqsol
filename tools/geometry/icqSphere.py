#!/usr/bin/env python

import vtk
import numpy

from icqShape import Shape
from icqUniformGrid import uniformGrid2D, uniformIndexGrid2D

class Sphere(Shape):

  def __init__(self, radius, origin,
               n_theta=32, n_phi=16):
    """
    Constructor
    @param radius radius
    @param origin center of the sphere
    @param n_theta number of theta cells
    @param n_phi number of azimuthal cells
    """

    Shape.__init__(self)
    
    # data structure holding array of points
    self.pointArray = vtk.vtkDoubleArray()
    self.pointArray.SetNumberOfComponents(3)
    
    # data structire holding points
    self.points = vtk.vtkPoints()
    self.points.SetData(self.pointArray)
    
    # data structure holding cell indices
    self.cellInds = vtk.vtkIdTypeArray()
    
    # data structure holding cell connectivity data
    self.cells = vtk.vtkCellArray()
    
    # list of numpy arrays, one element for each surface mesh
    pointList = []
    connectivityList = []
    
    # create outer, side mesh
    pointsOuter, connectivityOuter = self._createSideMesh(origin, radius,
                                                          n_theta=n_theta,
                                                          n_phi=n_phi)
    pointList.append(pointsOuter)
    connectivityList.append(connectivityOuter)
    
    # concatenate all the meshes
    totNumPoints = 0
    totNumCells = 0
    numMeshes = len(pointList)
    for i in range(numMeshes):
      cl = connectivityList[i]
      # increment the connectivity
      cl[:, 1:] += totNumPoints
      numPoints = pointList[i].shape[0]
      numCells = cl.shape[0]
      totNumPoints += numPoints
      totNumCells += numCells
    
    self.allPoints = numpy.zeros( (totNumPoints, 3), numpy.float64 )
    self.allConnectivity = numpy.zeros( (totNumCells, 4), numpy.int )
    
    totNumPoints = 0
    totNumCells = 0
    for i in range(numMeshes):
      pl = pointList[i]
      cl = connectivityList[i]
      numPoints = pl.shape[0]
      numCells = cl.shape[0]
      self.allPoints[totNumPoints: totNumPoints + numPoints, :] = pl
      self.allConnectivity[totNumCells: totNumCells + numCells, :] = cl
      totNumPoints += numPoints
      totNumCells += numCells

    # set the points and connectivity array
    self.pointArray.SetNumberOfTuples(totNumPoints)
    self.pointArray.SetVoidArray(self.allPoints, totNumPoints*3, 1)
    
    self.cellInds.SetNumberOfTuples(totNumCells)
    self.cellInds.SetVoidArray(self.allConnectivity, totNumCells*4, 1)
    
    self.cells.SetCells(totNumCells, self.cellInds)

    self.surfPolyData.SetPoints(self.points)
    self.surfPolyData.SetPolys(self.cells)

  def _createSideMesh(self, origin, radius, n_theta, n_phi):
    """
    Create side mesh
    @param origin origin of the sphere
    @param radius radius
    @param n_theta number of poloidal cells
    @param n_phi number of azimuthal cells
    @return points, connectivity
    """
    # u ~ theta, v ~ z
    uu, vv = uniformGrid2D([0., -numpy.pi/2.], [2*numpy.pi, numpy.pi/2.],
                           (n_theta, n_phi))
    xx = radius*numpy.cos(vv)*numpy.cos(uu)
    yy = radius*numpy.cos(vv)*numpy.sin(uu)
    zz = radius*numpy.sin(vv)
    xx += origin[0]
    yy += origin[1]
    zz += origin[2]
    numPoints = len(xx.flat)
    points = numpy.zeros( (numPoints, 3), numpy.float64 )
    points[:, 0] = xx.flat
    points[:, 1] = yy.flat
    points[:, 2] = zz.flat
      
    # cell connectivity
    numQuadCells = n_theta * n_phi
    numCells = 2 * numQuadCells
    iiTheta, jjPhi = uniformIndexGrid2D((0, 0), (n_theta, n_phi))
    bigII = iiTheta*(n_phi + 1) + jjPhi
      
    connectivity = numpy.zeros( (numCells, 4), numpy.int )
      
    # lower triangles
    connectivity[:numQuadCells, 1] = numpy.ravel(bigII[:-1, :-1])
    connectivity[:numQuadCells, 2] = numpy.ravel(bigII[1:, :-1])
    connectivity[:numQuadCells, 3] = numpy.ravel(bigII[1:, 1:])
    
    # upper triangles
    connectivity[numQuadCells:, 1] = numpy.ravel(bigII[:-1, :-1])
    connectivity[numQuadCells:, 2] = numpy.ravel(bigII[1:, 1:])
    connectivity[numQuadCells:, 3] = numpy.ravel(bigII[:-1, 1:])
    connectivity[:, 0] = 3 # triangles

    return points, connectivity

################################################################################
def test():

  sph = Sphere(radius=1.0, origin=(0., 0., 0.),
                 n_theta=8, n_phi=4)
  sph.save('sph.vtk', file_format='vtk', file_type='ascii')
  sph.show()

if __name__ == '__main__':
  test()

