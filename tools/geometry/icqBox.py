#!/usr/bin/env python

import vtk
import numpy

from icqShape import Shape
from icqUniformGrid import uniformGrid2D, uniformIndexGrid2D

class Box(Shape):

  def __init__(self, origin, lengths,
               n_x=2, n_y=2, n_z=2):
    """
    Constructor
    @param radius radius
    @param origin/low end of the box
    @param lengths lengths in x, y, and z
    @param n_x number of x cells
    @param n_y number of y cells
    @param n_z number of z cells
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
    
    # create side meshes
    for side in [(-1,0,0), (1,0,0), (0,-1,0), (0,1,0), (0,0,-1), (0,0,1)]:
        points, connectivity = self._createSideMesh(origin, lengths,
                                                    side=side,
                                                    n_x=n_x, n_y=n_y, n_z=n_z)
        pointList.append(points)
        connectivityList.append(connectivity)
    
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

  def _getDirection(self, vector):
    """
    Return the first index where vector is non-zero
    @return index
    """
    index = 0
    while index < 3 and vector[index] == 0:
        index += 1
    return index

  def _createSideMesh(self, origin, lengths, side, n_x, n_y, n_z):
    """
    Create side mesh
    @param origin origin of the box
    @param lengths lengths of the box in x, y, and z
    @param side normal vector of the side
    @param n_x number of x cells
    @param n_y number of y cells
    @param n_z number of z cells
    @return points, connectivity
    """

    # find the direction aligned to the normal vector 'side'
    index = self._getDirection(side)

    eVec = numpy.zeros((3,))
    indexU = (index - 1) % 3
    indexV = (indexU - 1) % 3
    eVec[indexV] = 1

    # the two directions tangential to the face, uHat cross vHat points out
    uHat = numpy.cross(eVec, side)
    vHat = numpy.cross(side, uHat)

    # vHat is always pointing in the positive direction but uHat may 
    # be pointing down or up, depending on the side
    lo = numpy.copy(origin)
    if side[index] == 1:
        lo[index] += lengths[index]
    if uHat[indexU] == -1:
        lo[indexU] += lengths[indexU]

    ns = (n_x, n_y, n_z)
    n_u, n_v = ns[indexU], ns[indexV]

    uu, vv = uniformGrid2D([0., 0.], [1., 1.], (n_u, n_v))

    xx = (uu*uHat[0] + vv*vHat[0]) * lengths[0]
    yy = (uu*uHat[1] + vv*vHat[1]) * lengths[1]
    zz = (uu*uHat[2] + vv*vHat[2]) * lengths[2]
    xx += lo[0]
    yy += lo[1]
    zz += lo[2]

    numPoints = len(xx.flat)
    points = numpy.zeros( (numPoints, 3), numpy.float64 )
    points[:, 0] = xx.flat
    points[:, 1] = yy.flat
    points[:, 2] = zz.flat
      
    # cell connectivity
    numQuadCells = n_u * n_v
    numCells = 2 * numQuadCells
    iiU, jjV = uniformIndexGrid2D((0, 0), (n_u, n_v))
    bigII = iiU*(n_v + 1) + jjV
      
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

  box = Box(origin=(0., 0., 0.), lengths=(0.5, 1., 2.), 
                 n_x=2, n_y=3, n_z=4)
  box.save('box.vtk', file_format='vtk', file_type='ascii')
  box.show()

if __name__ == '__main__':
  test()

