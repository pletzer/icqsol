# standard python modules
import re

# extensions
import vtk
import numpy

from icqBaseShape import BaseShape

class PrimitiveShape(BaseShape):

  def __init__(self,):
    """
    Constructor
    """
    BaseShape.__init__(self,)

    self.surfaceFuncs = []
    self.evalFunc = None

  def setSurfaceFunctions(self, funcs):
    """
    Set the surface functions
    @param funcs a list of 3-tuples, [(x(u,v), y(u,v), z(u,v)), ...] 
           with 0 <= u, v <= 1
    """
    self.surfaceFuncs = funcs

  def setEvaluateFunction(self, func):
    """
    Set the function of (x, y, z) that returns > 0 inside and 0 outside
    """
    self.evalFunc = func

  def evaluate(self, pts):
    """
    Evaluate the inside/outside function 
    @param pts array of points
    """
    return self.evalFunc(pts)

  def computeSurfaceMeshes(self, maxTriArea):
    """
    Compute the surface meshes
    @param maxTriArea maximum triangle area
    """
    
    if self.surfaceMeshes:
      return

    points = []
    pointCount = 0
    self.surfaceMeshes = []

    # iterate over the faces
    for faceId in range(len(self.surfaceFuncs)):

      # estimate for the u, v griod resolution
      n = self._getSurfaceNumberOfCells(faceId, maxTriArea)
      nu, nv = n, n
      nu1, nv1 = nu + 1, nv + 1

      # compute the nodes
      xx, yy, zz = self._getSurfacePointArrays(faceId, nu, nv)
      numPoints = len(xx.flat)
      pts = numpy.zeros( (numPoints, 3), numpy.float64 )
      pts[:, 0], pts[:, 1], pts[:, 2] = xx.flat, yy.flat, zz.flat

      # store the vertices
      points.append(pts)

      numQuadCells = nu * nv
      jj = numpy.outer( range(nv1), numpy.ones((nu1,), numpy.int) )
      ii = numpy.outer( numpy.ones((nv1,), numpy.int), range(nu1) )

      # big index
      bigII = nu1*jj + ii

      connect = numpy.zeros( (2*numQuadCells, 3), numpy.int )

      # compute the connectivity for the lower triangles     
      connect[:numQuadCells, 0] = numpy.ravel(bigII[:-1, :-1])
      connect[:numQuadCells, 1] = numpy.ravel(bigII[:-1, 1:])
      connect[:numQuadCells, 2] = numpy.ravel(bigII[1:, 1:])

      # compute the connectivity for the upper triangles     
      connect[numQuadCells:, 0] = numpy.ravel(bigII[:-1, :-1])
      connect[numQuadCells:, 1] = numpy.ravel(bigII[1:, 1:])
      connect[numQuadCells:, 2] = numpy.ravel(bigII[1:, :-1])

      # increment vertex index with the number of points so far
      connect += pointCount

      # store
      self.surfaceMeshes.append(connect)

      # update the number of points
      pointCount += numPoints

    # pointCount is the total number of points
    self.points = numpy.zeros( (pointCount, 3), numpy.float64 )

    # copy the points into a single array. Note that there will be 
    # duplicate vertices
    iBeg = 0
    for faceId in range(len(self.surfaceFuncs)):
      numPoints = points[faceId].shape[0]
      iEnd = iBeg + numPoints
      self.points[iBeg:iEnd, :] = points[faceId]
      iBeg += numPoints

  def _getSurfacePointArrays(self, faceId, nu, nv):
    """
    Get the vertex arrays
    @param faceId face face Id
    @param nu number of u cells
    @param nv number of v cells
    @return xx, yy, zz coordinates
    """

    du, dv = 1.0/float(nu), 1.0/float(nv)
    nu1, nv1 = nu + 1, nv + 1

    # 1D arrays
    u = numpy.array([du*i for i in range(nu1)])
    v = numpy.array([dv*j for j in range(nv1)])

    # 2D arrays
    uu = numpy.outer(u, numpy.ones((nv1,), numpy.float64))
    vv = numpy.outer(numpy.ones((nu1,), numpy.float64), v)

    # compute the vertices
    xx = self.surfaceFuncs[faceId][0](uu, vv)
    yy = self.surfaceFuncs[faceId][1](uu, vv)
    zz = self.surfaceFuncs[faceId][2](uu, vv)

    return xx, yy, zz

  def _getSurfaceCellAreaVectors(self, xx, yy, zz):
    """
    Get cell area vectors
    @param xx x vertices (2D array)
    @param yy y vertices (2D array)
    @param zz z vertices (2D array)
    @return the components of the cell areas for each cell
    """

    dxu = xx[1:, :-1].copy()
    dxv = xx[:-1, 1:].copy()
    dyu = yy[1:, :-1].copy()
    dyv = yy[:-1, 1:].copy()
    dzu = zz[1:, :-1].copy()
    dzv = zz[:-1, 1:].copy()

    dxu -= xx[:-1, :-1]
    dxv -= xx[:-1, :-1]
    dyu -= yy[:-1, :-1]
    dyv -= yy[:-1, :-1]
    dzu -= zz[:-1, :-1]
    dzv -= zz[:-1, :-1]

    ax = dyu*dzv - dyv*dzu
    ay = dzu*dxv - dzv*dxu
    az = dxu*dyv - dxv*dyu
    
    return ax, ay, az
   
  def _getSurfaceNumberOfCells(self, faceId, maxTriArea):
    """
    Estimate the number of nu and nv cells from supplied maximum triangle area
    @param faceId face Id
    @param maxTriArea maximum triangle area
    """

    # should not need a lot of resolution to estimate
    n = 5
    # assume same resolution in u and v
    nu, nv = n, n
    # get the grid
    xx, yy, zz = self._getSurfacePointArrays(faceId, nu, nv)
    # get the cell areas
    ax, ay, az = self._getSurfaceCellAreaVectors(xx, yy, zz)
    # mac cell area
    maxHalfArea = 0.5*max( numpy.sqrt(ax*ax + ay*ay + az*az).flat )
    # estimate the number of cells
    print '*** maxTriArea = {0} maxHalfArea = {1}'.format(maxTriArea, maxHalfArea)
    nOpt = int(numpy.sqrt(maxHalfArea/maxTriArea) * n + 0.5)

    return max(nOpt, 1)

################################################################################
def test():

  # create a cylinder terminated by two half-spheres
  radius = 1.0
  radius2 = radius**2
  length = 0.5

  #
  # define the surface of the cylinder
  #
  def xCylSurf(u, v):
    return radius*numpy.cos(2*numpy.pi*u)

  def yCylSurf(u, v):
    return length*v

  def zCylSurf(u, v):
    return radius*numpy.sin(2*numpy.pi*u)

  #
  # define the low end surface
  #
  def xLoSurf(u, v):
    return -radius*numpy.sin(numpy.pi*u)*numpy.cos(numpy.pi*v)

  def yLoSurf(u, v):
    return -radius*numpy.sin(numpy.pi*u)*numpy.sin(numpy.pi*v)

  def zLoSurf(u, v):
    return radius*numpy.cos(numpy.pi*u)

  #
  # define the high end surface
  #
  def xHiSurf(u, v):
    return -radius*numpy.sin(numpy.pi*u)*numpy.cos(numpy.pi*v)

  def yHiSurf(u, v):
    return length + radius*numpy.sin(numpy.pi*u)*numpy.sin(numpy.pi*v)

  def zHiSurf(u, v):
    return radius*numpy.cos(numpy.pi*u)

  def objectVolume(x, y, z):
    rho2 = x**2 + z**2
    r12 = rho2 + y**2
    r22 = rho2 + (y-length)**2
    res = y < length
    res &= y > 0
    res |= (r12 < radius2)
    res |= (r22 < radius2)
    return res

  s = PrimitiveShape()
  s.setSurfaceFunctions([(xCylSurf, yCylSurf, zCylSurf), 
                         (xLoSurf, yLoSurf, zLoSurf),
                         (xHiSurf, yHiSurf, zHiSurf)])  

  print 'number of cells for face 0: ', s._getSurfaceNumberOfCells(0, 0.1)
  print 'number of cells for face 1: ', s._getSurfaceNumberOfCells(1, 0.1)
  print 'number of cells for face 2: ', s._getSurfaceNumberOfCells(2, 0.1)

  s.computeSurfaceMeshes(0.1)
  s.computeSurfaceNormals()

  s.save('testPrimitiveShape.vtk')

if __name__ == '__main__': test()
