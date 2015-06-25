# standard python modules
import re

# extensions
import vtk
import numpy
from numpy import sin, cos, pi

from icqBaseShape import BaseShape

class PrimitiveShape(BaseShape):

  def __init__(self, angle=None, length=None, radius=None, origin=None, end_point=None ):
    """
    Constructor
    """
    BaseShape.__init__(self)
    self.angle = angle
    self.length = length
    self.radius = radius
    self.origin = origin
    self.end_point = end_point
    self.surfaceFuncs = []
    self.evalFunc = None
    self.surfaceFuncs = []
    self.evalFunc = None

  @property
  def origin_x(self):
    if self.origin is not None:
        return self.origin[ 0 ]
    return None

  @property
  def origin_y(self):
    if self.origin is not None:
      return self.origin[ 1 ]
    return None

  @property
  def origin_z(self):
    if self.origin is not None:
      return self.origin[ 2 ]
    return None

  @property
  def end_x(self):
    if self.end_point is not None:
      return self.end_point[ 0 ]
    return None

  @property
  def end_y(self):
    if self.end_point is not None:
      return self.end_point[ 1 ]
    return None

  @property
  def end_z(self):
    if self.end_point is not None:
      return self.end_point[ 2 ]
    return None

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
    nOpt = int(numpy.sqrt(maxHalfArea/maxTriArea) * n + 0.5)

    return max(nOpt, 1)


class Box( PrimitiveShape ):

    def __init__( self, origin, end_point ):
        PrimitiveShape.__init__( self, origin=origin, end_point=end_point )
        self.loBound = numpy.array( [ self.origin_x, self.origin_y, self.origin_z ] )
        self.hiBound = numpy.array( [ self.end_x, self.end_y, self.end_z ] )
        self.deltas = self.hiBound - self.loBound
        # Define the six faces of the box.
        self.surfaceFuncs = [ ( lambda u, v: self.loBound[0] * numpy.ones( u.shape ),
                                lambda u, v: self.loBound[1] + v * self.deltas[1],
                                lambda u, v: self.loBound[2] + u * self.deltas[2] ),
                              ( lambda u, v: self.hiBound[0] * numpy.ones( u.shape ),
                                lambda u, v: self.loBound[1] + u * self.deltas[1],
                                lambda u, v: self.loBound[2] + v * self.deltas[2] ),
                              ( lambda u, v: self.loBound[0] + u * self.deltas[0],
                                lambda u, v: self.loBound[1] * numpy.ones( u.shape ),
                                lambda u, v: self.loBound[2] + v * self.deltas[2] ),
                              ( lambda u, v: self.loBound[0] + v * self.deltas[0],
                                lambda u, v: self.hiBound[1] * numpy.ones( u.shape ),
                                lambda u, v: self.loBound[2] + u * self.deltas[2]),
                              ( lambda u, v: self.loBound[0] + v * self.deltas[0],
                                lambda u, v: self.loBound[1] + u * self.deltas[1],
                                lambda u, v: self.loBound[2] * numpy.ones( u.shape ) ),
                              ( lambda u, v: self.loBound[0] + u * self.deltas[0],
                                lambda u, v: self.loBound[1] + v * self.deltas[1],
                                lambda u, v: self.hiBound[2] * numpy.ones( u.shape ) ) ]
    def evalFunc( pts ):
        x, y, z = pts[ :, 0 ], pts[ :, 1 ], pts[ :, 2 ]
        res = x > self.loBound[0]
        res &= x < self.hiBound[0]
        res &= y > self.loBound[1]
        res &= y < self.hiBound[1]
        res &= z > self.loBound[2]
        res &= x < self.hiBound[2]
        return res


class Cone( PrimitiveShape ):

    def __init__( self, length, radius, origin ):
        PrimitiveShape.__init__( self, length=length, radius=radius, origin=origin )
        self.surfaceFuncs = [ ( lambda u, v: u * self.radius * cos( 2 * pi * v ) + self.origin_x,
                                lambda u, v: u * self.radius * sin( 2 * pi * v ) + self.origin_y,
                                lambda u, v: u * self.length + self.origin_z ),
                              ( lambda u, v: v * self.radius * cos( 2 * pi * u ) + self.origin_x,
                                lambda u, v: v * self.radius * sin( 2 * pi * u ) + self.origin_y,
                                lambda u, v: ( self.length + self.origin_z ) * numpy.ones( u.shape ) ) ]
        self.radiusSq = self.radius**2

    def evalFunc( pts ):
        xNorm = pts[ :, 0 ] - self.origin_x
        yNorm = pts[ :, 1 ] - self.origin_y
        zNorm = pts[ :, 2 ] - self.origin_z
        res = zNorm > 0
        res &= zNorm < self.length
        u = zNorm / self.length
        uRadius = u * self.radius
        res &= uRadius * uRadius - xNorm * xNorm - yNorm * yNorm > 0
        return res


class Cylinder( PrimitiveShape ):

    def __init__( self, length, radius, origin ):
        PrimitiveShape.__init__( self, length=length, radius=radius, origin=origin )
        self.surfaceFuncs = [ ( lambda u, v: self.radius * cos( 2 * pi * u ) + self.origin_x, 
                                lambda u, v: self.radius * sin( 2 * pi * u ) + self.origin_y, 
                                lambda u, v: self.length * ( v-0.5 ) + self.origin_z ),
                              ( lambda u, v: v * self.radius * cos( 2 * pi * u ) + self.origin_x,
                                lambda u, v: v * self.radius * sin( 2 * pi * u ) + self.origin_y,
                                lambda u, v: ( self.origin_z - 0.5 * self.length ) * numpy.ones( u.shape ) ),
                              ( lambda u, v: u * self.radius * cos( 2 * pi * v ) + self.origin_x,
                                lambda u, v: u * self.radius * sin( 2 * pi * v ) + self.origin_y,
                                lambda u, v: ( self.origin_z + 0.5 * self.length ) * numpy.ones( u.shape ) ) ]
        self.radiusSq = self.radius**2

    def evalFunc( pts ):
        res = ( self.radiusSq - ( pts[ :,0 ] - self.origin_x )**2 - ( pts[ :,1 ] - self.origin_y )**2 ) > 0
        zNorm = pts[ :,2 ] - self.origin_z
        res &= zNorm > -0.5 * self.length
        res &= zNorm < 0.5 * self.length
        return res


class Sphere( PrimitiveShape ):

    def __init__( self, radius, origin ):
        PrimitiveShape.__init__( self, radius=radius, origin=origin )

        self.surfaceFuncs = [ ( lambda u, v: self.radius * sin( pi * u ) * cos( 2 * pi * v ) + self.origin_x,
                                lambda u, v: self.radius * sin( pi * u ) * sin( 2 * pi * v ) + self.origin_y,
                                lambda u, v: self.radius * cos( pi * u ) + self.origin_z ) ]
        self.radiusSq = self.radius**2

    def evalFunc( pts ):
        xNorm = pts[ :, 0 ] - self.origin_x
        yNorm = pts[ :, 1 ] - self.origin_y
        zNorm = pts[ :, 2 ] - self.origin_z
        return self.radiusSq - xNorm**2 - yNorm**2 - zNorm**2 > 0

################################################################################
def test():

  # create a cylinder terminated by two half-spheres
  print "Testing creation of a cylinder terminated by 2 half spheres..."
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

  print "Testing creation of a cylinder..."
  radius = 1.0
  length = 0.5
  origin_tup = ( 0.0, 0.0, 0.0 )
  
  s = Cylinder( length, radius, origin_tup  )
  print 'number of cells for face 0: ', s._getSurfaceNumberOfCells(0, 0.1)
  print 'number of cells for face 1: ', s._getSurfaceNumberOfCells(1, 0.1)
  print 'number of cells for face 2: ', s._getSurfaceNumberOfCells(2, 0.1)
  
if __name__ == '__main__': test()
