# standard python modules
import re

# extensions
import vtk
import numpy

from icqBaseShape import BaseShape
from icqUniformGrid import uniformGrid
from icqTriangulation import Triangulation

class CompositeShape(BaseShape):

  def __init__(self,):
    """
    Constructor
    """
    BaseShape.__init__(self,)

    self.expression = ''
    self.argShapes = []

  def assemble(self, expression, shapes):
    """
    Assemble shapes into a composite object
    @param expression expression of the shape objects, $0 for shape 0, $1 for shape 1, etc. 
           operations include + (union), * (intersection), and - (removal)
    """
    self.expression = expression
    self.argShapes = shapes
    # replace $1 by self.argShape[0]
    for i in range(len(self.argShapes)):
      self.expression = re.sub('\${0}'.format(i), 
        'self.argShapes[{0}].evaluate(pts)'.format(i), self.expression)

  def evaluate(self, pts):
    """
    Evaluate the inside/outside function 
    @param pts array of points
    """
    return eval(self.expression)

  def computeSurfaceMeshes(self, maxTriArea):
    """
    Compute the surface meshes
    @param maxTriArea maximum triangle area
    """
    #
    # find all the points inside the volume
    #

    # move point outwards by a small amount to avoid floating
    # point comparison issues when determining whether a point is
    # in or out
    eps = 1.e-6

    numFaces = len(self.surfaceMeshes)
    insidePointList = []
    totalNumPoints = 0
    loBound = numpy.array([float('inf')] * 3)
    hiBound = numpy.array([-float('inf')] * 3)
    for shp in self.argShapes:
      shp.computeSurfaceMeshes(maxTriArea)
      shp.computeSurfaceNormals()
      verts = shp.getPoints()
      inPts = verts[ self.evaluate(verts) > 0 ]
      totalNumPoints += inPts.shape[0]
      insidePointList.append(inPts)
      xs, ys, zs = inPts[:, 0], inPts[:, 1], inPts[:, 2]
      loBound[0] = min(loBound[0], xs.min())
      loBound[1] = min(loBound[1], ys.min())
      loBound[2] = min(loBound[2], zs.min())
      hiBound[0] = max(hiBound[0], xs.max())
      hiBound[1] = max(hiBound[1], ys.max())
      hiBound[2] = max(hiBound[2], zs.max())

    # add more points by sampling the volume containing the boundary surfaces
    cellLength = numpy.sqrt(2*maxTriArea)
    ns = numpy.array([ max(1, int( (hiBound[i] - loBound[i])/cellLength + 0.5)) \
                      for i in range(3) ])
    xx, yy, zz = uniformGrid(loBound, hiBound, ns)
    npts = len(xx.flat)
    gridPts = numpy.zeros( (npts, 3), numpy.float64 )
    gridPts[:, 0] = xx.flat
    gridPts[:, 1] = yy.flat
    gridPts[:, 2] = zz.flat

    # only add the points inside
    gridPts = gridPts[ self.evaluate(gridPts) > 0 ]
    totalNumPoints += gridPts.shape[0]
    insidePointList.append(gridPts)

    # create a single list out of all these points
    insidePoints = numpy.zeros( (totalNumPoints, 3), numpy.float64 )
    iBeg = 0
    for pts in insidePointList:
      iEnd = iBeg + pts.shape[0]
      insidePoints[iBeg:iEnd, :] = pts
      iBeg = iEnd 

    # tessellate
    tri = Triangulation()
    tri.setInputPoints(insidePoints)
    tri.triangulate()
    ugrid = tri.delny.GetOutput()

    # remove the cells that are outside of the closed object
    pts = ugrid.GetPoints()
    cells = ugrid.GetCells()
    ptIds = vtk.vtkIdList()
    numCells = cells.GetNumberOfCells()
    centroidPoints = numpy.zeros( (numCells, 3,), numpy.float64 )
    cells.InitTraversal()
    for i in range(numCells):
      cell = cells.GetNextCell(ptIds)
      i0, i1, i2, i3 = ptIds.GetId(0), ptIds.GetId(1), ptIds.GetId(2), ptIds.GetId(3)
      p0, p1, p2, p3 = pts.GetPoint(i0), pts.GetPoint(i1), pts.GetPoint(i2), pts.GetPoint(i3)
      centroidPoints[i, :] += numpy.array(p0)
      centroidPoints[i, :] += numpy.array(p1)
      centroidPoints[i, :] += numpy.array(p2)
      centroidPoints[i, :] += numpy.array(p3)
      centroidPoints[i, :] *= 0.25

    validCells = self.evaluate(centroidPoints)

    # construct unstructured grid with cells that are fully
    # contained within the composite object
    ugrid2 = vtk.vtkUnstructuredGrid()
    ugrid2.SetPoints(pts)
    numCells2 = int(validCells.sum())
    ugrid2.Allocate(numCells2, 1)
    cells.InitTraversal()
    for i in range(numCells):
      cell = cells.GetNextCell(ptIds)
      if validCells[i]:
        ugrid2.InsertNextCell(vtk.VTK_TETRA, ptIds)

    # apply filter to extract boundary cell faces
    surf = vtk.vtkGeometryFilter()
    if vtk.VTK_MAJOR_VERSION >= 6:
      surf.SetInputData(ugrid2)
    else:
      surf.SetInput(ugrid2)
    surf.Update()

    # get the boundary cells
    # copy all the surface meshes into a single connectivity array
    # single surface in the case of a composite object
    pdata = surf.GetOutput()
    numCells = pdata.GetNumberOfPolys()
    cells = pdata.GetPolys()
    cells.InitTraversal()
    cellArr = numpy.zeros( (numCells, 3), numpy.int )
    for i in range(numCells):
      cell = cells.GetNextCell(ptIds)
      npts = ptIds.GetNumberOfIds()
      assert(npts == 3)
      cellArr[i, :] = ptIds.GetId(0), ptIds.GetId(1), ptIds.GetId(2)

    # set the points
    numPoints = pts.GetNumberOfPoints()
    self.points = numpy.zeros( (numPoints, 3), numpy.float64 )
    for i in range(numPoints):
      self.points[i, :] = pts.GetPoint(i)

    self.surfaceMeshes = [cellArr]

  def _getBounds(self, maxTriArea):
    """
    Compute the min/max corner points
    @param maxTriArea defines the resolution of the sampling grid
    """
    loBound = numpy.array([float('inf')] * 3)
    hiBound = numpy.array([-float('inf')] * 3)
    for shp in self.argShapes:
      points = shp.getPoints()
      if len(points) == 0:
        shp.computeSurfaceMeshes(maxTriArea=maxTriArea)
        points = shp.getPoints()
      loBound[0] = min(loBound[0], points[:, 0].min())
      loBound[1] = min(loBound[1], points[:, 1].min())
      loBound[2] = min(loBound[2], points[:, 2].min())
      hiBound[0] = max(hiBound[0], points[:, 0].min())
      hiBound[1] = max(hiBound[1], points[:, 1].min())
      hiBound[2] = max(hiBound[2], points[:, 2].min())
    return loBound, hiBound

  def _getUniformGridResolution(self, loBound, hiBound, maxTriArea):
    """
    Get the uniform grid's resolution
    @param loBound low corner of the grid
    @param hiBound high corner of the grid
    @param maxTriArea maximum surface triangle area
    """
    sr3 = numpy.sqrt(3.)
    nx = max(1, int(sr3*(hiBound[0] - loBound[0])**2/maxTriArea + 0.5))
    ny = max(1, int(sr3*(hiBound[1] - loBound[1])**2/maxTriArea + 0.5))
    nz = max(1, int(sr3*(hiBound[2] - loBound[2])**2/maxTriArea + 0.5))
    return nx, ny, nz



################################################################################
def test():

  from icqPrimitiveShape import PrimitiveShape
  from numpy import sin, cos, pi

  # create first sphere
  radius1 = 1.0
  origin1 = numpy.array([0., 0., 0.])
  surfaceFunctions1 = [(lambda u,v: radius1*sin(pi*u)*cos(2*pi*v) + origin1[0], 
                       lambda u,v: radius1*sin(pi*u)*sin(2*pi*v) + origin1[1], 
                       lambda u,v: radius1*cos(pi*u) + origin1[2])]
  def evalFunction1(pts):
    # make sure the boundary triangles are inside the shape
    eps = 1.e-6
    xNorm = pts[:, 0] - origin1[0]
    yNorm = pts[:, 1] - origin1[1]
    zNorm = pts[:, 2] - origin1[2]
    return radius1**2 + eps - xNorm**2 - yNorm**2 - zNorm**2 > 0

  s1 = PrimitiveShape()
  s1.setSurfaceFunctions(surfaceFunctions1)
  s1.setEvaluateFunction(evalFunction1)
  #s1.computeSurfaceMeshes(maxTriArea=0.1)

  # create second sphere
  radius2 = 1.0
  origin2 = numpy.array([2.0, 0.0, 0.0])
  surfaceFunctions2 = [(lambda u,v: radius2*sin(pi*u)*cos(2*pi*v) + origin2[0], 
                        lambda u,v: radius2*sin(pi*u)*sin(2*pi*v) + origin2[1],
                        lambda u,v: radius2*cos(pi*u) + origin2[2])]
  def evalFunction2(pts):
    # make sure the boundary triangles are inside the shape
    eps = 1.e-6
    xNorm = pts[:, 0] - origin2[0]
    yNorm = pts[:, 1] - origin2[1]
    zNorm = pts[:, 2] - origin2[2]
    return radius2**2 + eps - xNorm**2 - yNorm**2 - zNorm**2 > 0

  s2 = PrimitiveShape()
  s2.setSurfaceFunctions(surfaceFunctions2)
  s2.setEvaluateFunction(evalFunction2)
  s2.computeSurfaceMeshes(maxTriArea=0.1)
  s2.save('s2.vtk')

  # assemble
  cs = CompositeShape()
  #cs.assemble('$0 + $1', (s1, s2))
  cs.assemble('$0 + $1', (s1, s2,))
  cs.computeSurfaceMeshes(maxTriArea=0.05)
  cs.save('testCompositeShape.vtk')


if __name__ == '__main__': test()
