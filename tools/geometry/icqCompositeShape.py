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
  
    # compute cell normals
    self.computeSurfaceNormals()

    #
    # find all the points inside the volume
    #

    # move point out by a small amount to avoid floating point comparision
    # issues
    eps = 1.e-6
    ONE_THIRD = 1./3.

    print '*** self.surfaceMeshes = ', self.surfaceMeshes
    numFaces = len(self.surfaceMeshes)
    insidePointList = []
    totalNumPoints = 0
    loBound = numpy.array([float('inf')] * 3)
    hiBound = numpy.array([-float('inf')] * 3)
    for shp in self.argShapes:
      numFaces = len(shp.surfaceMeshes)
      for iFace in range(numFaces):
        print '*** iFace = ', iFace
        face = shp.surfaceMeshes[iFace]
        shp.computeSurfaceNormals()
        normals = shp.surfaceNormals[iFace]
        print '*** normals = ', normals
        print '*** face = ', face
        print '*** points = ', shp.points
        verts = shp.points[face] 
        displVerts = verts.copy()
        displVerts[:, 0, :] += eps*normals
        displVerts[:, 1, :] += eps*normals
        displVerts[:, 2, :] += eps*normals
        inPts = verts[self.evaluate(displVerts) > 0]
        totalNumPoints += inPts.shape[0]
        insidePointList.append(inPts)
        xs, ys, zs = inPts[:, 0], inPts[:, 1], inPts[:, 2]
        print '*** xs, ys, zs = ', xs, ys, zs
        loBound[0] = min(loBound[0], xs.min())
        loBound[1] = min(loBound[1], ys.min())
        loBound[2] = min(loBound[2], zs.min())
        hiBound[0] = max(hiBound[0], xs.max())
        hiBound[1] = max(hiBound[1], ys.max())
        hiBound[2] = max(hiBound[2], zs.max())

    # add more points by sampling the volume containing the boundary surfaces
    cellLength = numpy.sqrt(2*maxTriArea)
    ns = numpy.array([ int( (hiBound[i] - loBound[i])/cellLength + 0.5) for i in range(3) ])
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

    # create a single list of points
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

    # remove the cells that are outside the closed object
    pts = ugrid.GetPoints()
    cells = ugrid.GetCells()
    ptIds = vtk.vtkIdList()
    print '*** cells = ', cells
    numCells = cells.GetNumberOfCells()
    volumeMesh = numpy.zeros( (numCells, 4), numpy.int )
    centroidPoints = numpy.zeros( (numCells, 3,), numpy.float64 )
    for i in range(numCells):
      cell = cells.GetCell(i, ptIds)
      i0, i1, i2, i3 = ptIds.GetId(0), ptIds.GetId(1), ptIds.GetId(2), ptIds.GetId(3)
      p0, p1, p2, p3 = pts.GetPoint(i0), pts.GetPoint(i1), pts.GetPoint(i2), pts.GetPoint(i3)
      centroidPoints[i, :] += numpy.array(p0)
      centroidPoints[i, :] += numpy.array(p1)
      centroidPoints[i, :] += numpy.array(p2)
      centroidPoints[i, :] += numpy.array(p3)
      centroidPoints[i, :] *= 0.25
      volumeMesh[i, :] = i0, i1, i2, i3

    validCells = self.evaluate(centroidPoints)
    print '*** validCells = ', validCells

    ugrid2 = vtk.vtkUnstructuredGrid()
    ugrid2.SetPoints(pts)
    numCells2 = int(validCells.sum())
    ugrid2.Allocate(numCells2, 1)
    for i in range(numCells):
      if validCells[i]:
        cell = cells.GetCell(i, ptIds)
        ugrid2.InsertNextCell(vtk.VTK_TETRA, 4, volumeMesh[i, :])

    # apply filter to extract boundary cell faces

    print '*** ugrid2 = ', ugrid2

    surf = vtk.vtkGeometryFilter()
    if vtk.VTK_MAJOR_VERSION >= 6:
      surf.SetInputData(ugrid2)
    else:
      surf.SetInput(ugrid2)
    surf.Update()

    print '*** surf = ', surf

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(surf.GetOutputPort())

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    ren = vtk.vtkRenderer()
    win = vtk.vtkRenderWindow()
    iren = vtk.vtkRenderWindowInteractor()

    ren.AddActor(actor)
    win.AddRenderer(ren)
    iren.SetRenderWindow(win)
    win.SetSize(500, 500)

    iren.Initialize()
    win.Render()
    iren.Start()


    # get the boundary cells
    # copy all the surface meshes into a single connectivity array
    numCells = surf.GetNumberOfCells()
    cellArr = numpy.zeros( (numCells, 3), numpy.int )
    for i in range(numCells):
      cell = surf.GetCell(i)
      ptIds = cell.GetPointIds()
      npts = ptIds.GetNumberOfIds()
      pInds = [0 for j in range(npts)]
      for j in range(npts):
        pInds[j] = int(ptIds.GetId(j))
      cellArr[i, :] = pInds

    # single surface in the case of a composite object
    self.surfaceMeshes = [cellArr]

  def computeSurfaceMeshes2(self, maxTriArea):
    """
    Compute the surface meshes
    @param maxTriArea maximum triangle area
    """

    # get the bounds
    loBound, hiBound = self._getBounds(maxTriArea)

    # get the optimal number of cells in each direction
    nx, ny, nz = self._getUniformGridResolution(loBound, hiBound, maxTriArea)
    nx1, ny1, nz1 = nx + 1, ny + 1, nz + 1

    x = loBound[0] + (hiBound[0] - loBound[0])*numpy.array([i for i in range(nx1)])/nx
    y = loBound[1] + (hiBound[1] - loBound[1])*numpy.array([j for j in range(ny1)])/ny
    z = loBound[2] + (hiBound[2] - loBound[2])*numpy.array([i for k in range(nz1)])/nz
    xxx = numpy.outer(numpy.ones((nz1, ny1)), x)
    yyy = numpy.outer(numpy.outer(numpy.ones((nz1,)), y), numpy.ones((nx1,)))
    zzz = numpy.outer(z, numpy.ones((ny1, nx1)))
    verts = numpy.zeros( (nz1*ny1*nx1, 3), numpy.float64 )
    verts[:, 0] = xxx.flat
    verts[:, 1] = yyy.flat
    verts[:, 2] = zzz.flat

    insideVals = self.evaluate(verts)

    # build the VTK contour pipeline
    xyz = vtk.vtkDoubleArray()
    pts = vtk.vtkPoints()
    data = vtk.vtkDoubleArray()
    grd = vtk.vtkStructuredGrid()
    cont = vtk.vtkContourFilter()

    numPoints = nx1*ny1*nz1
    xyz.SetNumberOfTuples(numPoints)
    xyz.SetNumberOfComponents(3)
    xyz.SetVoidArray(verts, numPoints*3, 1)

    pts.SetData(xyz)

    data.SetNumberOfTuples(numPoints)
    data.SetNumberOfComponents(1)
    data.SetVoidArray(insideVals, numPoints*3, 1)

    grd.SetDimensions(nz1, ny1, nx1)
    grd.SetPoints(pts)
    grd.GetPointData().SetScalars(data)

    if vtk.VTK_MAJOR_VERSION >= 6:
      cont.SetInputData(grd)
    else:
      cont.SetInput(grd)
    cont.SetValue(0, 0.0) # or should it be 0.5?
    cont.Update()

    surfPolyData = cont.GetOutput()
    
    # gather the boundary points
    surfPoints = surfPolyData.GetPoints()
    numSurfPoints = surfPoints.GetNumberOfPoints()
    self.points = numpy.zeros( (numSurfPoints, 3), numpy.float64 )
    for i in range(numSurfPoints):
      self.points[i, :] = surfPoints.GetPoint(i)

    # gather the boundary triangles
    numSurfCells = surfPolyData.GetNumberOfCells()
    cellArr = numpy.zeros( (numSurfCells, 3), numpy.int )
    for i in range(numSurfCells):
      cell = surfPolyData.GetCell(i)
      ptIds = cell.GetPointIds()
      npts = ptIds.GetNumberOfIds()
      pInds = [0 for j in range(npts)]
      for j in range(npts):
        pInds[j] = int(ptIds.GetId(j))
      cellArr[i, :] = pInds

    # single surface mesh in this case
    self.surfaceMeshes = [cellArr]


  def _getBounds(self, maxTriArea):
    """
    Compute the min/max corner points
    @param maxTriArea may need to compute the individual shape surface meshes and use maxTriArea
                      to do so
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
  radius1 = 0.5
  origin1 = numpy.array([0.1, 0.2, 0.3])
  surfaceFunctions1 = [(lambda u,v: radius1*sin(pi*u)*cos(2*pi*v) + origin1[0], 
                       lambda u,v: radius1*sin(pi*u)*sin(2*pi*v) + origin1[1], 
                       lambda u,v: radius1*cos(pi*u) + origin1[2])]
  def evalFunction1(pts):
    xNorm = pts[:, 0] - origin1[0]
    yNorm = pts[:, 1] - origin1[1]
    zNorm = pts[:, 2] - origin1[2]
    return radius1*radius1 - xNorm**2 - yNorm**2 - zNorm**2 > 0

  s1 = PrimitiveShape()
  s1.setSurfaceFunctions(surfaceFunctions1)
  s1.setEvaluateFunction(evalFunction1)
  s1.computeSurfaceMeshes(maxTriArea=0.1)

  # create second sphere
  radius2 = 0.5
  origin2 = numpy.array([0.9, 0.8, 0.6])
  surfaceFunctions2 = [(lambda u,v: radius2*sin(pi*u)*cos(2*pi*v) + origin2[0], 
                       lambda u,v: radius2*sin(pi*u)*sin(2*pi*v) + origin2[1], 
                       lambda u,v: radius2*cos(pi*u) + origin2[2])]
  def evalFunction2(pts):
    xNorm = pts[:, 0] - origin2[0]
    yNorm = pts[:, 1] - origin2[1]
    zNorm = pts[:, 2] - origin2[2]
    return radius2*radius2 - xNorm**2 - yNorm**2 - zNorm**2 > 0

  # create another sphere
  s2 = PrimitiveShape()
  s2.setSurfaceFunctions(surfaceFunctions2)
  s2.setEvaluateFunction(evalFunction2)
  s2.computeSurfaceMeshes(maxTriArea=0.1)

  # assemble
  cs = CompositeShape()
  cs.assemble('$0 + $1', (s1, s2))
  cs.computeSurfaceMeshes(maxTriArea=0.01)
  cs.save('testCompositeShape.vtk')


if __name__ == '__main__': test()
