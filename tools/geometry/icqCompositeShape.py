# standard python modules
import re

# extensions
import vtk
import numpy

from icqBaseShape import BaseShape

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

    # get the bounds
    loBound, hiBound = self._getBounds()

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

    data.SetNumberOfTuples(numPoints)
    data.SetNumberOfComponents(1)
    data.SetVoidArray(insideVals, numPoints*3, 1)

    grd.SetDimensions(nz1, ny1, nx1)
    grd.SetPoints(pts)
    grd.GetPointData().SetScalars(data)

    cont.SetInputConnection(grd.GetOutputPort())
    cont.SetValue(0, 0.0) # or should it be 0.5?
    cont.Update()

    pdata = cont.GetOutput()
    print pdata

  def _getBounds(self,):
    """
    Compute the min/max corner points
    """
    loBound = numpy.array([float('inf')] * 3)
    hiBound = numpy.array([-float('inf')] * 3)
    for shp in self.argShapes:
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
  origin2 = numpy.array([0.4, 0.5, 0.6])
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


if __name__ == '__main__': test()
