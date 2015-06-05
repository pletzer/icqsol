#!/usr/bin/env python

# standard python modules
import re

# extensions
import vtk
import numpy

# local includes
from icqSphere import Sphere
from icqCylinder import Cylinder
from icqBox import Box

class Geometry:

  def __init__(self):
    """
    Constructor
   """
    self.object = vtk.vtkImplicitBoolean()

    self.sampleFunc = vtk.vtkSampleFunction()
    self.sampleFunc.SetImplicitFunction(self.object)
    self.sampleFunc.ComputeNormalsOff()

    self.surface = vtk.vtkContourFilter()
    self.surface.SetValue(0, 0.0)
    if vtk.VTK_MAJOR_VERSION >= 6:
      self.surface.SetInputConnection(self.sampleFunc.GetOutputPort())
    else:
      self.surface.SetInput(self.sampleFunc.GetOutput())

    self.surfPolyData = None

    self.objectList = []
    self.loBound = numpy.array([float('inf')] * 3)
    self.hiBound = numpy.array([-float('inf')] * 3)
    self.boundsSet = False

  def __iadd__(self, otherObj):
    """
    += operator or union
    @param otherObj other shape
    """
    self.objectList.append(otherObj)
    self.object.SetOperationTypeToUnion()
    self.object.AddFunction(otherObj)
    return self

  def __isub__(self, otherObj):
    """
    -= operator or remove
    @param otherObj other shape
    """
    self.objectList.append(otherObj)
    self.object.SetOperationTypeToDifference()
    self.object.AddFunction(otherObj)
    return self

  def __imul__(self, otherObj):
    """
    *= operator or intersect
    @param otherObj instances of type Shape
    """
    self.objectList.append(otherObj)
    self.object.SetOperationTypeToIntersection()
    self.object.AddFunction(otherObj)
    return self

  def computeBoundarySurface(self, nx, ny, nz):
    """
    Discretize the boundary 
    @param nx number of cells in the x direction
    @param ny number of cells in the y direction
    @param nz number of cells in the z direction
    @return {'points': array of points, 'cells': array of cells}
    """

    for o in self.objectList:
      objLo, objHi = o.getBounds()
      self.loBound = numpy.minimum(self.loBound, objLo)
      self.hiBound = numpy.maximum(self.hiBound, objHi)

    self.boundsSet = True

    # set the bounds for th esampling function
    self.sampleFunc.SetModelBounds(self.loBound[0], self.hiBound[0], \
                                   self.loBound[1], self.hiBound[1], \
                                   self.loBound[2], self.hiBound[2])


    self.sampleFunc.SetSampleDimensions(nx, ny, nz)
    self.surface.Update()
    self.surfPolyData = self.surface.GetOutput()

  def getBoundarySurface(self):
    """
    Discretize the boundary 
    @return {'points': array of points, 'cells': array of cells}
    """

    # gather the boundary vertices
    points = self.surfPolyData.GetPoints()
    numPoints = points.GetNumberOfPoints()
    pointArr = numpy.zeros( (numPoints, 3), numpy.float32 )
    for i in range(numPoints):
      pointArr[i, :] = points.GetPoint(i)

    # gather the triangles
    numCells = self.surfPolyData.GetNumberOfCells()
    cellArr = numpy.zeros( (numCells, 3), numpy.int )
    for i in range(numCells):
      cell = self.surfPolyData.GetCell(i)
      ptIds = cell.GetPointIds()
      npts = ptIds.GetNumberOfIds()
      pInds = [0 for j in range(npts)]
      for j in range(npts):
        pInds[j] = int(ptIds.GetId(j))
      cellArr[i, :] = pInds

    return {'points': pointArr, 'cells': cellArr}

  def show(self, windowSizeX=600, windowSizeY=400):
    """
    Show the boundary surface
    @param windowSizeX number of pixels in x
    @param windowSizeY number of pixels in y
    """
    # create a rendering window and renderer
    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    renWin.SetSize(windowSizeX, windowSizeY)
 
    # create a renderwindowinteractor
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)
 
    # mapper
    mapper = vtk.vtkPolyDataMapper()
    if vtk.VTK_MAJOR_VERSION >= 6:
      mapper.SetInputConnection(self.surface.GetOutputPort())
    else:
      mapper.SetInput(self.surface.GetOutput())
    mapper.ScalarVisibilityOff()
 
    # actor
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(1,1,1)
    
    #
    # add axes
    #

    axes = [vtk.vtkArrowSource(), vtk.vtkArrowSource(), vtk.vtkArrowSource()]
    axesColrs = [(1., 0., 0.,), (0., 1., 0.,), (0., 0., 1.,)]
    axesTransf = [vtk.vtkTransform(), vtk.vtkTransform(), vtk.vtkTransform()]

    
    if self.boundsSet:
      # scale
      for i in range(3):
        factor = self.hiBound[i] - self.loBound[i]
        axesTransf[i].Scale(factor, factor, factor)

      # translate to loBounds
      for at in axesTransf:
        at.Translate(self.loBound)

    # rotate the y and z arrows (initially along x). Order of operations is
    # last first (rotation before translation before scaling)
    axesTransf[1].RotateZ(90.0) 
    axesTransf[2].RotateY(-90.0)

    axesTPD = [vtk.vtkTransformPolyDataFilter(), vtk.vtkTransformPolyDataFilter(), vtk.vtkTransformPolyDataFilter()]
    axesMappers = [vtk.vtkPolyDataMapper(), vtk.vtkPolyDataMapper(), vtk.vtkPolyDataMapper()]
    axesActors = [vtk.vtkActor(), vtk.vtkActor(), vtk.vtkActor()]
    for i in range(3):
      axesTPD[i].SetInputConnection(axes[i].GetOutputPort())
      axesTPD[i].SetTransform(axesTransf[i])
      axesMappers[i].SetInputConnection(axesTPD[i].GetOutputPort())
      axesActors[i].SetMapper(axesMappers[i])
      axesActors[i].GetProperty().SetColor(axesColrs[i])
      ren.AddActor(axesActors[i])

    # assign actor to the renderer
    ren.AddActor(actor)
 
    # enable user interface interactor
    iren.Initialize()
    renWin.Render()
    iren.Start()

  def applyPrefixExpression(self, prefixExpr, argNameVals):
    """
    Apply prefix expression 
    @param prefixExpression, e.g. (* (+ s1 s2) s3)
    @param argNameVals dictionary for argument name and corresponding values
    """

    # map betwen operator and union or intersection
    opMap = {
      '+': self.__iadd__,
      '*': self.__imul__,
      '-': self.__isub__,
    }

    expr = self.squeeze(prefixExpr)

    # build the geometry recursively by reducing the expression until there are no 
    # more parentheses, starting from the right most nested ()
    while len(expr) > 0:

      posBeg = expr.rfind("(")
      posEnd = expr.find(")", posBeg)

      subExpr = expr[posBeg + 1: posEnd]
      tokens = subExpr.split()

      op = opMap[tokens[0]]
      args = tokens[1:]

      # apply the operation on all the arguments
      for a in args:
        argVal = argNameVals[a]
        op.__call__(argVal)

      # remove subExpr from expr
      expr = self.squeeze(expr[:posBeg] + expr[posEnd + 1:])

  def squeeze(self, expr):
    """
    Squeeze all the spaces out
    @param expr expression
    @return new expression
    """
    expr = re.sub(r'\s+', ' ', expr)
    expr = re.sub(r'\(\s*', '(', expr)
    expr = re.sub(r'\(\s*', '(', expr)
    expr = re.sub(r'^\s*', '', expr)
    expr = re.sub(r'\s*$', '', expr)
    return expr


###############################################################################

def testConstructiveGeometry():

  geom = Geometry()
  geom += Sphere(radius=0.6, origin=(0.1, 0.2, 0.3))
  geom *= Box(bxLo=(0.1, 0.2, 0.3), bxHi=(1.0, 1.0, 1.0))
  geom -= Cylinder(radius=0.5, origin=(0.3, 0.4, 0.5), length=0.6)
  geom.computeBoundarySurface(100, 100, 100)
  print geom.getBoundarySurface()
  geom.show()

def testApplyPrefixExpression():

  s = Sphere(radius=0.6, origin=(0.1, 0.2, 0.3))
  b = Box(bxLo=(0.1, 0.2, 0.3), bxHi=(1.0, 1.0, 1.0))
  c = Cylinder(radius=0.5, origin=(0.3, 0.4, 0.5), length=0.6)

  geom = Geometry()
  geom.applyPrefixExpression('(- (* (+ s) b) c)', 
    {'s': s, 'b': b, 'c': c})
  geom.computeBoundarySurface(100, 100, 100)
  geom.show()

if __name__ == '__main__': 
  #testConstructiveGeometry()
  testApplyPrefixExpression()







