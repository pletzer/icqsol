#!/usr/bin/env python

# standard python modules
import re

# extensions
import vtk
import numpy

class Shape:

  """
  Base class for shapes
  """

  def __init__(self):
    """
    Constructor
    """

    # implicit function, either boolean or non-boolean
    self.func = None

    # low/high corners of hte box containing the object. 
    self.loBound = numpy.array([float('inf')] * 3)
    self.hiBound = numpy.array([-float('inf')] * 3)

    # this is to compute the boundary points
    self.sampleFunc = vtk.vtkSampleFunction()
    self.sampleFunc.ComputeNormalsOff()

    self.surface = vtk.vtkContourFilter()
    self.surface.SetValue(0, 0.0)
    if vtk.VTK_MAJOR_VERSION >= 6:
      self.surface.SetInputConnection(self.sampleFunc.GetOutputPort())
    else:
      self.surface.SetInput(self.sampleFunc.GetOutput())

    self.surfPolyData = None

    # store all the transformations 
    self.transf = vtk.vtkTransform()
    # I prefer to apply the rotations after creation of the object
    self.transf.PostMultiply()

  def updateBounds(self, loBound, hiBound):
    """
    Update the domain bounds
    @param loBound low corner 
    @param hiBound high corner 
    """
    self.loBound = numpy.minimum(self.loBound, numpy.array(loBound))
    self.hiBound = numpy.maximum(self.hiBound, numpy.array(hiBound))

  def getBounds(self):
    """
    Get the domain bounds
    @return bounds
    """
    return self.loBound, self.hiBound

  def rotate(self, axis, angleDeg):
    """
    Rotate object around axis
    @param axis axis index (0 for x, 1 for y, and 2 for z)
    @param angleDeg angle in degrees (counterclockwise is positive)
    """
    # must have created object
    if self.func == None:
      # a no-operation
      return

    self.func.SetTransform(self.transf)

    if axis == 0:
      self.transf.RotateX(angleDeg)

    elif axis == 1:
      self.transf.RotateY(angleDeg)

    elif axis == 2:
      self.transf.RotateZ(angleDeg)

  def __add__(self, otherShape):
    """
    + operator or union
    @param otherShape other shape
    @return composite shape
    """
    res = Shape()

    res.updateBounds(self.loBound, self.hiBound)
    res.updateBounds(otherShape.loBound, otherShape.hiBound)

    res.func = vtk.vtkImplicitBoolean()
    res.func.SetOperationTypeToUnion()
    res.func.AddFunction(self.func)
    res.func.AddFunction(otherShape.func)

    return res

  def __mul__(self, otherShape):
    """
    * operator or intersect
    @param otherShape other shape
    @return composite shape
    """
    res = Shape()

    res.updateBounds(self.loBound, self.hiBound)
    res.updateBounds(otherShape.loBound, otherShape.hiBound)

    res.func = vtk.vtkImplicitBoolean()
    res.func.SetOperationTypeToIntersection()
    res.func.AddFunction(self.func)
    res.func.AddFunction(otherShape.func)

    return res

  def __sub__(self, otherShape):
    """
    - operator or remove
    @param otherShape other shape
    @return composite shape
    """
    res = Shape()

    res.updateBounds(self.loBound, self.hiBound)
    res.updateBounds(otherShape.loBound, otherShape.hiBound)

    res.func = vtk.vtkImplicitBoolean()
    res.func.SetOperationTypeToDifference()
    res.func.AddFunction(self.func)
    res.func.AddFunction(otherShape.func)

    return res

  def computeBoundarySurface(self, nx, ny, nz):
    """
    Discretize the boundary 
    @param nx number of cells in the x direction
    @param ny number of cells in the y direction
    @param nz number of cells in the z direction
    @return {'points': array of points, 'cells': array of cells}
    """

    self.sampleFunc.SetImplicitFunction(self.func)

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

  def show(self, windowSizeX=600, windowSizeY=400, filename=''):
    """
    Show the boundary surface or write to file
    @param windowSizeX number of pixels in x
    @param windowSizeY number of pixels in y
    @param filename write to a file if this keyword is present and a 
                    non-empty string
    """
    # create a rendering window and renderer
    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    renWin.SetSize(windowSizeX, windowSizeY)
 
    # create a renderwindowinteractor
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)

    # camera
    camera = vtk.vtkCamera()
    camera.SetFocalPoint(self.hiBound)
    center = 0.5*(self.loBound + self.hiBound)
    camera.SetPosition(center + self.hiBound - self.loBound)
    camera.Zoom(1.0)
    #camera.ParallelProjectionOn()
    ren.SetActiveCamera(camera)
 
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

    for a in axes:
      a.SetShaftRadius(0.01)
      a.SetTipLength(0.2)
      a.SetTipRadius(0.03)

    for at in axesTransf:
      at.PostMultiply()
    
    # rotate the y and z arrows (initially along x)
    axesTransf[1].RotateZ(90.0) 
    axesTransf[2].RotateY(-90.0)

    # scale
    for i in range(3):
      factor = self.hiBound[i] - self.loBound[i]
      scale = [1., 1., 1.]
      scale[i] = factor
      axesTransf[i].Scale(scale)

    # translate to loBounds
    for at in axesTransf:
      at.Translate(self.loBound)

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

    # write to file
    writer = None
    renderLarge = vtk.vtkRenderLargeImage()
    if filename:
      if filename.lower().find('.png') > 0:
        writer = vtk.vtkPNGWriter()
      elif filename.lower().find('.jp') > 0:
        writer = vtk.vtkJPEGWriter()
      elif filename.lower().find('.tiff') > 0:
        writer = vtk.vtkTIFFWriter()
      if writer:
        renderLarge.SetInput(ren)
        renderLarge.SetMagnification(1)
        renderLarge.Update()
        writer.SetFileName(filename)
        writer.SetInputConnection(renderLarge.GetOutputPort())
        writer.Write()
    else:
      # fire up interactor
      iren.Initialize()
      renWin.Render()
      iren.Start()

###############################################################################

def testConstructiveGeometry():

  from icqSphere import Sphere
  from icqBox import Box
  from icqCylinder import Cylinder

  s1 = Sphere(radius=0.7, origin=(0., 0., 0.))
  s = Sphere(radius=0.2, origin=(0.1, 0.2, 0.3))
  b = Box(loBound=(0.1, 0.2, 0.3), hiBound=(1.0, 1.0, 1.0))
  c = Cylinder(radius=0.5, origin=(0.3, 0.4, 0.5), length=0.6)

  geom = c*b - s - s1

  geom.computeBoundarySurface(100, 100, 100)
  print geom.getBoundarySurface()
  geom.show()

if __name__ == '__main__': 
  testConstructiveGeometry()







