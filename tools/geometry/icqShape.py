#!/usr/bin/env python

"""
@brief A base class for constructing shapes
@author pletzer@psu.edu
"""

# standard python modules

# extensions
import vtk
import numpy
from csg.geom import Vector, Vertex, Polygon
from csg.core import CSG

class Shape:

  """
  Base class for shapes
  """

  def __init__(self, csg=None):
    """
    Constructor
    @param csg instance of csg.core.CSG
    """
    self.csg = csg

  def toPolygons(self):
    return  self.csg.toPolygons()

  def fromPolygons(self, polys):
    return Shape(csg=CSG.fromPolygons(polys))

  def rotate(self, axis=(1., 0., 0.), angleDeg=0.0):
    """
    Rotate along axis
    @param axis rotation axis
    @param angleDeg angle in degrees
    @return new Shape
    """
    res = Shape()
    res.csg = self.csg.clone()
    res.csg.rotate(axis, angleDeg)
    return res

  def translate(self, disp=(0., 0., 0.)):
    """
    Translate
    @param disp displacement
    @return new Shape
    """
    res = Shape()
    res.csg = self.csg.clone()
    res.csg.translate(disp)
    return res

  def __add__(self, other):
    """
    Union
    @param other Shape instance
    @return composite shape
    """
    return Shape(self.csg + other.csg)

  def __sub__(self, other):
    """
    Removal 
    @param other Shape instance
    @return composite shape
    """
    return Shape(self.csg - other.csg)

  def __mul__(self, other):
    """
    Intersection
    @param other Shape instance
    @return composite shape
    """
    return Shape(self.csg * other.csg)

  def toVTKPolyData(self):
    """
    Convert the data to a VTK polydata object
    """

    verts, polys, count = self.csg.toVerticesAndPolygons()

    self.points = vtk.vtkPoints()
    numPoints = len(verts)
    self.points.SetNumberOfPoints(numPoints)
    for i in range(numPoints):
      self.points.SetPoint(i, verts[i])

    pdata = vtk.vtkPolyData()
    pdata.SetPoints(self.points)
    numCells = len(polys)
    pdata.Allocate(numCells, 1)
    ptIds = vtk.vtkIdList()
    for i in range(numCells):
      npts = len(polys[i])
      ptIds.SetNumberOfIds(npts)
      for j in range(npts):
        ptIds.SetId(j, polys[i][j])
      pdata.InsertNextCell(vtk.VTK_POLYGON, ptIds)

    # may be we should also return points?
    return pdata

  @classmethod
  def load(self, file_name):
    """
    Load geometry from file
    @param file_name file name, suffix should be .ply or .vtk
    @return a Shape object
    """

    reader = None
    # select the reader according to the suffix
    if file_name.lower().find('.ply') >= 0:
      reader = vtk.vtkPLYReader()
    else:
      reader = vtk.vtkPolyDataReader()
    reader.SetFileName(file_name)
    # read
    reader.Update()

    # vtkPolyData
    pdata = reader.GetOutput()

    # store the cell connectivity as CSG polygons
    numCells = pdata.GetNumberOfPolys()
    cells = pdata.GetPolys()
    cells.InitTraversal()
    ptIds = vtk.vtkIdList()
    polygons = []
    for i in range(numCells):
      cells.GetNextCell(ptIds)
      npts = ptIds.GetNumberOfIds()
      verts = []
      for j in range(npts):
        pointIndex = ptIds.GetId(j)
        pt = pdata.GetPoint(pointIndex)
        v = Vertex(Vector(pt[0], pt[1], pt[2]))
        verts.append(v)
      polygons.append(verts)

    # instantiate the shape
    csg = CSG.fromPolygons(polygons)
    return Shape(csg=csg)
        
  def save(self, file_name, file_format, file_type):
    """
    Save the shape in file
    @param file_name file name
    @param file_format file format, currently either VTK or PLY
    @param file_type either 'ascii' or 'binary'
    """
    writer = None
    if file_format.lower() == 'ply':
      writer = vtk.vtkPLYWriter()
    else:
      writer = vtk.vtkPolyDataWriter()
    
    writer.SetFileName(file_name)
    if file_type.lower() == 'ascii':
      writer.SetFileTypeToASCII()
    else:
      writer.SetFileTypeToBinary()
    
    pdata = self.toVTKPolyData()
    if vtk.VTK_MAJOR_VERSION >= 6:
      writer.SetInputData(pdata)
    else:
      writer.SetInput(pdata)
    
    writer.Write()
    writer.Update()

  def debug(self):
    """
    Debug output of this object
    """

    points, polys, count = self.csg.toVerticesAndPolygons()

    numPoints = len(points)
    numCells = len(polys)

    print 'Number of points: ', numPoints
    for i in range(numPoints):
      p = points[i]
      print '{} {:>20} {:>20} {:>20}'.format(i, p[0], p[1], p[2])

    print 'Number of cells: ', numCells
    for i in range(numCells):
      c = polys[i]
      np = len(c)
      print 'cell index: {:>4} number of points: {:>3} point indices: '.format(i, np),
      for j in range(np):
        print '{:>8} '.format(c[j]), 
      print

  def show(self, windowSizeX=600, windowSizeY=400, filename=''):
    """
    Show the boundary surface or write image to file
    @param windowSizeX number of pixels in x
    @param windowSizeY number of pixels in y
    @param filename write to a file if this keyword is present and a 
                    non-empty string
    """

    pdata = self.toVTKPolyData()

    # create a rendering window and renderer
    try:
      ren = vtk.vtkRenderer()
      renWin = vtk.vtkRenderWindow()
      iren = vtk.vtkRenderWindowInteractor()

      camera = vtk.vtkCamera()
      mapper = vtk.vtkPolyDataMapper()
      actor = vtk.vtkActor()

      axes = [vtk.vtkArrowSource(), vtk.vtkArrowSource(), vtk.vtkArrowSource()]
      axesTransf = [vtk.vtkTransform(), vtk.vtkTransform(), vtk.vtkTransform()]
      axesTPD = [vtk.vtkTransformPolyDataFilter(),
                 vtk.vtkTransformPolyDataFilter(),
                 vtk.vtkTransformPolyDataFilter()]
      axesMappers = [vtk.vtkPolyDataMapper(), vtk.vtkPolyDataMapper(), vtk.vtkPolyDataMapper()]
      axesActors = [vtk.vtkActor(), vtk.vtkActor(), vtk.vtkActor()]

      renderLarge = vtk.vtkRenderLargeImage()
    except:
      print 'WARNING: Cannot call show method -- likely missing VTK components'
      return

    renWin.AddRenderer(ren)
    renWin.SetSize(windowSizeX, windowSizeY)
 
    # create a renderwindowinteractor
    iren.SetRenderWindow(renWin)

    # camera
    xmin, xmax, ymin, ymax, zmin, zmax = pdata.GetBounds()
    lo = numpy.array([xmin, ymin, zmin])
    hi = numpy.array([xmax, ymax, zmax])
    camera.SetFocalPoint(hi)
    center = 0.5*(lo + hi)
    camera.SetPosition(center + hi - lo)
    camera.Zoom(1.0)
    #camera.ParallelProjectionOn()
    ren.SetActiveCamera(camera)
 
    # mapper
    if vtk.VTK_MAJOR_VERSION >= 6:
      mapper.SetInputData(pdata)
    else:
      mapper.SetInput(pdata)
    mapper.ScalarVisibilityOff()

    # actor
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(1,1,1)
    
    #
    # add axes
    #

    axesColrs = [(1., 0., 0.,), (0., 1., 0.,), (0., 0., 1.,)]

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
      factor = hi[i] - lo[i]
      scale = [1., 1., 1.]
      scale[i] = factor
      axesTransf[i].Scale(scale)

    # translate to loBounds
    for at in axesTransf:
      at.Translate(lo)

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
  b = Box(origin=(0.1, 0.2, 0.3), lengths=(1.1, 1.2, 1.3))
  c = Cylinder(radius=0.5, origin=(0.3, 0.4, 0.5), lengths=(1.0, 0.0, 0.0))

  geom = c*b - s - s1

  geom.save('geom.ply', file_format='ply', file_type='binary')

  geom2 = Shape.load('geom.ply')

  geom2.show()

if __name__ == '__main__': 
  testConstructiveGeometry()







