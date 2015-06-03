#!/usr/bin/env python

import vtk
import numpy

class Geometry:

  def __init__(self, bxLo, bxHi):
    """
    Constructor
    @param bxLo low corner of the box
    @param bxHi high corner of the box
    """
    self.bxLo = bxLo
    self.bxHi = bxHi
    self.object = vtk.vtkImplicitBoolean()

    self.sampleFunc = vtk.vtkSampleFunction()
    self.sampleFunc.SetImplicitFunction(self.object)
    self.sampleFunc.SetModelBounds(self.bxLo[0], self.bxHi[0], \
                                   self.bxLo[1], self.bxHi[1], \
                                   self.bxLo[2], self.bxHi[2])
    self.sampleFunc.ComputeNormalsOff()

    self.surface = vtk.vtkContourFilter()
    self.surface.SetValue(0, 0.0)
    self.surface.SetInput(self.sampleFunc.GetOutput())

    self.surfPolyData = None

  def __iadd__(self, otherObj):
    """
    += operator or union
    @param otherObj instances of type Shape
    """
    self.object.SetOperationTypeToUnion()
    self.object.AddFunction(otherObj)
    return self

  def __isub__(self, otherObj):
    """
    -= operator or remove
    @param otherObj instances of type Shape
    """
    self.object.SetOperationTypeToDifference()
    self.object.AddFunction(otherObj)
    return self

  def __imul__(self, otherObj):
    """
    *= operator or intersect
    @param otherObj instances of type Shape
    """
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
    mapper.SetInput(self.surface.GetOutput())
    mapper.ScalarVisibilityOff()
 
    # actor
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(1,1,1)

    # assign actor to the renderer
    ren.AddActor(actor)
 
    # enable user interface interactor
    iren.Initialize()
    renWin.Render()
    iren.Start()

###############################################################################

def test():

  from sqSphere import Sphere

  geom = Geometry(bxLo = (0., 0., 0.), bxHi = (1., 1., 1.))
  geom += Sphere(radius=0.6, origin=(0.1, 0.2, 0.3))
  geom += Sphere(radius=0.5, origin=(0.8, 0.7, 0.6))
  geom.computeBoundarySurface(10, 10, 10)
  print geom.getBoundarySurface()
  geom.show()

if __name__ == '__main__': test()







