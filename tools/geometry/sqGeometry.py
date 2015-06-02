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


  def __iadd__(self, otherObj):
    """
    Add object (create union)
    @param otherObj instances of type Shape
    """
    self.object.SetOperationTypeToUnion()
    self.object.AddFunction(otherObj)

  def __isub__(self, otherObj):
    """
    Remove object
    @param otherObj instances of type Shape
    """
    self.object.SetOperationTypeToDifference()
    self.object.AddFunction(otherObj)

  def __imul__(self, otherObj):
    """
    Intersect with other object
    @param otherObj instances of type Shape
    """
    self.object.SetOperationTypeToIntersection()
    self.object.AddFunction(otherObj)

  def getBoundary(self, nx, ny, nz):
    """
    Discretize the boundary 
    @param nx number of cells in the x direction
    @param ny number of cells in the y direction
    @param nz number of cells in the z direction
    @return {'points': array of points, 'cells': array of cells}
    """

    sampleFunc = vtk.vtkSampleFunction()
    sampleFunc.SetImplicitFunction(self.object)
    sampleFunc.SetModelBounds(self.bxLo[0], self.bxHi[0], \
                              self.bxLo[1], self.bxHi[1], \
                              self.bxLo[2], self.bxHi[2])
    sampleFunc.SetSampleDimensions(nx, ny, nz)
    sampleFunc.ComputeNormalsOff()

    surface = vtk.vtkContourFilter()
    surface.SetValue(0, 0.0)
    surface.SetInput( sampleFunc.GetOutput() )
    surface.Update()

    polydata = surfaceGetOutput()

    # gather the boundary vertices
    points = polydata.GetPoints()
    numPoints = points.GetNumberOfPoints()
    pointArr = numpy.zeros( (numPoints, 3), numpy.float32 )
    for i in range(numPoints):
      pointArr[i, :] = points.GetPoint(i)

    # gather the triangles
    numCells = polydata.GetNumberOfCells()
    cellArr = numpy.zeros( (numCells, 3), numpy.int )
    for i in range(numCells):
      cell = polydata.GetCell(i)
      ptIds = cell.GetPointIds()
      npts = ptIds.GetNumberOfIds()
      pInds = [0 for j in range(npts)]
      for j in range(npts):
        pInds[j] = int(ptIds.GetId(j))
      cellArr[i, :] = pInds

    return {'points': pointArr, 'cells': cellArr}
   

###############################################################################

def test():

  geom = Geometry(bxLo = (0., 0.), bxHi = (1., 1.))
  geom += Sphere(radius=0.4, origin=(0.1, 0.2, 0.3))
  geom += Sphere(radius=0.6, origin=(0.4, 0.3, 0.2))

if __name__ == '__main__': test()







