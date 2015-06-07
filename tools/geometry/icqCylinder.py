#!/usr/bin/env python

import vtk
import numpy

from icqShape import Shape

class Cylinder(Shape):

  def __init__(self, radius, origin, length):
    """
    Constructor
    @param radius radius
    @param origin center of the cylinder in the x, y directions
    @param length length of the cylinder in the z direction
    """

    Shape.__init__(self)

    # infinite cylinder in the y direction
    self.cyl = vtk.vtkCylinder()
    self.cyl.SetRadius(radius)

    # rotate the cylinder so the axis is along z (initially along y)
    self.cylTransform = vtk.vtkGeneralTransform()
    self.cylTransform.PostMultiply()
    self.cylTransform.Translate(-numpy.array(origin))
    self.cylTransform.RotateX(90.0)
    self.cyl.SetTransform(self.cylTransform)


    # cut the cylinder
    self.planeLo = vtk.vtkPlane()
    self.planeLo.SetOrigin(origin[0], origin[1], origin[2] - 0.5*length)
    self.planeLo.SetNormal(0., 0., -1.)
    self.planeHi = vtk.vtkPlane()
    self.planeHi.SetOrigin(origin[0], origin[1], origin[2] + 0.5*length)
    self.planeHi.SetNormal(0., 0., 1.)

    # combine
    self.func = vtk.vtkImplicitBoolean()
    self.func.SetOperationTypeToIntersection()
    self.func.AddFunction(self.cyl)
    self.func.AddFunction(self.planeLo)
    self.func.AddFunction(self.planeHi)

    self.loBound = numpy.array([origin[0] - radius, 
                                origin[1] - radius, origin[2] - 0.5*length])
    self.hiBound = numpy.array([origin[0] + radius, 
                                origin[1] + radius, origin[2] + 0.5*length])

################################################################################
def test():

  cyl = Cylinder(radius=0.3, origin=(0., 0., 0.5), length=1.2)
  print cyl.getBounds()
  cyl.computeBoundarySurface(100, 100, 100)
  cyl.show()

if __name__ == '__main__':
  test()

