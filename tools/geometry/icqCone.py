#!/usr/bin/env python

import vtk
import numpy

from icqShape import Shape

class Cone(Shape):

  def __init__(self, angle, origin, length):
    """
    Constructor
    @param angle half angle in degrees
    @param origin location of the cone vertex
    @param length length of the cone axis
    """

    # call base class constructor
    Shape.__init__(self)

    # infinite cone in the x direction
    self.func = vtk.vtkImplicitBoolean()

    self.cone = vtk.vtkCone()
    self.cone.SetAngle(angle)

    self.coneTransform = vtk.vtkTransform()
    self.coneTransform.PostMultiply()

    # move the cone to the origin location
    self.coneTransform.Translate(-numpy.array(origin))


    # apply the transform    # make rotation axis z (initially along x)
    self.coneTransform.RotateY(-90.0)

    self.cone.SetTransform(self.coneTransform)

    # cut the cone
    self.planeLo = vtk.vtkPlane()
    self.planeLo.SetOrigin(origin[0], origin[1], origin[2])
    self.planeLo.SetNormal(0., 0., -1.)
    self.planeHi = vtk.vtkPlane()
    self.planeHi.SetOrigin(origin[0], origin[1], origin[2] + length)
    self.planeHi.SetNormal(0., 0., 1.)

    # combine
    self.func.SetOperationTypeToIntersection()
    self.func.AddFunction(self.cone)
    self.func.AddFunction(self.planeLo)
    self.func.AddFunction(self.planeHi)

    maxRadius = length * numpy.tan(numpy.pi*angle/180.)
    self.loBound = numpy.array([origin[0] - maxRadius, 
                                origin[1] - maxRadius, 
                                origin[2]])
    self.hiBound = numpy.array([origin[0] + maxRadius, 
                                origin[1] + maxRadius, \
                                origin[2] + length])

####################################################################################################

def test():

  cone = Cone(angle=70.0, origin=(0., 0., 0.2), length = 1.0)
  print cone.getBounds()
  cone.computeBoundarySurface(100, 100, 100)
  cone.show()

if __name__ == '__main__':
  test()
