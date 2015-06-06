#!/usr/bin/env python

import vtk
import numpy

class Cone(vtk.vtkImplicitBoolean):

  def __init__(self, angle, origin, length):
    """
    Constructor
    @param angle angle in degrees
    @param origin location of the cone vertex
    @param length length of the cone axis
    """
    # infinite cone in the x direction
    self.cone = vtk.vtkCone()
    self.cone.SetAngle(angle)

    # rotate the cylinder so the axis is along z
    self.coneTransform = vtk.vtkGeneralTransform()
    self.coneTransform.RotateY(-90.0)

    # move the code to the origin location
    self.coneTransform.Translate(origin)

    # apply the transform
    self.cone.SetTransform(self.coneTransform)

    # cut the cone
    self.planeLo = vtk.vtkPlane()
    self.planeLo.SetOrigin(origin[0], origin[1], origin[2])
    self.planeLo.SetNormal(0., 0., -1.)
    self.planeHi = vtk.vtkPlane()
    self.planeHi.SetOrigin(origin[0], origin[1], origin[2] + length)
    self.planeHi.SetNormal(0., 0., 1.)

    # combine
    self.SetOperationTypeToIntersection()
    self.AddFunction(self.cone)
    self.AddFunction(self.planeLo)
    self.AddFunction(self.planeHi)

    maxRadius = length * numpy.tan(numpy.pi*angle/180.)
    self.loBounds = numpy.array([origin[0] - maxRadius, 
                                 origin[1] - maxRadius, 
                                 origin[2]])
    self.hiBounds = numpy.array([origin[0] + maxRadius, 
                                 origin[1] + maxRadius, \
                                 origin[2] + length])

  def getBounds(self): 
    """
    Get min/max bounds
    @return low bound, hi bound
    """
    return self.loBounds, self.hiBounds

####################################################################################################

def test():

  from icqGeometry import Geometry

  cone = Cone(angle=70.0, origin=(0., 0., 0.), length = 1.0)
  print cone.getBounds()

  geom = Geometry()
  geom += cone
  geom.computeBoundarySurface(100, 100, 100)
  geom.show()

if __name__ == '__main__':
  test()
