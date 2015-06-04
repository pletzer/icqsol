#!/usr/bin/env python

import vtk
import numpy

class Cylinder(vtk.vtkImplicitBoolean):

  def __init__(self, radius, origin, length):
    """
    Constructor
    @param radius radius
    @param origin center of the cylinder
    @param length length of the cylinder
    """
    self.cyl = vtk.vtkCylinder()
    self.cyl.SetRadius(radius)
    self.cyl.SetCenter(origin)

    self.cylTransform = vtk.vtkGeneralTransform()
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
    self.SetOperationTypeToIntersection()
    self.AddFunction(self.cyl)
    self.AddFunction(self.planeLo)
    self.AddFunction(self.planeHi)

