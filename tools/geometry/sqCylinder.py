#!/usr/bin/env python

import vtk

class Cylinder(vtk.vtkCylinder):

  EPS = 1.2345678e-10

  def __init__(self, radius, origin, length):
    """
    Constructor
    @param radius radius
    @param origin center of the cylinder
    @param length length of the cylinder
    """

    self.SetRadius(radius)
    self.SetCenter(origin)






