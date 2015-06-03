#!/usr/bin/env python

import vtk
import numpy

class Cylinder(vtk.vtkCylinder):

  def __init__(self, radius, origin):
    """
    Constructor
    @param radius radius
    @param origin center of the cylinder
    """
    self.SetRadius(radius)
    self.SetCenter(origin)







