#!/usr/bin/env python

import vtk
import numpy

class Sphere(vtk.vtkSphere):

  def __init__(self, radius, origin):
    """
    Constructor
    @param radius radius
    @param origin origin of the sphere
    """
    self.SetRadius(radius)
    self.SetCenter(origin)







