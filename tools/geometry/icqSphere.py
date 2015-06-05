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
    self.SetCenter(numpy.array(origin))
    self.loBounds = numpy.array(origin) - radius
    self.hiBounds = numpy.array(origin) + radius

  def getBounds(self): 
    """
    Get min/max bounds
    @return low bound, hi bound
    """
    return self.loBounds, self.hiBounds






