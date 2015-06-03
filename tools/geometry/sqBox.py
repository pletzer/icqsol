#!/usr/bin/env python

import vtk
import numpy

class Box(vtk.vtkBox):

  def __init__(self, bxLo, bxHi):
    """
    Constructor
    @param bxLo low end of the box
    @param bxHi high end of the box
    """
    self.SetBounds(bxLo[0], bxHi[0], \
                   bxLo[1], bxHi[1], \
                   bxLo[2], bxHi[2])






