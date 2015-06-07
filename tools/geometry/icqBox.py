#!/usr/bin/env python

import vtk
import numpy

from icqShape import Shape

class Box(Shape):

  def __init__(self, loBound, hiBound):
    """
    Constructor
    @param loBound low end of the box
    @param hiBound high end of the box
    """

    # call the base class constructor
    Shape.__init__(self)

    # that's why it's a box
    self.func = vtk.vtkBox()
    self.func.SetBounds(loBound[0], hiBound[0], \
                        loBound[1], hiBound[1], \
                        loBound[2], hiBound[2])

    self.loBound = numpy.array(loBound)
    self.hiBound = numpy.array(hiBound)

################################################################################
def test():

  box = Box(loBound=(0.1, 0.2, 0.3), hiBound=(0.9, 0.8, 0.7))
  print box.getBounds()
  box.computeBoundarySurface(100, 100, 100)
  box.show()

if __name__ == '__main__':
  test()



