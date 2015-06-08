#!/usr/bin/env python

import vtk
import numpy

from icqShape import Shape

class Sphere(Shape):

  def __init__(self, radius, origin):
    """
    Constructor
    @param radius radius
    @param origin origin of the sphere
    """

    # call base class constructor
    Shape.__init__(self)

    self.func = vtk.vtkSphere()
    self.func.SetRadius(radius)
    ori = numpy.array(origin)
    self.func.SetCenter(ori)
    self.updateBounds(ori - radius, ori + radius)

################################################################################
def test():

  sphere = Sphere(radius=1.0, origin=(0., 0., 0.2))
  print sphere.getBounds()
  sphere.computeBoundarySurface(100, 100, 100)
  sphere.show()

if __name__ == '__main__':
  test()

