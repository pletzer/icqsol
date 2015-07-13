#!/usr/bin/env python

import vtk
import numpy

from icqShape import Shape
from csg.core import CSG
from csg.geom import Vector

class Cylinder(Shape):

  def __init__(self, radius, origin, lengths,
               n_theta=32):
    """
    Constructor
    @param radius radius
    @param origin center of low end disk
    @param lengths lengths of the cylinder along each axis
    @param n_theta number of theta cells
    """

    ori = Vector(origin[0], origin[1], origin[2])
    end = Vector(origin[0] + lengths[0],
                 origin[1] + lengths[1],
                 origin[2] + lengths[2])
    shp = CSG.cylinder(start = ori, 
                       end = end,
                       radius = radius,
                       slices = n_theta)

    Shape.__init__(self, csg=shp)

################################################################################
def test():

  cyl = Cylinder(radius=1.0, origin=(0., 0., 0.), 
    lengths=(1., 0., 0.), n_theta=8)
  cyl.save('cyl.vtk', file_format='vtk', file_type='ascii')
  cyl.show()

if __name__ == '__main__':
  test()

