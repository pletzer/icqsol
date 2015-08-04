#!/usr/bin/env  python

from icqShape import Shape
from csg.core import CSG


class Box(Shape):

    def __init__(self, origin, lengths):
        """
        Constructor
        @param  radius  radius
        @param  origin/low  end  of  the  box
        @param  lengths  lengths  in  x,  y,  and  z
        """
        center = [origin[i] + 0.5*lengths[i] for i in range(len(origin))]
        radius = [0.5*le for le in lengths]
        shp = CSG.cube(center=center, radius=radius)

        Shape.__init__(self, csg=shp)

###############################################################################


def test():

    box = Box(origin=(0.,  0.,  0.), lengths=(0.5,  1.,  2.),)
    box.save('box.vtk', file_format='vtk', file_type='ascii')
    box.show()

if __name__ == '__main__':
    test()
