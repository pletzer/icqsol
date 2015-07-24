#!/usr/bin/env python

from icqShape import Shape
from csg.core import CSG


class Sphere(Shape):

    def __init__(self, radius, origin,
                 n_theta=16, n_phi=8):
        """
        Constructor
        @param radius radius
        @param origin center of the sphere
        @param n_theta number of theta cells
        @param n_phi number of azimuthal cells
        """

        shp = CSG.sphere(center=origin,
                         radius=radius,
                         slices=n_theta,
                         stacks=n_phi)

        Shape.__init__(self, csg=shp)

##############################################################################


def test():

    sph = Sphere(radius=1.0, origin=(0., 0., 0.),
                 n_theta=8, n_phi=4)
    sph.save('sph.vtk', file_format='vtk', file_type='ascii')
    sph.show()

if __name__ == '__main__':
    test()
