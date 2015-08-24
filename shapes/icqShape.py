#!/usr/bin/env python

"""
@brief A base class for constructing shapes
@author pletzer@psu.edu
"""

from csg.core import CSG
from csg.geom import Vector

DEFAULTS = dict(origin=[0.0, 0.0, 0.0],
                lengths=[1.0, 1.0, 1.0],
                radius=[0.5, 0.5, 0.5],
                angle=90.0,
                n_theta=16,
                n_phi=8)


class Shape(object):
    """
    Base class for shapes
    """
    def __init__(self, csg=None):
        """
        Constructor
        @param csg instance of csg.core.CSG
        """
        self.csg = csg

    def debug(self):
        """
        Print/Debug output of this object
        """
        points, polys, count = self.csg.toVerticesAndPolygons()

        numPoints = len(points)
        numCells = len(polys)

        print 'Number of points: ', numPoints
        for i in range(numPoints):
            p = points[i]
            print '{} {:>20} {:>20} {:>20}'.format(i, p[0], p[1], p[2])

        print 'Number of cells: ', numCells
        for i in range(numCells):
            c = polys[i]
            np = len(c)
            print 'cell: {:>4} num points: {:>3} point: '.format(i, np),
            for j in range(np):
                print '{:>8} '.format(c[j]),
            print


class Box(Shape):

    def __init__(self, origin, lengths):
        """
        Constructor
        @param  origin/low  end  of  the  box
        @param  lengths  lengths  in  x,  y,  and  z
        """
        center = [origin[i] + 0.5*lengths[i] for i in range(len(origin))]
        radius = [0.5*le for le in lengths]
        self.csg = CSG.cube(center=center, radius=radius)


class Cone(Shape):

    def __init__(self, radius, origin, lengths, n_theta=16):
        """
        Constructor
        @param radius radius
        @param origin location of the focal point
        @param lengths lengths of the cone
        @param n_theta number of theta cells
        """

        ori = Vector(origin[0], origin[1], origin[2])
        end = Vector(origin[0] + lengths[0],
                     origin[1] + lengths[1],
                     origin[2] + lengths[2])
        self.csg = CSG.cone(start=ori,
                            end=end,
                            radius=radius,
                            slices=n_theta)


class Cylinder(Shape):

    def __init__(self, radius, origin, lengths, n_theta=16):
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
        self.csg = CSG.cylinder(start=ori,
                                end=end,
                                radius=radius,
                                slices=n_theta)


class Sphere(Shape):

    def __init__(self, radius, origin, n_theta=16, n_phi=8):
        """
        Constructor
        @param radius radius
        @param origin center of the sphere
        @param n_theta number of theta cells
        @param n_phi number of azimuthal cells
        """
        self.csg = CSG.sphere(center=origin,
                              radius=radius,
                              slices=n_theta,
                              stacks=n_phi)


def CompositeShape(shape_tuples=[], expression=''):
    """
    @param shape_tuples list of (variable_name, shape) pairs
    @param expression expression involving +, -, and * operations.
    """
    for i in range(len(shape_tuples)):
        exec(shape_tuples[i][0] + ' = shape_tuples[i][1]')
    return eval(expression)
