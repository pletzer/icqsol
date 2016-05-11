#!/usr/bin/env python

"""
@brief A base class for constructing shapes
@author alexander@gokliya.net
"""

from __future__ import print_function
from csg.core import CSG
from csg.geom import Vector
import numpy

DEFAULTS = dict(origin=[0.0, 0.0, 0.0],
                lengths=[1.0, 1.0, 1.0],
                radius=1.0,
                angle=90.0,
                n_theta=16,
                n_phi=8)

def Box(origin, lengths):
    """
    Create box
    @param  origin/low  end  of  the  box
    @param  lengths  lengths  in  x,  y,  and  z
    """
    center = [origin[i] + 0.5*lengths[i] for i in range(len(origin))]
    radius = [0.5*le for le in lengths]
    return CSG.cube(center=center, radius=radius)


def Cone(radius, origin, lengths, n_theta=16):
    """
    Create cone
    @param radius radius
    @param origin location of the focal point
    @param lengths lengths of the cone
    @param n_theta number of theta cells
    """
    ori = Vector(origin[0], origin[1], origin[2])
    end = Vector(origin[0] + lengths[0],
                 origin[1] + lengths[1],
                 origin[2] + lengths[2])
    return CSG.cone(start=ori,
                    end=end,
                    radius=radius,
                    slices=n_theta)


def Cylinder(radius, origin, lengths, n_theta=16):
    """
    Create cylinder
    @param radius radius
    @param origin center of low end disk
    @param lengths lengths of the cylinder along each axis
    @param n_theta number of theta cells
    """
    ori = Vector(origin[0], origin[1], origin[2])
    end = Vector(origin[0] + lengths[0],
                 origin[1] + lengths[1],
                 origin[2] + lengths[2])
    return CSG.cylinder(start=ori,
                        end=end,
                        radius=radius,
                        slices=n_theta)


def Sphere(radius, origin, n_theta=16, n_phi=8):
    """
    Create sphere
    @param radius radius
    @param origin center of the sphere
    @param n_theta number of theta cells
    @param n_phi number of azimuthal cells
    """
    return CSG.sphere(center=origin,
                      radius=radius,
                      slices=n_theta,
                      stacks=n_phi)


def CompositeShape(shape_tuples=[], expression=''):
    """
    @param shape_tuples list of (variable_name, shape) pairs
    @param expression expression involving +, -, and * operations.
    """
    for i in range(len(shape_tuples)):
        varName = shape_tuples[i][0]
        cmd = '{0} = shape_tuples[{1}][1]'.format(varName, i)
        exec(cmd)
    return eval(expression)
