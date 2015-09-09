#!/usr/bin/env python

"""
Test for determining whether a point is inside a shape.
@author pletzer@psu.edu
"""

import numpy
import argparse
from icqsol.shapes.icqInside import Inside
from icqsol.shapes.icqShapeManager import VtkShapeManager

parser = argparse.ArgumentParser(description='Test point inside shape')

parser.add_argument('--n_theta', type=int, dest='n_theta', default=8,  help='Set number of poloidal cells')
parser.add_argument('--n_phi', type=int, dest='n_phi', default=4, help='Set number of azimuthal cells')
args = parser.parse_args()

shape_mgr = VtkShapeManager('POLYDATA')
shp = shape_mgr.createShape('sphere', radius=1.0, origin=(0., 0., 0.), n_theta=args.n_theta, n_phi=args.n_phi)

inside = Inside(shp)

pt = numpy.array([0., 0., 0.])
# assert(inside.isInside(pt, 0.) == 1)

pt = numpy.array([1.01, 0., 0.])
# assert(inside.isInside(pt, 0.) == -1)

pt = numpy.array([0., 0.9999999999, 0.])
assert(inside.isInside(pt, 0.0) == 1)
# assert(inside.isInside(pt, 0.01) == 0)

pt = numpy.array([0., 1.0000000001, 0.])
# assert(inside.isInside(pt, 0.) == -1)
# assert(inside.isInside(pt, 0.01) == 0)
