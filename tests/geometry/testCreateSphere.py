#!/usr/bin/env python

"""
Test creation of sphere
@author pletzer@psu.edu
"""

import argparse
from icqsol.shapes.icqShapeManager import ShapeManager

parser = argparse.ArgumentParser(description='Create sphere')

parser.add_argument('--output', type=str, dest='output', default='',
                    help='Set output file, suffix should be .vtk or .ply')
parser.add_argument('--radius', type=float, dest='radius', default=1.0,
                    help='Set radius')
parser.add_argument('--origin', type=str, dest='origin', default="0., 0., 0.",
                    help='Set origin, comma separated list of three floats')
parser.add_argument('--n_theta', type=int, dest='n_theta', default=8,
                    help='Set number of poloidal cells')
parser.add_argument('--n_phi', type=int, dest='n_phi', default=4,
                    help='Set number of azimuthal cells')
args = parser.parse_args()

file_format = 'vtk'
file_type = 'ascii'
if args.output.find('.ply') > 0:
    file_format = 'ply'

shape_mgr = ShapeManager(file_format)
s = shape_mgr.createShape('sphere',
                          radius=args.radius,
                          origin=eval(args.origin),
                          n_theta=args.n_theta,
                          n_phi=args.n_phi)

if args.output:
    shape_mgr.saveShape(shape=s, args.output, file_type)
else:
    shape_mgr.show(s)
