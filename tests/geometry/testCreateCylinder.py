#!/usr/bin/env python

"""
Test creation of cylinders
@author pletzer@psu.edu
"""

import argparse
from icqsol.shapes.shape_manager import ShapeManager

parser = argparse.ArgumentParser(description='Create cylinder')

parser.add_argument('--output', type=str, dest='output', default='',
                    help='Set output file, suffix should be .vtk or .ply')
parser.add_argument('--radius', type=float, dest='radius', default=1.0,
                    help='Set radius')
parser.add_argument('--lengths', type=str, dest='lengths',
                    default="1.0, 0.0, 0.0",
                    help='Set lengths')
parser.add_argument('--origin', type=str, dest='origin', default="0.,0.,0.",
                    help='Set origin, comma separated list of three floats')
parser.add_argument('--n_theta', type=int, dest='n_theta', default=32,
                    help='Set number of poloidal cells')
args = parser.parse_args()

shape_mgr = ShapeManager()
s = shape_mgr.createShape('cylinder',
                          radius=args.radius,
                          origin=eval(args.origin),
                          lengths=eval(args.lengths),
                          n_theta=args.n_theta)

if args.output:
    file_format = 'vtk'
    file_type = 'ascii'
    if args.output.find('.ply') > 0:
        file_format = 'ply'
    shape_mgr.save(s, args.output, file_format, file_type)
else:
    shape_mgr.show(s)
