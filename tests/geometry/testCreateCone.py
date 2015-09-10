#!/usr/bin/env python

"""
Test creation of cone
@author pletzer@psu.edu
"""

import argparse
from icqsol.shapes.icqShapeManager import ShapeManager

parser = argparse.ArgumentParser(description='Create cone')

parser.add_argument('--output', type=str, dest='output', default='',
                    help='Set output file, suffix should be .vtk or .ply')
parser.add_argument('--radius', type=float, dest='radius', default=1.0,
                    help='Set radius')
parser.add_argument('--lengths', type=str, dest='lengths',
                    default="1.0, 0.0, 0.0",
                    help='Set lengths')
parser.add_argument('--origin', type=str, dest='origin', default="0.,0.,0.",
                    help='Set origin, comma separated list of three floats')
parser.add_argument('--n_theta', type=int, dest='n_theta', default=8,
                    help='Set number of poloidal cells')
args = parser.parse_args()

file_format = 'vtk'
if args.output.find('.ply') > 0:
    file_format = 'ply'

file_type = 'ascii'

if file_format == 'vtk':
    shape_mgr = ShapeManager(file_format=file_format, vtk_dataset_type='POLYDATA')
else:
    shape_mgr = ShapeManager(file_format=file_format)

s = shape_mgr.createShape('cone',
                          origin=eval(args.origin),
                          lengths=eval(args.lengths),
                          radius=args.radius,
                          n_theta=args.n_theta)

if args.output:
    shape_mgr.saveShape(shape=s, file_name=args.output, file_type=file_type)
else:
    shape_mgr.show(s)
