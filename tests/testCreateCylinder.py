#!/usr/bin/env python

"""
Test creation of cylinders
@author alexander.net
"""

import argparse
from icqsol.shapes.icqShapeManager import ShapeManager
from icqsol import util

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

# Get the format of the input - either vtk or ply.
file_format = util.getFileFormat(args.output)

if file_format == util.VTK_FORMAT:
    shape_mgr = ShapeManager(file_format=file_format, vtk_dataset_type=util.POLYDATA)
else:
    shape_mgr = ShapeManager(file_format=file_format)

s = shape_mgr.createShape('cylinder',
                          origin=eval(args.origin),
                          lengths=eval(args.lengths),
                          radius=args.radius,
                          n_theta=args.n_theta)

if args.output:
    shape_mgr.saveShape(shape=s, file_name=args.output, file_type=util.ASCII)
else:
    shape_mgr.showShape(s)
