#!/usr/bin/env python

"""
Test creation of box
@author alexander@gokliya.net
"""

import argparse
from icqsol.shapes.icqShapeManager import ShapeManager
from icqsol import util

parser = argparse.ArgumentParser(description='Create box')

parser.add_argument('--output', type=str, dest='output', default='',
                    help='Set output file, suffix should be .vtk or .ply')
parser.add_argument('--lengths', type=str, dest='lengths',
                    default="1., 1., 1.",
                    help='Set box lengths')
parser.add_argument('--origin', type=str, dest='origin',
                    default="0., 0., 0.",
                    help='Set origin, 3 floats')
args = parser.parse_args()

# Get the format of the input - either vtk or ply.
file_format = util.getFileFormat(args.output)

if file_format == util.VTK_FORMAT:
    shape_mgr = ShapeManager(file_format=file_format, vtk_dataset_type=util.POLYDATA)
else:
    shape_mgr = ShapeManager(file_format=file_format)

s = shape_mgr.createShape('box', origin=eval(args.origin), lengths=eval(args.lengths))

if args.output:
    shape_mgr.saveShape(shape=s, file_name=args.output, file_type=util.ASCII)
else:
    shape_mgr.showShape(s)
