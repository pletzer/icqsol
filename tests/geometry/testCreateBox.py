#!/usr/bin/env python

"""
Test creation of box
@author pletzer@psu.edu
"""

import argparse
from icqsol.shapes.icqShapeManager import PlyShapeManager, VtkShapeManager

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

file_format = 'vtk'
file_type = 'ascii'
if args.output.find('.ply') > 0:
    file_format = 'ply'

if file_format == 'vtk':
    shape_mgr = VtkShapeManager('POLYDATA')
else:
    shape_mgr = PlyShapeManager()

s = shape_mgr.createShape('box', origin=eval(args.origin), lengths=eval(args.lengths))

if args.output:
    shape_mgr.saveShape(shape=s, file_name=args.output, file_type=file_type)
else:
    shape_mgr.show(s)
