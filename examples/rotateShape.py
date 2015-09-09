#!/usr/bin/env python

"""
Rotate shape
"""

import argparse
import time
import sys
import re

from icqsol.shapes.icqShapeManager import PlyShapeManager, VtkShapeManager

# time stamp
tid = re.sub(r'\.', '', str(time.time()))

parser = argparse.ArgumentParser(description='Translate shape.')

parser.add_argument('--input', dest='input', default='',
                    help='List of input files (PLY or VTK)')

parser.add_argument('--angle', dest='angle', type=float, default=0.0,
                    help='Specify rotation angle in degrees')

parser.add_argument('--axis', dest='axis', default="0., 0., 1.",
                    help='Specify rotation axis (3 floating point numbers)')

parser.add_argument('--ascii', dest='ascii', action='store_true',
                    help='Save data in ASCII format (default is binary)')

parser.add_argument('--output',
                    dest='output',
                    default='createCompositeShape-{0}.vtk'.format(tid),
                    help='Output file.')

args = parser.parse_args()

if not args.input:
    print 'ERROR: must specify one input file with --input <file>'
    sys.exit(3)

file_format = 'vtk'
file_type = 'binary'
if args.ascii:
    file_type = 'ascii'
if args.output.lower().find('.ply') >= 0:
    file_format = 'ply'

if file_format == 'vtk':
    shape_mgr = VtkShapeManager('POLYDATA')
else:
    shape_mgr = PlyShapeManager()

shp = shape_mgr.loadAsShape(args.input)
axis = eval(args.axis)
shape_mgr.rotateShape(shp, angleDeg=args.angle, axis=axis)
if args.output:
    shape_mgr.saveShape(shape=shp, file_name=args.output, file_type=file_type)
