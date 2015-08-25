#!/usr/bin/env python

"""
Rotate shape
"""

import argparse
import time
import sys
import re

from icqsol.shapes.icqShapeManager import ShapeManager

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

shape_mgr = ShapeManager()
shp = shape_mgr.load(args.input)
axis = eval(args.axis)
shape_mgr.rotateShape(shp, angleDeg=args.angle, axis=axis)
if args.output:
    file_format = 'vtk'
    file_type = 'binary'
    if args.ascii:
        file_type = 'ascii'
    if args.output.lower().find('.ply') >= 0:
        file_format = 'ply'
    shape_mgr.save(shp, args.output, file_format=file_format, file_type=file_type)
