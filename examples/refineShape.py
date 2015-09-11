#!/usr/bin/env python

"""
Refine a shape by subdividing the polygons
"""

import argparse
import time
import sys
import re

from icqsol.shapes.icqShapeManager import ShapeManager

# time stamp
tid = re.sub(r'\.', '', str(time.time()))

parser = argparse.ArgumentParser(description='Refine a shape.')

parser.add_argument('--input', dest='input', default='',
                    help='Input file (PLY or VTK)')

parser.add_argument('--refine', dest='refine', type=int, default=1,
                    help='Number of subdivisions')

parser.add_argument('--ascii', dest='ascii', action='store_true',
                    help='Save data in ASCII format (default is binary)')

parser.add_argument('--output', dest='output',
                    default='createCompositeShape-{0}.vtk'.format(tid),
                    help='Output file.')

args = parser.parse_args()

if args.refine <= 0:
    print 'ERROR: refine must be a positive, integer number: --refine #'
    sys.exit(2)

if not args.input:
    print 'ERROR: must specify input file: --input <file>'
    sys.exit(3)

# Build the shape manager.
shape_mgr = ShapeManager()
shape_mgr.setReader(file_format='vtk', vtk_dataset_type='POLYDATA')
if args.input.lower().find('.ply') >= 0:
    shape_mgr.setReader(file_format='ply')

# Read the file.
s = shape_mgr.loadAsShape(args.input)

# Refine.
for i in range(args.refine):
    s = s.refine()

# Save.
file_type = 'binary'
if args.ascii:
    file_type = 'ascii'
shape_mgr.setWriter(file_format='vtk', vtk_dataset_type='POLYDATA')
shape_mgr.saveShape(shape=s, file_name=args.output, file_type=file_type)

