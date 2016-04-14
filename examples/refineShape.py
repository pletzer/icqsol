#!/usr/bin/env python

"""
Refine a shape by subdividing the polygons
"""

import argparse
import time
import sys
import re

from icqsol.shapes.icqShapeManager import ShapeManager
from icqsol import util

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

file_format = util.getFileFormat(args.input)

# Build the shape manager.
shape_mgr = ShapeManager()

if file_format == util.PLY_FORMAT:
    shape_mgr.setReader(file_format=util.PLY_FORMAT)
else:
    # We have a VTK file, so Get the dataset type.
    vtk_dataset_type = util.getVtkDatasetType(args.input)
    shape_mgr.setReader(file_format=util.VTK_FORMAT, vtk_dataset_type=vtk_dataset_type)

# Read the file.
s = shape_mgr.loadAsShape(args.input)

# Refine.
for i in range(args.refine):
    s = s.refine()

# Save.
if args.ascii:
    file_type = util.ASCII
else:
    file_type = util.BINARY

# Always produce VTK POLYDATA.
shape_mgr.setWriter(file_format=util.VTK_FORMAT, vtk_dataset_type=util.POLYDATA)
shape_mgr.saveShape(shape=s, file_name=args.output, file_type=file_type)

