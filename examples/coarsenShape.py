#!/usr/bin/env python

"""
Coarsen a shape by merging polygons
"""

from __future__ import print_function
import argparse
import time
import sys
import re

from icqsol.shapes.icqShapeManager import ShapeManager
from icqsol import util

# time stamp
tid = re.sub(r'\.', '', str(time.time()))

parser = argparse.ArgumentParser(description='Coarsen a shape.')

parser.add_argument('--input', dest='input', default='',
                    help='Input file (PLY or VTK)')

parser.add_argument('--min_cell_area', dest='min_cell_area', type=float, default=1.0,
                    help='Minimum cell area')

parser.add_argument('--ascii', dest='ascii', action='store_true',
                    help='Save data in ASCII format (default is binary)')

parser.add_argument('--output', dest='output',
                    default='createCompositeShape-{0}.vtk'.format(tid),
                    help='Output file.')

args = parser.parse_args()

if args.min_cell_area <= 0:
    print('ERROR: min_cell_area must be a positive, integer number: --min_cell_area #')
    sys.exit(2)

if not args.input:
    print('ERROR: must specify input file: --input <file>')
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
s = shape_mgr.loadAsVtkPolyData(args.input)

# Coarsen.
s = shape_mgr.coarsenVtkPolyData(s, min_cell_area=args.min_cell_area)

# Save.
if args.ascii:
    file_type = util.ASCII
else:
    file_type = util.BINARY

# Always produce VTK POLYDATA.
shape_mgr.setWriter(file_format=util.VTK_FORMAT, vtk_dataset_type=util.POLYDATA)

# Save the output.
shape_mgr.saveVtkPolyData(vtk_poly_data=s, file_name=args.output, file_type=file_type)

