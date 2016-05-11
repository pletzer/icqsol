#!/usr/bin/env python

"""
Scale shape
"""

from __future__ import print_function
import argparse
import time
import sys
import re
import operator

from icqsol.shapes.icqShapeManager import ShapeManager
from icqsol import util

# time stamp
tid = re.sub(r'\.', '', str(time.time()))

parser = argparse.ArgumentParser(description='scale shape.')

parser.add_argument('--input', dest='input', default='',
                    help='List of input files (PLY or VTK)')

parser.add_argument('--scale', dest='scale', default='1, 1, 1',
                    help='Specify the scaling vector as three floats')

parser.add_argument('--ascii', dest='ascii', action='store_true',
                    help='Save data in ASCII format (default is binary)')

parser.add_argument('--output', dest='output',
                    default='createCompositeShape-{0}.vtk'.format(tid),
                    help='Output file.')

args = parser.parse_args()

if not args.scale:
    print('ERROR: must specify --scale <float, float, float>', end="\n")
    sys.exit(2)

if not args.input:
    print('ERROR: must specify one input file with --input <file>', end="\n")
    sys.exit(3)

# Get the format of the input - either vtk or ply.
file_format = util.getFileFormat(args.input)

if args.ascii:
    file_type = util.ASCII
else:
    file_type = util.BINARY

if file_format == util.VTK_FORMAT:
    # We have a VTK file, so Get the dataset type.
    vtk_dataset_type = util.getVtkDatasetType(args.input)
    shape_mgr = ShapeManager(file_format=file_format, vtk_dataset_type=vtk_dataset_type)
else:
    shape_mgr = ShapeManager(file_format=file_format)

pdata = shape_mgr.loadAsVtkPolyData(args.input)
scaleVec = eval(args.scale)
if reduce(operator.and_, [fact > 0.0 for fact in scaleVec], True):
     # All scaling factors must be > 0
     shape_mgr.scaleVtkPolyData(pdata, scaleVec)

if args.output:
    shape_mgr.saveVtkPolyData(pdata, file_name=args.output, file_type=file_type)
