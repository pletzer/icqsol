#!/usr/bin/env python

"""
Translate shape
"""

import argparse
import time
import sys
import re

from icqsol.shapes.icqShapeManager import ShapeManager
from icqsol import util

# time stamp
tid = re.sub(r'\.', '', str(time.time()))

parser = argparse.ArgumentParser(description='Translate shape.')

parser.add_argument('--input', dest='input', default='',
                    help='List of input files (PLY or VTK)')

parser.add_argument('--translate', dest='translation',
                    help='Specify the translation vector as three floats')

parser.add_argument('--ascii', dest='ascii', action='store_true',
                    help='Save data in ASCII format (default is binary)')

parser.add_argument('--output', dest='output',
                    default='createCompositeShape-{0}.vtk'.format(tid),
                    help='Output file.')

args = parser.parse_args()

if not args.translation:
    print 'ERROR: must specify --translate <float, float, float>'
    sys.exit(2)

if not args.input:
    print 'ERROR: must specify one input file with --input <file>'
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
transVec = eval(args.translation)
shape_mgr.translateVtkPolyData(pdata, transVec)
if args.output:
    shape_mgr.saveVtkPolyData(pdata, file_name=args.output, file_type=file_type)
