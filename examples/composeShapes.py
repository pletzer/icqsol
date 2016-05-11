#!/usr/bin/env python

"""
Combine shape objects to create a composite object
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

parser = argparse.ArgumentParser(description='Compose shapes.')

parser.add_argument('--shapeTuples', dest='shapeTuples', nargs='+', default=[],
                    help='List of tuples <expression var>, <input file>')

parser.add_argument('--compose', dest='expression',
                    help='Expression with +, -,, and * on shapes A, B...')

parser.add_argument('--ascii', dest='ascii', action='store_true',
                    help='Save data in ASCII format (default is binary)')

parser.add_argument('--output', dest='output',
                    default='createCompositeShape-{0}.vtk'.format(tid),
                    help='Output file.')

args = parser.parse_args()

if not args.expression:
    print('ERROR: must specify --compose <expression>')
    sys.exit(2)

if len(args.shapeTuples) == 0:
    print('ERROR: must specify shape tuples: --shapeTuples <var, file>...')
    sys.exit(3)

shape_tuples = []

shape_mgr = ShapeManager()

for shapeTuple in args.shapeTuples:
    expression_var, input_file = re.sub(r'\s*', '', shapeTuple).split(',')
    if util.isVtkFile(input_file):
        vtk_dataset_type = util.getVtkDatasetType(input_file)
        shape_mgr.setReader(file_format=util.VTK_FORMAT, vtk_dataset_type=vtk_dataset_type)
    else:
        shape_mgr.setReader(file_format=util.PLY_FORMAT)
    s = shape_mgr.loadAsShape(input_file)
    shape_tuples.append((expression_var, s))

if args.ascii:
    fileType = util.ASCII
else:
    fileType = util.BINARY

compositeShape = shape_mgr.composeShapes(shape_tuples, args.expression)

if args.output:
    if util.isVtkFile(args.output):
        shape_mgr.setWriter(file_format=util.VTK_FORMAT, vtk_dataset_type=util.POLYDATA)
    else:
        shape_mgr.setWriter(file_format=util.PLY_FORMAT)
    shape_mgr.saveShape(shape=compositeShape, file_name=args.output, file_type=fileType)
