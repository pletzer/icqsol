#!/usr/bin/env python

"""
Combine shape objects to create a composite object
"""

import argparse
import time
import sys
import re

from icqsol.shapes.icqShapeManager import PlyShapeManager, VtkShapeManager

def isVtkFile(file_name):
    return file_name.lower().endswith('.vtk')

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
    print 'ERROR: must specify --compose <expression>'
    sys.exit(2)

if len(args.shapeTuples) == 0:
    print 'ERROR: must specify shape tuples: --shapeTuples <var, file>...'
    sys.exit(3)

shape_tuples = []

for shapeTuple in args.shapeTuples:
    expression_var, input_file = re.sub(r'\s*', '', shapeTuple).split(',')
    if isVtkFile(input_file):
        shape_mgr = VtkShapeManager('POLYDATA')
    else:
        shape_mgr = PlyShapeManager()
    s = shape_mgr.loadAsShape(input_file)
    shape_tuples.append((expression_var, s))

fileType = 'binary'
if args.ascii:
    fileType = 'ascii'

if isVtkFile(args.output):
    shape_mgr = VtkShapeManager('POLYDATA')
else:
    shape_mgr = PlyShapeManager()

compositeShape = shape_mgr.composeShapes(shape_tuples, args.expression)

if args.output:

    shape_mgr.saveShape(shape=compositeShape, file_name=args.output, file_type=fileType)
