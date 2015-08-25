#!/usr/bin/env python

"""
Combine shape objects to create a composite object
"""

import argparse
import time
import sys
import re

from icqsol.shapes.icqShapeManager import ShapeManager

# time stamp
tid = re.sub(r'\.', '', str(time.time()))

parser = argparse.ArgumentParser(description='Compose shapes.')

parser.add_argument('--input', dest='input', nargs='+', default=[],
                    help='List of input files (PLY or VTK)')

parser.add_argument('--compose', dest='expression',
                    help='Expression with +, -,, and * on shapes $0, $1...')

parser.add_argument('--ascii', dest='ascii', action='store_true',
                    help='Save data in ASCII format (default is binary)')

parser.add_argument('--output', dest='output',
                    default='createCompositeShape-{0}.vtk'.format(tid),
                    help='Output file.')

args = parser.parse_args()

if not args.expression:
    print 'ERROR: must specify --compose <expression>'
    sys.exit(2)

if len(args.input) == 0:
    print 'ERROR: must specify input files: --input <file1> <file2> ...'
    sys.exit(3)

argShapes = []
shape_mgr = ShapeManager()
for inputFile in args.input:
    s = shape_mgr.load(inputFile)
    argShapes.append(s)

expr = args.expression
for i in range(len(argShapes)):
    expr = re.sub(r'\${}'.format(i), 'argShapes[{}]'.format(i), expr)
compositeShape = eval(expr)

if args.output:
    fileFormat = 'vtk'
    fileType = 'binary'
    if args.ascii:
        fileType = 'ascii'
    if args.output.lower().find('.ply') >= 0:
        fileFormat = 'ply'
    shape_mgr.save(compositeShape,
                   args.output, file_format=fileFormat, file_type=fileType)
