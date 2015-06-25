#!/usr/bin/env python

"""
Combine multiple shape objects
"""

import argparse
import time
import sys
import re
import numpy
from numpy import cos, sin, pi

from icqsol.tools.geometry.icqCompositeShape import CompositeShape
from icqsol.tools.geometry.icqBaseShape import BaseShape

# time stamp
tid = re.sub(r'\.', '', str(time.time()))

parser = argparse.ArgumentParser(description='Create composite shape.')

parser.add_argument('--input', dest='input', nargs='+', default=[],
  help='List of input files (PLY or VTK)')

parser.add_argument('--compose', dest='expression',
	help='An expression containing + (union) - (removal) and * (intersection) operations. Shape variable names are $0, $1, ...')

parser.add_argument('--maxTriangleArea', dest='maxTriangleArea', default=0.1, type=float, 
  help='Max area of triangles')

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
  print 'ERROR: must specify at least one input file with --input <file1> <file2> ...'
  sys.exit(3)

shp = CompositeShape()
argShapes = []
for inputFile in args.input:
  s = BaseShape()
  s.load(inputFile)
  argShapes.append(s)

shp.compose(args.expression, argShapes)
shp.computeSurfaceMeshes(args.maxTriangleArea)

if args.output:
  fileFormat = 'vtk'
  fileType = 'binary'
  if args.ascii:
    fileType = 'ascii'
  if args.output.lower().find('.ply') > 0:
    fileFormat = 'ply'
  shp.save(args.output, fileFormat=fileFormat, fileType=fileType)

