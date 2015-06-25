#!/usr/bin/env python

"""
Rotate shape
"""

import argparse
import time
import sys
import re
import numpy

from icqsol.tools.geometry.icqBaseShape import BaseShape

# time stamp
tid = re.sub(r'\.', '', str(time.time()))

parser = argparse.ArgumentParser(description='Translate shape.')

parser.add_argument('--input', dest='input', default='',
  help='List of input files (PLY or VTK)')

parser.add_argument('--angle', dest='angle', type=float, default=0.0,
	help='Specify rotation angle in degrees')

parser.add_argument('--origin', dest='origin', default="0.,0.,0.",
	help='Specify rotation point as three floating point numbers')

parser.add_argument('--axis', dest='axis', default="0.,0.,1.",
	help='Specify rotation axis as three floating point numbers')

parser.add_argument('--output', dest='output', 
  default='createCompositeShape-{0}.vtk'.format(tid), 
	help='Output file.')

args = parser.parse_args()

if not args.input:
  print 'ERROR: must specify one input file with --input <file>'
  sys.exit(3)

shp = BaseShape()
origin = eval(args.origin)
axis = eval(args.axis)
shp.rotate(origin=origin, angleDeg=args.angle, axis=axis)
if args.output:
  shp.save(args.output)

