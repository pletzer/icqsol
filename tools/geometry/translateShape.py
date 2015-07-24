#!/usr/bin/env python

"""
Translate shape
"""

import argparse
import time
import sys
import re

from icqsol.tools.geometry.icqBaseShape import BaseShape

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

shp = BaseShape()
shp.load(args.input)
transVec = eval(args.translation)
shp.translate(transVec)
if args.output:
    fileFormat = 'vtk'
    fileType = 'binary'
    if args.ascii:
        fileType = 'ascii'
    if args.output.lower().find('.ply') >= 0:
        fileFormat = 'ply'
    shp.save(args.output, fileFormat=fileFormat, fileType=fileType)
