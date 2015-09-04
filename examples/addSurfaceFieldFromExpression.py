#!/usr/bin/env python

"""
Apply a surface field to a shape
"""
import argparse
import time
import sys
import re

from icqsol.shapes.icqShapeManager import ShapeManager

# time stamp
tid = re.sub(r'\.', '', str(time.time()))

description = 'Apply a surface field to a shape'
parser = argparse.ArgumentParser(description=description)

parser.add_argument('--input', dest='input', default='',
                    help='Input file (PLY or VTK)')

parser.add_argument('--expression', dest='expression', default='scalar_field',
                    help='Expression of x, y, z, and t')

parser.add_argument('--name', dest='name',
                    help='Set the name of the field')

parser.add_argument('--times', dest='times', default='',
                    help='Comma separated list of time values')

parser.add_argument('--ascii', dest='ascii', action='store_true',
                    help='Save data in ASCII format (default is binary)')

parser.add_argument('--output', dest='output',
                    default='addSurfaceFieldFromExpression-{0}.vtk'.format(tid),
                    help='VTK Output file.')

args = parser.parse_args()

if not args.expression:
    print 'ERROR: must specify --expression <expression>'
    sys.exit(2)

if not args.input:
    print 'ERROR: must specify input file: --input <file>'
    sys.exit(3)

# make sure the field name contains no spaces
args.name = re.sub('\s', '_', args.name)

shape_mgr = ShapeManager()
shp = shape_mgr.load(args.input)
times = [0.0]
if args.times:
    times = eval(args.times)
pdata = shape_mgr.addSurfaceFieldFromExpression(shp, args.name, args.expression, times)

if args.output:
    if args.ascii:
        file_type = 'ascii'
    else:
        file_type = 'binary'
    shape_mgr.save(file_name=args.output, file_format='vtk', file_type=file_type, vtk_poly_data=pdata)
