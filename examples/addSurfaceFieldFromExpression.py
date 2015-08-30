#!/usr/bin/env python

"""
Apply a surface field to a shape
"""

import argparse
import time
import sys
import re
from numpy import linspace
from math import sin, cos, tan, log, exp, pi, asin, acos, atan, atan2, e

import vtk
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
pdata = shape_mgr.shapeToVTKPolyData(shp)
points = pdata.GetPoints()
numPoints = points.GetNumberOfPoints()
data = vtk.vtkDoubleArray()
data.SetName(args.name)
times = [0.0]
if args.times:
    times = eval(args.times)
numTimes = len(times)
data.SetNumberOfComponents(numTimes)
data.SetNumberOfTuples(numPoints)
for i in range(numPoints):
    x, y, z = points.GetPoint(i)
    for j in range(numTimes):
        t = times[j]
        fieldValue = eval(args.expression)
        data.SetComponent(i, j, fieldValue)
pdata.GetPointData().SetScalars(data)

if args.output:
    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(args.output)
    if args.ascii:
        writer.SetFileTypeToASCII()
    else:
        writer.SetFileTypeToBinary()
    if vtk.VTK_MAJOR_VERSION >= 6:
        writer.SetInputData(pdata)
    else:
        writer.SetInput(pdata)
    writer.Write()
    writer.Update()
