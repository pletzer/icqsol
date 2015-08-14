#!/usr/bin/env python

"""
Apply a surface field to a shape
"""

import argparse
import time
import sys
import re
from math import exp, log, sin, cos, pi, tan, acos, asin, atan2, atan, pi
from numpy import linspace

import vtk
from icqsol.tools.geometry.icqShape import Shape

# time stamp
tid = re.sub(r'\.', '', str(time.time()))

parser = argparse.ArgumentParser(description='Apply a surface field to a shape')

parser.add_argument('--input', dest='input', default='',
                    help='Input files (PLY or VTK)')

parser.add_argument('--expression', dest='expression',
                    help='Expression of x, y, z, and t')

parser.add_argument('--times', dest='times', default='0.0',
                    help='Comma separated list of time values')

parser.add_argument('--ascii', dest='ascii', action='store_true',
                    help='Save data in ASCII format (default is binary)')

parser.add_argument('--output', dest='output',
                    default='createSurfaceFieldFromExpression-{0}.vtk'.format(tid),
                    help='VTK Output file.')

args = parser.parse_args() 

if not args.expression:
    print 'ERROR: must specify --expression <expression>'
    sys.exit(2)

if not args.input:
    print 'ERROR: must specify input file: --input <file>'
    sys.exit(3)

shp = Shape.load(args.input)
pdata = shp.toVTKPolyData()
points = pdata.GetPoints()
numPoints = points.GetNumberOfPoints()
data = vtk.vtkDoubleArray()
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
