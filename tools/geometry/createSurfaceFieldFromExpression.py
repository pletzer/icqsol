#!/usr/bin/env python

"""
Apply a surface field to a shape
"""

import argparse
import time
import sys
import re
from math import exp, log, sin, cos, pi, tan, acos, asin, atan2, atan, pi

import vtk
from icqsol.tools.geometry.icqShape import Shape

# time stamp
tid = re.sub(r'\.', '', str(time.time()))

parser = argparse.ArgumentParser(description='Apply a surface field to a shape')

parser.add_argument('--input', dest='input', default='',
                    help='Input files (PLY or VTK)')

parser.add_argument('--expression', dest='expression',
                    help='Expression of the coordinates x, y, and z')

parser.add_argument('--ascii', dest='ascii', action='store_true',
                    help='Save data in ASCII format (default is binary)')

parser.add_argument('--output', dest='output',
                    default='createSurfaceFieldFromExpression-{0}.vtk'.format(tid),
                    help='Output file, the suffix (vtk or ply) determines the format.')

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
data = vtk.vtkFloatArray()
data.SetNumberOfComponents(1)
data.SetNumberOfTuples(numPoints)
for i in range(numPoints):
    x, y, z = points.GetPoint(i)
    fieldValue = eval(args.expression)
    data.SetTuple1(i, fieldValue)
pdata.GetPointData().SetScalars(data)

if args.output:
    writer = vtk.vtkPolyDataWriter()
    if args.output.lower().find('.ply') >= 0:
        writer = vtk.vtkPLYWriter()
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

