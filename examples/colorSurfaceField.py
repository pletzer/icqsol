#!/usr/bin/env python

"""
Color a surface field
"""

import argparse
import time
import os
import re
import sys

import numpy
import vtk
from icqsol.tools.color.icqColorMap import ColorMap

# time stamp
tid = re.sub(r'\.', '', str(time.time()))

parser = argparse.ArgumentParser(description='Color surface field')

parser.add_argument('--input', dest='input', default='',
                    help='VTK input file')

parser.add_argument('--colormap', dest='colormap', default='hot',
                    help='Colormap ("hot", "cold", or "blackbody")')

parser.add_argument('--name', dest='name', default='',
                    help='Set the name of the field')

parser.add_argument('--ascii', dest='ascii', action='store_true',
                    help='Save data in ASCII format (default is binary)')

parser.add_argument('--output', dest='output',
                    default='colorSurfaceField-{0}.vtk'.format(tid),
                    help='VTK Output file.')

args = parser.parse_args()

if not args.input:
    print 'ERROR: must specify input file: --input <file>'
    sys.exit(3)

if not os.path.exists(args.input):
    print 'ERROR: file {} does not exist'.format(args.input)
    sys.exit(2)

# read the data
reader = vtk.vtkPolyDataReader()
reader.SetFileName(args.input)
reader.Update()

pdataInput = reader.GetOutput()
pointData = pdataInput.GetPointData()
numComps = pointData.GetNumberOfComponents()
numPoints = pdataInput.GetPoints().GetNumberOfPoints()
numArrays = pointData.GetNumberOfArrays()

if not args.name and numArrays > 1: 
    print 'ERROR: more than one array, must specify --name NAME'
    sys.exit(1)

# min/max field values
array = pointData.GetScalars(args.name)
fmin, fmax = array.GetRange()

# build the colormap
colormap = ColorMap(fmin, fmax)
colrMethod = eval('colormap.' + args.colormap)

# create another poly data with the same points/cells
pdataOutput = vtk.vtkPolyData()
pdataOutput.SetPoints(pdataInput.GetPoints())
pdataOutput.SetPolys(pdataInput.GetPolys())

# color the points
rArray = vtk.vtkUnsignedCharArray()
gArray = vtk.vtkUnsignedCharArray()
bArray = vtk.vtkUnsignedCharArray()
rArray.SetNumberOfComponents(numComps)
gArray.SetNumberOfComponents(numComps)
bArray.SetNumberOfComponents(numComps)
rArray.SetNumberOfTuples(numPoints)
gArray.SetNumberOfTuples(numPoints)
bArray.SetNumberOfTuples(numPoints)
rs = numpy.zeros((numComps,), numpy.uint8)
gs = numpy.zeros((numComps,), numpy.uint8)
bs = numpy.zeros((numComps,), numpy.uint8)

for i in range(numPoints):
    fs = array.GetTuple(i)
    for j in range(numComps):
        rs[j], gs[j], bs[j] = colrMethod(fs[j])
    rArray.SetTuple(i, rs)
    gArray.SetTuple(i, gs)
    bArray.SetTuple(i, bs)

pointDataOut = pdataOutput.GetPointData()
for arr in rArray, gArray, bArray:
    pointDataOut.AddArray(arr)
pointDataOut.GetArray(0).SetName('red')
pointDataOut.GetArray(1).SetName('green')
pointDataOut.GetArray(2).SetName('blue')

if args.output:
    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(args.output)
    if args.ascii:
        writer.SetFileTypeToASCII()
    else:
        writer.SetFileTypeToBinary()
    if vtk.VTK_MAJOR_VERSION >= 6:
        writer.SetInputData(pdataOutput)
    else:
        writer.SetInput(pdataOutput)
    writer.Write()
    writer.Update()
