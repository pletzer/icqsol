#!/usr/bin/env python

"""
Refine a shape by subdividing the polygons
"""

import argparse
import time
import sys
import re

from icqsol.tools.geometry.icqShape import Shape
import vtk

# time stamp
tid = re.sub(r'\.', '', str(time.time()))

parser = argparse.ArgumentParser(description='Refine a shape.')

parser.add_argument('--input', dest='input', default='',
                    help='Input file (PLY or VTK)')

parser.add_argument('--refine', dest='refine', type=int,
                    help='Number of subdivisions')

parser.add_argument('--ascii', dest='ascii', action='store_true',
                    help='Save data in ASCII format (default is binary)')

parser.add_argument('--output', dest='output',
                    default='createCompositeShape-{0}.vtk'.format(tid),
                    help='Output file.')

args = parser.parse_args()

if args.refine <= 0:
    print 'ERROR: refine must be a positive, integer number: --refine #'
    sys.exit(2)

if not args.input:
    print 'ERROR: must specify input file: --input <file>'
    sys.exit(3)

shp = Shape.load(args.input)
pdata = shp.toVTKPolyData()
densifyFilter = vtk.vtkDensifyPolyData()
densifyFilter.SetNumberOfSubdivisions(args.refine)
if vtk.VTK_MAJOR_VERSION >= 6:
    densifyFilter.SetInputData(pdata)
else:
    densifyFilter.SetInput(pdata)
densifyFilter.Update()

outPdata = densifyFilter.GetOutput()

if args.output:
    writer = vtk.vtkPolyDataWriter()
    if args.output.lower().find('.ply') >= 0:
        writer = vtk.vtk.PolyDataWriter()
    writer.SetFileName(args.output)
    writer.SetFileTypeToBinary()
    if args.ascii:
        writer.SetFileTypeToASCII()
    if vtk.VTK_MAJOR_VERSION >= 6:
        writer.SetInputData(outPdata)
    else:
        writer.SetInput(outPdata)
    writer.Write()
    writer.Update()