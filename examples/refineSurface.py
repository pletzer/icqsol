#!/usr/bin/env python

"""
Refine a surface by adding points along cell edges
"""
from __future__ import print_function
import argparse
import time
import sys
import re

from icqsol.shapes.icqShapeManager import ShapeManager
from icqsol import util

# time stamp
tid = re.sub(r'\.', '', str(time.time()))

description = 'Refine a surface by adding points along cell edges'
parser = argparse.ArgumentParser(description=description)

parser.add_argument('--input', dest='input', default='',
                    help='Input file (PLY or VTK)')

parser.add_argument('--maxedge', dest='maxedge', type=float,
                    default=float('inf'),
                    help='Maximum edge length')

parser.add_argument('--ascii', dest='ascii', action='store_true',
                    help='Save data in ASCII format (default is binary)')

parser.add_argument('--output', dest='output',
                    default='addSurfaceFieldFromExpressionToShape-{0}.vtk'.format(tid),
                    help='VTK Output file.')

args = parser.parse_args()

if not args.input:
    print('ERROR: must specify input file: --input <file>')
    sys.exit(3)

# Get the format of the input - either vtk or ply.
file_format = util.getFileFormat(args.input)

if file_format == util.PLY_FORMAT:
    shape_mgr = ShapeManager(file_format=util.PLY_FORMAT)
else:
    # We have a VTK file, so Get the dataset type.
    vtk_dataset_type = util.getVtkDatasetType(args.input)
    shape_mgr = ShapeManager(file_format=util.VTK_FORMAT,
                             vtk_dataset_type=vtk_dataset_type)

pdata_input = shape_mgr.loadAsVtkPolyData(args.input)

pdata = shape_mgr.refineVtkPolyData(pdata_input, max_edge_length=args.maxedge)

if args.output:
    # Always produce VTK POLYDATA.
    shape_mgr.setWriter(file_format=util.VTK_FORMAT,
                        vtk_dataset_type=util.POLYDATA)
    if args.ascii:
        file_type = util.ASCII
    else:
        file_type = util.BINARY
    shape_mgr.saveVtkPolyData(vtk_poly_data=pdata,
                              file_name=args.output,
                              file_type=file_type)
