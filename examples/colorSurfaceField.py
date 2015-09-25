#!/usr/bin/env python

"""
Color a surface field
"""

import argparse
import time
import os
import re
import sys

from icqsol.shapes.icqShapeManager import ShapeManager
from icqsol import util

# time stamp
tid = re.sub(r'\.', '', str(time.time()))

parser = argparse.ArgumentParser(description='Color surface field')

parser.add_argument('--input', dest='input', default='',
                    help='VTK input file')

parser.add_argument('--colormap', dest='colormap', default='hot',
                    help='Colormap ("gnu", "hot", "cold", or "blackbody")')

parser.add_argument('--name', dest='name', default='',
                    help='Set the name of the field')

parser.add_argument('--component', dest='component', type=int, default=0,
                    help='Set the component of the field')

parser.add_argument('--ascii', dest='ascii', action='store_true',
                    help='Save data in ASCII format (default is binary)')

parser.add_argument('--display', dest='display', action='store_true',
                    help='Display colored geometry')

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

file_format = util.getFileFormat(args.input)

if file_format != util.VTK_FORMAT:
    print 'ERROR: file {} must be VTK format'.format(args.input)
    sys.exit(2)

vtk_dataset_type = util.getVtkDatasetType(args.input)

if vtk_dataset_type not in util.VTK_DATASET_TYPES:
    print 'ERROR: invalid VTK dataset type {}'.format(vtk_dataset_type)
    sys.exit(2)

shape_mgr = ShapeManager(file_format=util.VTK_FORMAT, vtk_dataset_type=vtk_dataset_type)
pDataInput = shape_mgr.loadAsVtkPolyData(args.input)
pDataColored = shape_mgr.colorSurfaceField(pDataInput, args.colormap,
                                           field_name=args.name,
                                           field_component=args.component)
if args.display:
    shape_mgr.showVtkPolyData(pDataColored)

if args.output:
    if args.ascii:
        file_type = util.ASCII
    else:
        file_type = util.BINARY
    shape_mgr.saveVtkPolyData(vtk_poly_data=pDataColored, file_name=args.output, file_type=file_type)
