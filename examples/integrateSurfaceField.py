#!/usr/bin/env python

"""
Integrate a surface field
"""

from __future__ import print_function
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

parser.add_argument('--name', dest='name', default='',
                    help='Set the name of the field')

parser.add_argument('--component', dest='component', type=int, default=0,
                    help='Set the component of the field')

parser.add_argument('--ascii', dest='ascii', action='store_true',
                    help='Save data in ASCII format (default is binary)')

args = parser.parse_args()

if not args.input:
    print('ERROR: must specify input file: --input <file>')
    sys.exit(3)

if not os.path.exists(args.input):
    print('ERROR: file {0} does not exist'.format(args.input))
    sys.exit(2)

file_format = util.getFileFormat(args.input)

if file_format != util.VTK_FORMAT:
    print('ERROR: file {0} must be VTK format'.format(args.input))
    sys.exit(2)

vtk_dataset_type = util.getVtkDatasetType(args.input)

if vtk_dataset_type not in util.VTK_DATASET_TYPES:
    print('ERROR: invalid VTK dataset type {0}'.format(vtk_dataset_type))
    sys.exit(2)

shape_mgr = ShapeManager(file_format=util.VTK_FORMAT, vtk_dataset_type=vtk_dataset_type)
pDataInput = shape_mgr.loadAsVtkPolyData(args.input)
integral = shape_mgr.integrateSurfaceField(pDataInput,
                                           field_name=args.name,
                                           field_component=args.component)

print('integral = {0}'.format(integral))
