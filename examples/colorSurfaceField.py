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

# TODO: Enhance this to read the filoe and discover the vtk_dataset_type.
shape_mgr = ShapeManager(file_format='vtk', vtk_dataset_type='POLYDATA')
pDataInput = shape_mgr.loadAsVtkPolyData(args.input)
pDataColored = shape_mgr.colorSurfaceField(pDataInput, args.colormap, field_name=args.name)

if args.output:
    if args.ascii:
        file_type = 'ascii'
    else:
        file_type = 'binary'
    shape_mgr.saveVtkPolyData(vtk_poly_data=pDataColored, file_name=args.output, file_type=file_type)
