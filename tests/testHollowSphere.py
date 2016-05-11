#!/usr/bin/env python

"""
Test creation of sphere
@author alexander@gokliya.net
"""

from __future__ import print_function
import argparse
from icqsol.shapes.icqShapeManager import ShapeManager
from icqsol import util

parser = argparse.ArgumentParser(description='Create sphere')

parser.add_argument('--output', type=str, dest='output', default='',
                    help='Set output file, suffix should be .vtk or .ply')
parser.add_argument('--n_theta', type=int, dest='n_theta', default=8,
                    help='Set number of poloidal cells')
parser.add_argument('--n_phi', type=int, dest='n_phi', default=4,
                    help='Set number of azimuthal cells')
args = parser.parse_args()

# Get the format of the input - either vtk or ply.
file_format = util.getFileFormat(args.output)

if file_format == util.VTK_FORMAT:
    shape_mgr = ShapeManager(file_format=file_format, vtk_dataset_type='POLYDATA')
else:
    shape_mgr = ShapeManager(file_format=file_format)

bigSphere = shape_mgr.createShape('sphere',
                                  origin=(0., 0., 0.),
                                  radius=2.0,
                                  n_theta=args.n_theta,
                                  n_phi=args.n_phi)
                                  
smallSphere = shape_mgr.createShape('sphere',
                                    origin=(0., 0., 0.),
                                    radius=1.0,
                                    n_theta=2*args.n_theta,
                                    n_phi=2*args.n_phi)
                                    
geom = bigSphere - smallSphere

if args.output:
    shape_mgr.saveShape(shape=geom, file_name=args.output, file_type=util.ASCII)
else:
    shape_mgr.showShape(geom)
