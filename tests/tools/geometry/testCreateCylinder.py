#!/usr/bin/env python

"""
Test creation of cylinders
@author pletzer@psu.edu
"""

import argparse
from icqsol.tools.geometry.icqCylinder import Cylinder

parser = argparse.ArgumentParser(description='Create cylinder')

parser.add_argument('--output', type=str, dest='output', default='',
                   help='Set output file')
parser.add_argument('--radius', type=float, dest='radius', default=1.0,
                   help='Set radius')
parser.add_argument('--length', type=float, dest='length', default=1.0,
                   help='Set length')
parser.add_argument('--origin', type=str, dest='origin', default="0.,0.,0.",
                   help='Set origin, comma separated list of three floats')
parser.add_argument('--n_rho', type=int, dest='n_rho', default=1,
                   help='Set number of radial cells')
parser.add_argument('--n_theta', type=int, dest='n_theta', default=32,
                   help='Set number of poloidal cells')
parser.add_argument('--n_z', type=int, dest='n_z', default=1,
                   help='Set number of z cells')
args = parser.parse_args()

s = Cylinder(radius=args.radius, origin=eval(args.origin), length=args.length,
             n_rho=args.n_rho, n_theta=args.n_theta, n_z=args.n_z)

if args.output:
	file_format = 'vtk'
        file_type = 'ascii'
        if args.output.find('.ply') > 0:
    		file_format = 'ply'
	s.save(args.output, file_format, file_type)
else:
	s.show()
