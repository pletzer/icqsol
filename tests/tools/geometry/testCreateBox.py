#!/usr/bin/env python

"""
Test creation of box
@author pletzer@psu.edu
"""

import argparse
from icqsol.tools.geometry.icqBox import Box

parser = argparse.ArgumentParser(description='Create box')

parser.add_argument('--output', type=str, dest='output', default='',
                   help='Set output file, the suffix (either .vtk or .ply) determines the format')
parser.add_argument('--lengths', type=str, dest='lengths', default="1.,1.,1.",
                   help='Set box lengths')
parser.add_argument('--origin', type=str, dest='origin', default="0.,0.,0.",
                   help='Set origin, comma separated list of three floats')
parser.add_argument('--n_x', type=int, dest='n_x', default=1,
                   help='Set number of x cells')
parser.add_argument('--n_y', type=int, dest='n_y', default=1,
                   help='Set number of y cells')
parser.add_argument('--n_z', type=int, dest='n_z', default=1,
                   help='Set number of z cells')
args = parser.parse_args()

s = Box(origin=eval(args.origin), lengths=eval(args.lengths),
      n_x=args.n_x, n_y=args.n_y, n_z=args.n_z)

if args.output:
	file_format = 'vtk'
        file_type = 'ascii'
        if args.output.find('.ply') > 0:
    		file_format = 'ply'
	s.save(args.output, file_format, file_type)
else:
	s.show()
