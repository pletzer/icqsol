#!/usr/bin/env python

"""
Create a primitive shape of a certain type
"""

import argparse
import time
import sys
import re
import numpy
from numpy import cos, sin, pi

from icqsol.tools.geometry.icqPrimitiveShape import Box, Cone, Cylinder, Sphere


options = {
	'sphere': {
    'radius': 1.0,
    'origin': numpy.array([0.,0., 0.]),
  },
  'cylinder': {
    'radius': 1.0,
    'origin': numpy.array([0.,0., 0.]),
    'length': 1.0,
  },
  'cone': {
    'origin': numpy.array([0.,0., 0.]),
    'radius': 1.0,
    'length': 1.0,
  },
  'box': {
    'loBound': numpy.array([0., 0., 0.]),
    'hiBound': numpy.array([1., 1., 1.]),
  }, 
}

# time stamp
tid = re.sub(r'\.', '', str(time.time()))

parser = argparse.ArgumentParser(description='Create primitive shape.')

parser.add_argument('--maxArea', dest='maxArea', type=float, default=0.01,
  help='Max cell area estimate for surface discretization')

parser.add_argument('--type', dest='type',
	help='Type, currently either "sphere", "cylinder", "cone", or "box"')

parser.add_argument('--options', dest='options', nargs='*', default=[],
  help='Options to be passed to the shape constructor.')

parser.add_argument('--list', dest='list', action='store_true',
  help='List of options to be passed to the shape constructor of type shape.')

parser.add_argument('--check', dest='check', action='store_true',
  help='Check validity of supplied options, by default errors will ignored.')


parser.add_argument('--output', dest='output', 
  default='createPrimitiveShape-{0}.vtk'.format(tid), 
	help='Output file.')

args = parser.parse_args()

if not args.type:
  print 'ERROR: must specify --type'
  sys.exit(2)

# set the options
for optNameValue in args.options:
  optName, optValue = optNameValue.split('=')
  if args.check and optName not in options[args.type]:
    print 'ERROR: invalid option named "{0}"'.format(optName)
    sys.exit(3)
  options[args.type][optName] = eval(optValue)

# set the surface and volume function
surfaceFunctions = []
evalFunction = None
optDic = options[args.type]

if args.type == 'sphere':
  radius = optDic['radius']
  origin = optDic['origin']
  shp = Sphere( radius, origin )

elif args.type == 'cylinder':
  radius = optDic['radius']
  origin = optDic['origin']
  length = optDic['length']
  # the cylindrical side followed by the two end disks
  # u x v normal should point out
  shp = Cylinder( length, radius, origin )

elif args.type == 'cone':
  radius = optDic['radius']
  origin = optDic['origin']
  length = optDic['length']
  # origin is the focal point of the cone
  # cone expands in the z direction
  # surface of the cone followed by the end disk
  shp = Cone( length, radius, origin )

elif args.type == 'box':
  loBound = numpy.array(optDic['loBound'])
  hiBound = numpy.array(optDic['hiBound'])
  shp = Box( loBound, hiBound )

else:
  print 'ERROR: unknown shape'
  sys.exit(1)

# display list of options (if asking for it)
if args.list:
  for optName, optVal in options[args.type].items():
    print '{:>10} --> {:>20}'.format(optName, optVal)

shp.computeSurfaceMeshes(args.maxArea)

if args.output:
  fileFormat = 'vtk'
  fileType = 'binary'
  if args.output.lower().find('.ply') > 0:
    fileFormat = 'ply'
  shp.save(args.output, fileFormat=fileFormat, fileType=fileType)

