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

from icqsol.tools.geometry.icqPrimitiveShape import PrimitiveShape


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

parser.add_argument('--list', dest='list', 
  help='List of options to be passed to the shape constructor of type shape.')

parser.add_argument('--output', dest='output', 
  default='createPrimitiveShape-{0}'.format(tid), 
	help='Output file.')

args = parser.parse_args()

if not args.type:
  print 'ERROR: must specify --type'
  sys.exit(2)

# set the options
for optNameValue in args.options:
  optName, optValue = optNameValue.split('=')
  options[args.type][optName] = eval(optValue)

# set the surface and volume function
surfaceFunctions = []
evalFunction = None
optDic = options[args.type]

if args.type == 'sphere':
  radius = optDic['radius']
  origin = optDic['origin']
  surfaceFunctions = [(lambda u,v: radius*sin(pi*u)*cos(2*pi*v) + origin[0], 
                       lambda u,v: radius*sin(pi*u)*sin(2*pi*v) + origin[1], 
                       lambda u,v: radius*cos(pi*u) + origin[2])]
  radiusSq = radius**2
  def evalFunction(x, y, z):
    xNorm = x - origin[0]
    yNorm = y - origin[1]
    zNorm = z - origin[2]
    return radiusSq - xNorm**2 - yNorm**2 - zNorm**2 > 0

elif args.type == 'cylinder':
  radius = optDic['radius']
  origin = optDic['origin']
  length = optDic['length']
  # the cylindrical side followed by the two end disks
  # u x v normal should point out
  surfaceFunctions = [(lambda u,v: radius*cos(2*pi*u) + origin[0], 
                       lambda u,v: radius*sin(2*pi*u) + origin[1], 
                       lambda u,v: length*(v-0.5) + origin[2]),
                      (lambda u,v: v*radius*cos(2*pi*u) + origin[0],
                       lambda u,v: v*radius*sin(2*pi*u) + origin[1],
                       lambda u,v: (origin[2] - 0.5*length)*numpy.ones(u.shape)),
                      (lambda u,v: u*radius*cos(2*pi*v) + origin[0],
                       lambda u,v: u*radius*sin(2*pi*v) + origin[1],
                       lambda u,v: (origin[2] + 0.5*length)*numpy.ones(u.shape))]
  radiusSq = radius**2
  def evalFunction(x, y, z):
    res = (radiusSq - (x-origin[0])**2 - (y-origin[1])**2) > 0
    res &= (z - origin[2]) < 0.5*length
    return res

elif args.type == 'cone':
  radius = optDic['radius']
  origin = optDic['origin']
  length = optDic['length']
  # origin is the focal point of the cone
  # cone expands in the z direction
  # surface of the cone followed by the end disk
  surfaceFunctions = [(lambda u,v: u*radius*cos(2*pi*v) + origin[0], 
                       lambda u,v: u*radius*sin(2*pi*v) + origin[1], 
                       lambda u,v: u*length + origin[2]),
                      (lambda u,v: v*radius*cos(2*pi*u) + origin[0],
                       lambda u,v: v*radius*sin(2*pi*u) + origin[1],
                       lambda u,v: length + origin[2])]
  radiusSq = radius**2
  def evalFunction(x, y, z):
    xNorm = x - origin[0]
    yNorm = y - origin[1]
    zNorm = z - origin[2]
    res = zNorm > 0
    res &= zNorm < length
    u = zNorm/length
    uRadius = u*radius
    res &= uRadius*uRadius - xNorm*xNorm - yNorm*yNorm > 0
    return res

elif args.type == 'box':
  loBound = numpy.array(optDic['loBound'])
  hiBound = numpy.array(optDic['hiBound'])
  deltas = hiBound - loBound
  # the six faces of the box
  surfaceFunctions = [(lambda u,v: loBound[0], 
                       lambda u,v: loBound[1]+v*deltas[1], 
                       lambda u,v: loBound[2]+u*deltas[2])
                      (lambda u,v: hiBound[0], 
                       lambda u,v: loBound[1]+u*deltas[1], 
                       lambda u,v: loBound[2]+v*deltas[2]), 
                      (lambda u,v: loBound[0]+u*deltas[0], 
                       lambda u,v: loBound[1], lambda u,v: 
                       loBound[2]+v*deltas[2]), 
                      (lambda u,v: loBound[0]+v*deltas[0], 
                       lambda u,v: hiBound[1], 
                       lambda u,v: loBound[2]+u*deltas[2]), 
                      (lambda u,v: loBound[0]+v*deltas[0], 
                       lambda u,v: loBound[1]+u*deltas[1], 
                       lambda u,v: loBound[2]), 
                      (lambda u,v: loBound[0]+u*deltas[0], 
                       lambda u,v: loBound[1]+v*deltas[1], 
                       lambda u,v: hiBound[2])]
  def evalFunction(x, y, z):
    res = x > loBound[0]
    res &= x < hiBound[0]
    res &= y > loBound[1]
    res &= y < hiBound[1]
    res &= z > loBound[2]
    res &= x < hiBound[2]
    return res

else:
  print 'ERROR: unknown shape'
  sys.exit(1)

shp = PrimitiveShape()
shp.setSurfaceFunctions(surfaceFunctions)
shp.setEvaluateFunction(evalFunction)

shp.computeSurfaceMeshes(args.maxArea)

if args.output:
  shp.save(args.output)

