#!/usr/bin/env python
"""
Create a primitive shapes
"""
import argparse
import time
import sys
import re
import numpy
from icqsol.shapes.icqShapeManager import ShapeManager

options = {
    'sphere': {
        'radius': 1.0,
        'origin': numpy.array([0., 0., 0.]),
        'n_theta': 8,
        'n_phi': 16,
    },
    'cylinder': {
        'radius': 1.0,
        'origin': numpy.array([0.0, 0.0, 0.0]),
        'lengths': numpy.array([1.0, 0.0, 0.0]),
        'n_theta': 16,
    },
    'cone': {
        'origin': numpy.array([0.0, 0.0, 0.0]),
        'radius': 1.0,
        'lengths': numpy.array([1.0, 0.0, 0.0]),
        'n_theta': 16,
    },
    'box': {
        'origin': numpy.array([0.0, 0.0, 0.0]),
        'lengths': numpy.array([1.0, 1.0, 1.0]),
    },
}

# time stamp
tid = re.sub(r'\.', '', str(time.time()))

parser = argparse.ArgumentParser(description='Create primitive shape.')

parser.add_argument(
    '--type', dest='type',
    help='Type, currently either "sphere", "cylinder", "cone", or "box"')

parser.add_argument('--options', dest='options', nargs='*', default=[],
                    help='Options to be passed to the shape constructor.')

parser.add_argument(
    '--list', dest='list', action='store_true',
    help='List of options to be passed to the shape constructor.')

parser.add_argument(
    '--check', dest='check', action='store_true',
    help='Check validity of supplied options, by default errors will ignored.')

parser.add_argument('--ascii', dest='ascii', action='store_true',
                    help='Save data in ASCII format (default is binary)')

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

fileFormat = 'vtk'
fileType = 'binary'
if args.ascii:
    fileType = 'ascii'
if args.output.lower().find('.ply') > 0:
    fileFormat = 'ply'

if fileFormat == 'vtk':
    # TODO: Enhance this to read the filoe and discover the vtk_dataset_type.
    shape_mgr = ShapeManager(file_format=fileFormat, vtk_dataset_type='POLYDATA')
else:
    shape_mgr = ShapeManager(file_format=fileFormat)

# set the surface and volume function
surfaceFunctions = []
evalFunction = None
optDic = options[args.type]

if args.type == 'sphere':
    radius = optDic['radius']
    origin = optDic['origin']
    n_theta = optDic['n_theta']
    n_phi = optDic['n_phi']
    shp = shape_mgr.createShape('sphere', radius=radius, origin=origin,
                                n_theta=n_theta, n_phi=n_phi)

elif args.type == 'cylinder':
    radius = optDic['radius']
    origin = optDic['origin']
    lengths = optDic['lengths']
    n_theta = optDic['n_theta']
    # the cylindrical side followed by the two end disks
    shp = shape_mgr.createShape('cylinder', origin=origin, radius=radius,
                                lengths=lengths, n_theta=n_theta)

elif args.type == 'cone':
    radius = optDic['radius']
    origin = optDic['origin']
    lengths = optDic['lengths']
    n_theta = optDic['n_theta']
    # origin is axis location where radius is max
    shp = shape_mgr.createShape('cone', origin=origin, radius=radius,
                                lengths=lengths, n_theta=n_theta)

elif args.type == 'box':
    origin = numpy.array(optDic['origin'])
    lengths = numpy.array(optDic['lengths'])
    shp = shape_mgr.createShape('box', origin=origin, lengths=lengths)

else:
    print 'ERROR: unknown shape'
    sys.exit(1)

# display list of options (if asking for it)
if args.list:
    for optName, optVal in options[args.type].items():
        print '{:>10} --> {:>20}'.format(optName, optVal)

if args.output:
    shape_mgr.saveShape(shape=shp, file_name=args.output, file_type=fileType)
