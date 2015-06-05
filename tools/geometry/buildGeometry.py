#!/usr/bin/nev python

import argparse
import time
import re

from icqsol.tools.geometry.icqGeometry import Geometry
from icqsol.tools.geometry.icqBox import Box
from icqsol.tools.geometry.icqSphere import Sphere
from icqsol.tools.geometry.icqCylinder import Cylinder

# time stamp
tid = re.sub(r'\.', '', time.time())

parser = argparse.ArgumentParser(description='Build geometry.')

parser.add_argument('--create', dest='createExprs', nargs='+', 
	help='Create shape, e.g. "s = Cylinder(radius=0.1, origin=[0., 0., 0.])"')

parser.add_argument('--assemble', dest='assembleExpr', nargs=1, 
	help='Assemble geometry, e.g. "(s + s2)*s3 - s4"')

parser.add_argument('--output', dest='output', default='buildGeometry-{0}'.format(tid), 
	help='Output file.')

args = parser.parse_args()

# instantiate the shapes
for shapeExpr in args.createExprs:
	exec(shapeExpr)

# apply transformations
# TO DO 

# assemble the shapes
geom = Geometry()
geom.apply(assembleExpr)

# save
geom

