#!/usr/bin/nev python

import argparse
import time
import re

from icqsol.tools.geometry.icqBox import Box
from icqsol.tools.geometry.icqSphere import Sphere
from icqsol.tools.geometry.icqCylinder import Cylinder
from icqsol.tools.common.icqPLYWriter import PLYWriter

# time stamp
tid = re.sub(r'\.', '', str(time.time()))

parser = argparse.ArgumentParser(description='Build geometry.')

parser.add_argument('--create', dest='createExprs', nargs='+', 
	help='Create shape, e.g. "c = Cylinder(radius=0.1, origin=[0., 0., 0.], length=1.0)"')

parser.add_argument('--transforms', dest='transforms', nargs='*', 
  help='Apply transformations to each object in sequence, e.g. c.rotateX(25.0)')

parser.add_argument('--assemble', dest='assembleExpr',
	help='Assemble geometry, e.g. "(c1 + s2)*s3 - s4"')

parser.add_argument('--output', dest='output', default='buildGeometry-{0}'.format(tid), 
	help='Output file.')

args = parser.parse_args()

# instantiate the shapes
for shapeExpr in args.createExprs:
  exec(shapeExpr)

# apply transformations
if args.transforms:
  for transf in args.transforms:
    exec(transf)

# assemble the shapes
geom = eval(args.assembleExpr)

# save
if args.output:

  geom.computeBoundarySurface(100, 100, 100)
  data = geom.getBoundarySurface()

  pw = PLYWriter(args.output)
  pw.setVertices(data['points'])
  pw.setTriangles(data['cells'])
  pw.write()




