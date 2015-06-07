#!/usr/bin/env python

"""
Geometry of a screw
@author pletzer@psu.edu
"""

from icqsol.tools.geometry.icqGeometry import Geometry
from icqsol.tools.geometry.icqCone import Cone
from icqsol.tools.geometry.icqCylinder import Cylinder
from icqsol.tools.geometry.icqBox import Box

head = Cone(angle = 60.0, origin = (0., 0., 3.7), length = 1.0)
shaft = Cylinder(radius=0.5, origin=(0., 0., 2.0), length = 4.0)
notch1 = Box(bxLo = (-1., -0.1, 4.4), bxHi = (1., 0.1, 5.0))
notch2 = Box(bxLo = (-0.1, -1., 4.4), bxHi = (0.1, 1., 5.0))

geom = Geometry()
geom += head
geom -= notch1
geom -= notch2
geom += shaft

geom.computeBoundarySurface(100, 100, 100)
geom.show()
