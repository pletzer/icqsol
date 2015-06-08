#!/usr/bin/env python

"""
Test creation of a screw object
@author pletzer@psu.edu
"""

from icqsol.tools.geometry.icqShape import Shape
from icqsol.tools.geometry.icqCone import Cone
from icqsol.tools.geometry.icqCylinder import Cylinder
from icqsol.tools.geometry.icqBox import Box

head = Cone(angle = 60.0, origin = (0., 0., 3.7), length = 1.0)
shaft = Cylinder(radius=0.5, origin=(0., 0., 2.0), length = 4.0)
notch1 = Box(loBound = (-1., -0.1, 4.4), hiBound = (1., 0.1, 5.0))
notch2 = Box(loBound = (-0.1, -1., 4.4), hiBound = (0.1, 1., 5.0))

geom = head + shaft - notch1 - notch2

geom.computeBoundarySurface(100, 100, 100)
geom.show(filename='testScrew.png')
