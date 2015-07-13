#!/usr/bin/env python

"""
Test creation of a bolt object
@author pletzer@psu.edu
"""

from icqsol.tools.geometry.icqShape import Shape
from icqsol.tools.geometry.icqCone import Cone
from icqsol.tools.geometry.icqCylinder import Cylinder
from icqsol.tools.geometry.icqBox import Box

head = Cone(radius=0.7, origin = (0., 0., 3.7), lengths = (0., 0., 1.0))
shaft = Cylinder(radius=0.5, origin=(0., 0., 2.0), lengths = (0., 0., 4.0))
notch1 = Box(origin=(-1., -0.1, 4.4), lengths=(2., 0.2, 0.6))
notch2 = Box(origin=(-0.1, -1., 4.4), lengths=(0.2, 2., 0.6))

geom = head + shaft - notch1 - notch2

geom.rotate(axis=(0., 0., 1.), angleDeg=45.0)

geom.show(filename='testBolt.png')
