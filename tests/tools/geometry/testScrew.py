#!/usr/bin/env python

"""
Test simple ring object
@author pletzer@psu.edu
"""

from icqsol.tools.geometry.icqGeometry import Geometry
from icqsol.tools.geometry.icqCone import Cone
from icqsol.tools.geometry.icqCylinder import Cylinder

head = Cone(angle = 75.0, origin = (0., 0., 1.), length=0.1)
shaft = Cylinder(radius=0.1, origin=(0.1, 0.5, 0.5), length=1.0)

geom = Geometry()
geom += head
#geom += shaft
geom.computeBoundarySurface(100, 100, 100)
geom.show()
