#!/usr/bin/env python

"""
Test creation of simple ring shape object
@author pletzer@psu.edu
"""

from icqsol.tools.geometry.icqCylinder import Cylinder

c1 = Cylinder(radius=0.4, origin=(0.5, 0.5, 0.5), length=0.5)
c2 = Cylinder(radius=0.2, origin=(0.5, 0.5, 0.5), length=0.5)

geom = c1 - c2
geom.computeBoundarySurface(100, 100, 100)
geom.show()
