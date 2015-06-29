#!/usr/bin/env python

"""
Test creation of simple ring shape object
@author pletzer@psu.edu
"""

from icqsol.tools.geometry.icqCylinder import Cylinder

c1 = Cylinder(radius=0.4, origin=(0., 0., 0.), length=0.5, n_z=12)
c2 = Cylinder(radius=0.2, origin=(0., 0., -0.1), length=1.0, n_z=30)
geom = c1 + c2
geom.show() #filename='testRing.png')
