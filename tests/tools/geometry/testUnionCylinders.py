#!/usr/bin/env python

"""
Test union operation
@author pletzer@psu.edu
"""

from icqsol.tools.geometry.icqCylinder import Cylinder

s1 = Cylinder(radius=1.0, origin=(0., 0., 0.5), length=0.5)
s2 = Cylinder(radius=0.5, origin=(0., 0., 0.), length=1.0)
geom = s1 + s2
geom.show() #filename='testUnionCylinders.png')
