#!/usr/bin/env python

"""
Test intersection operation
@author pletzer@psu.edu
"""

from icqsol.tools.geometry.icqCylinder import Cylinder

s1 = Cylinder(radius=1.0, origin=(0., 0., 0.5), length=0.5)
# need to make the cylinder to subtract a little longer to avoid
# floating point issues
s2 = Cylinder(radius=0.5, origin=(0., 0., 0.), length=1.1)
geom = s1 - s2
geom.show() #filename='testSubtractionCylinders.png')
