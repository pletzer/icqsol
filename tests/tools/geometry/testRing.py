#!/usr/bin/env python

"""
Test creation of simple ring object
@author pletzer@psu.edu
"""

from icqsol.tools.geometry.icqCylinder import Cylinder

c1 = Cylinder(radius=0.4, origin=(0., 0., 0.), 
              lengths=(0.5, 0., 0.))
c2 = Cylinder(radius=0.2, origin=(0., 0., -0.1), 
              lengths=(1.0, 0., 0.))
geom = c1 - c2
geom.show() #filename='testRing.png')
