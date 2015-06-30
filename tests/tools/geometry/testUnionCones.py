#!/usr/bin/env python

"""
Test union operation
@author pletzer@psu.edu
"""

from icqsol.tools.geometry.icqCone import Cone

s1 = Cone(radius=1.0, origin=(0., 0., 0.5), length=0.5,
          n_rho=100, n_theta=32, n_z=100)
s2 = Cone(radius=0.5, origin=(0., 0., 0.), length=1.0,
          n_rho=100, n_theta=32, n_z=100)
geom = s1 + s2
geom.show() #filename='testUnionCones.png')
