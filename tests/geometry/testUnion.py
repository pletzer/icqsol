#!/usr/bin/env python

"""
Test union operation
@author pletzer@psu.edu
"""

from icqsol.tools.geometry.icqSphere import Sphere

s1 = Sphere(radius=1.0, origin=(0., 0., 0.), n_theta=16, n_phi=8)
s2 = Sphere(radius=1.0, origin=(2., 0., 0.), n_theta=16, n_phi=8)
geom = s1 + s2
geom.show(filename='testUnion.png')
