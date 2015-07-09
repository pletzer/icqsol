#!/usr/bin/env python

"""
Test union operation
@author pletzer@psu.edu
"""

from icqsol.tools.geometry.icqSphere import Sphere

s1 = Sphere(radius=1.0, origin=(0., 0., 0.))
s2 = Sphere(radius=1.0, origin=(2., 0., 0.))
geom = s1 + s2
geom.show(filename='testUnion.png')
