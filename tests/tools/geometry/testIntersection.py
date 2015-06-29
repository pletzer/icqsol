#!/usr/bin/env python

"""
Test intersection operation
@author pletzer@psu.edu
"""

from icqsol.tools.geometry.icqSphere import Sphere

s1 = Sphere(radius=1.0, origin=(0., 0., 0.))
s2 = Sphere(radius=1.0, origin=(1., 0., 0.))
geom = s1 * s2
geom.show() #filename='testIntersection.png')
