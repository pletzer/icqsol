#!/usr/bin/env python

"""
Test union operation
@author pletzer@psu.edu
"""

from icqsol.shapes.shape_manager import ShapeManager

shape_mgr = ShapeManager()
s1 = shape_mgr.createShape('sphere', radius=1.0, origin=(0., 0., 0.),
                           n_theta=16, n_phi=8)
s2 = shape_mgr.createShape('sphere', radius=1.0, origin=(2., 0., 0.),
                           n_theta=16, n_phi=8)
geom = s1 + s2
shape_mgr.show(geom, filename='testUnion.png')
