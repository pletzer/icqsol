#!/usr/bin/env python

"""
Test union operation
@author pletzer@psu.edu
"""

from icqsol.shapes.shape_manager import ShapeManager

shape_mgr = ShapeManager()
s1 = shape_mgr.createShape('cylinder', radius=1.0, origin=(0., 0., 0.5),
                           length=0.5)
s2 = shape_mgr.createShape('cylinder', radius=0.5, origin=(0., 0., 0.),
                           length=1.0)
geom = s1 + s2
shape_mgr.show(geom, filename='testUnionCylinders.png')
