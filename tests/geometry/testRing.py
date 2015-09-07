#!/usr/bin/env python

"""
Test creation of simple ring object
@author pletzer@psu.edu
"""

from icqsol.shapes.icqShapeManager import ShapeManager

shape_mgr = ShapeManager('vtk', 'POLYDATA')
c1 = shape_mgr.createShape('cylinder', radius=0.4, origin=(0., 0., 0.),
                           lengths=(0.5, 0., 0.))
c2 = shape_mgr.createShape('cylinder', radius=0.2, origin=(0., 0., -0.1),
                           lengths=(1.0, 0., 0.))
geom = c1 - c2
shape_mgr.show(geom, filename='testRing.png')
