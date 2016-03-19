#!/usr/bin/env python

"""
Test union operation
@author alexander@gokliya.net
"""

from icqsol.shapes.icqShapeManager import ShapeManager
from icqsol import util

shape_mgr = ShapeManager(file_format=util.VTK_FORMAT, vtk_dataset_type=util.POLYDATA)
s1 = shape_mgr.createShape('cylinder', radius=1.0, origin=(0., 0., 0.5), length=0.5)
s2 = shape_mgr.createShape('cylinder', radius=0.5, origin=(0., 0., 0.), length=1.0)
geom = s1 + s2
shape_mgr.show(geom, filename='testUnionCylinders.png')
