#!/usr/bin/env python

"""
Test creation of simple ring object
@author alexander@gokliya.net
"""

from __future__ import print_function
from icqsol.shapes.icqShapeManager import ShapeManager
from icqsol import util

shape_mgr = ShapeManager(file_format=util.VTK_FORMAT, vtk_dataset_type=util.POLYDATA)
c1 = shape_mgr.createShape('cylinder', radius=0.4, origin=(0., 0., 0.), lengths=(0.5, 0., 0.))
c2 = shape_mgr.createShape('cylinder', radius=0.2, origin=(0., 0., -0.1), lengths=(1.0, 0., 0.))
geom = c1 - c2
shape_mgr.showShape(geom, filename='testRing.png')
