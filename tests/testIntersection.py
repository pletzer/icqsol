#!/usr/bin/env python

"""
Test intersection operation
@author alexander.net
"""

from icqsol.shapes.icqShapeManager import ShapeManager
from icqsol import util

shape_mgr = ShapeManager(file_format=util.VTK_FORMAT, vtk_dataset_type=util.POLYDATA)
s1 = shape_mgr.createShape('sphere', radius=1.0, origin=(0., 0., 0.))
s2 = shape_mgr.createShape('sphere', radius=1.0, origin=(1., 0., 0.))
geom = s1 * s2
shape_mgr.showShape(geom, filename='testIntersection.png')
