#!/usr/bin/env python

"""
Test union operation
@author alexander.net
"""

from icqsol.shapes.icqShapeManager import ShapeManager
from icqsol import util

shape_mgr = ShapeManager(file_format=util.VTK_FORMAT, vtk_dataset_type=util.POLYDATA)
s1 = shape_mgr.createShape('sphere', radius=1.0, origin=(0., 0., 0.), n_theta=16, n_phi=8)
s2 = shape_mgr.createShape('sphere', radius=1.0, origin=(2., 0., 0.), n_theta=16, n_phi=8)
geom = s1 + s2
shape_mgr.showShape(geom, filename='testUnion.png')
