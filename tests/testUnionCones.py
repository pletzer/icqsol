#!/usr/bin/env python

"""
Test union operation
@author alexander@gokliya.net
"""

from __future__ import print_function
from icqsol.shapes.icqShapeManager import ShapeManager
from icqsol import util

shape_mgr = ShapeManager(file_format=util.VTK_FORMAT, vtk_dataset_type=util.POLYDATA)
s1 = shape_mgr.createShape('cone', radius=1.0, origin=(0., 0., 0.5),
                           length=0.5, n_rho=100, n_theta=32, n_z=100)
s2 = shape_mgr.createShape('cone', radius=0.5, origin=(0., 0., 0.),
                           length=1.0, n_rho=100, n_theta=32, n_z=100)
geom = s1 + s2
shape_mgr.show(geom, filename='testUnionCones.png')
