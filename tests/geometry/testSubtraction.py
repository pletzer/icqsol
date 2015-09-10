#!/usr/bin/env python

"""
Test intersection operation
@author pletzer@psu.edu
"""

from icqsol.shapes.icqShapeManager import ShapeManager

shape_mgr = ShapeManager(file_format='vtk', vtk_dataset_type='POLYDATA')
s1 = shape_mgr.createShape('sphere', radius=1.0, origin=(0., 0., 0.))
s2 = shape_mgr.createShape('sphere', radius=1.0, origin=(1., 0., 0.))
geom = s1 - s2
shape_mgr.show(geom, filename='testSubtraction.png')
