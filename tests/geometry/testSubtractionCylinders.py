#!/usr/bin/env python

"""
Test intersection operation
@author pletzer@psu.edu
"""

from icqsol.shapes.icqShapeManager import ShapeManager

shape_mgr = ShapeManager(file_format='vtk', vtk_dataset_type='POLYDATA')
s1 = shape_mgr.createShape('cylinder', radius=1.0, origin=(0., 0., 0.5), length=0.5)

# need to make the cylinder to subtract a little longer to avoid
# floating point issues
s2 = shape_mgr.createShape('cylinder', radius=0.5, origin=(0., 0., 0.), length=1.1)
geom = s1 - s2
shape_mgr.show(geom, filename='testSubtractionCylinders.png')
