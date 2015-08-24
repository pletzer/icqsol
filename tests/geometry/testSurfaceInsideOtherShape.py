#!/usr/bin/env python
"""
Test intersection operation
@author pletzer@psu.edu
"""
from icqsol.shapes.icqShapeManager import ShapeManager

shape_mgr = ShapeManager()
# box-box containment
a = shape_mgr.createShape('box', origin=(0., 0., 0.), lengths=(1., 1., 0.2))
b = shape_mgr.createShape('box', origin=(0.5, 0.5, 0.), lengths=(1., 1., 0.2))
s = shape_mgr.getBoundarySurfaceInsideShape(a, b)

shape_mgr.save(a, 'testInsideOtherShape_BoxBox_a.vtk', file_format='vtk',
               file_type='ascii')
shape_mgr.save(b, 'testInsideOtherShape_BoxBox_b.vtk', file_format='vtk',
               file_type='ascii')
shape_mgr.save(s, 'testInsideOtherShape_BoxBox_result.vtk',
               file_format='vtk', file_type='ascii')

# box - sphere
a = shape_mgr.createShape('box', origin=(0., 0., 0.), lengths=(1., 1., 0.2))
b = shape_mgr.createShape('sphere', origin=(0.5, 0.5, 0.), radius=0.5)
s = shape_mgr.getBoundarySurfaceInsideShape(a, b)

shape_mgr.save(a, 'testInsideOtherShape_BoxSphere_a.vtk', file_format='vtk',
               file_type='ascii')
shape_mgr.save(b, 'testInsideOtherShape_BoxSphere_b.vtk', file_format='vtk',
               file_type='ascii')
shape_mgr.save(s, 'testInsideOtherShape_BoxSphere_result.vtk',
               file_format='vtk', file_type='ascii')
