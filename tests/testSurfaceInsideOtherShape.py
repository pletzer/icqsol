#!/usr/bin/env python
"""
Test intersection operation
@author pletzer@psu.edu
"""
from icqsol.shapes.icqShapeManager import ShapeManager
from icqsol import util

shape_mgr = ShapeManager(file_format=util.VTK_FORMAT, vtk_dataset_type=util.POLYDATA)
# box-box containment
a = shape_mgr.createShape('box', origin=(0., 0., 0.), lengths=(1., 1., 0.2))
b = shape_mgr.createShape('box', origin=(0.5, 0.5, 0.), lengths=(1., 1., 0.2))
s = shape_mgr.getBoundarySurfaceInsideShape(a, b)

shape_mgr.saveShape(shape=a, file_name='testInsideOtherShape_BoxBox_a.vtk', file_type=util.ASCII)
shape_mgr.saveShape(shape=b, file_name='testInsideOtherShape_BoxBox_b.vtk', file_type=util.ASCII)
shape_mgr.saveShape(shape=s, file_name='testInsideOtherShape_BoxBox_result.vtk', file_type=util.ASCII)

# box - sphere
a = shape_mgr.createShape('box', origin=(0., 0., 0.), lengths=(1., 1., 0.2))
b = shape_mgr.createShape('sphere', origin=(0.5, 0.5, 0.), radius=0.5)
s = shape_mgr.getBoundarySurfaceInsideShape(a, b)

shape_mgr.saveShape(shape=a, file_name='testInsideOtherShape_BoxSphere_a.vtk', file_type=util.ASCII)
shape_mgr.saveShape(shape=b, file_name='testInsideOtherShape_BoxSphere_b.vtk', file_type=util.ASCII)
shape_mgr.saveShape(shape=s, file_name='testInsideOtherShape_BoxSphere_result.vtk', file_type=util.ASCII)
