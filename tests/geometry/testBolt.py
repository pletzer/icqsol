#!/usr/bin/env python

"""
Test creation of a bolt object
@author pletzer@psu.edu
"""

from icqsol.shapes.icqShapeManager import ShapeManager
from icqsol import util

shape_mgr = ShapeManager(file_format=util.VTK_FORMAT, vtk_dataset_type=util.POLYDATA)
shaft = shape_mgr.createShape('cylinder', origin=[0., 0., 0.], lengths=[1., 0., 0.], radius=0.1, n_theta=32)
head = shape_mgr.createShape('cone', origin=[-0.06, 0., 0.], lengths=[0.14, 0., 0.], radius=0.25)
notch1 = shape_mgr.createShape('box', origin=[-0.06, -0.015, -0.15], lengths=[0.03, 0.03, 0.30])
notch2 = shape_mgr.cloneShape(notch1)
shape_mgr.rotateShape(notch2, axis=(1., 0., 0.), angleDeg=90.0)

geom = head + shaft - notch1 - notch2

shape_mgr.saveShape(shape=geom, file_name='testBolt.vtk', file_type=util.ASCII)
shape_mgr.show(geom, filename='testBolt.png')
