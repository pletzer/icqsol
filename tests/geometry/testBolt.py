#!/usr/bin/env python

"""
Test creation of a bolt object
@author pletzer@psu.edu
"""

from icqsol.shapes.shape_manager import ShapeManager

shape_mgr = ShapeManager()
shaft = shape_mgr.createShape('cylinder', origin=[0., 0., 0.],
                              lengths=[1., 0., 0.], radius=0.1, n_theta=32)
head = shape_mgr.createShape('cone', origin=[-0.06, 0., 0.],
                             lengths=[0.14, 0., 0.], radius=0.25)
notch1 = shape_mgr.createShape('box', origin=[-0.06, -0.015, -0.15],
                               lengths=[0.03, 0.03, 0.30])
notch2 = shape_mgr.cloneShape(notch1)
shape_mgr.rotateShape(notch2, axis=(1., 0., 0.), angleDeg=90.0)

geom = head + shaft - notch1 - notch2

shape_mgr.save(geom, 'testBolt.vtk', file_format='vtk', file_type='ascii')
shape_mgr.show(geom, filename='testBolt.png')
