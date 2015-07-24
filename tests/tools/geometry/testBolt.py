#!/usr/bin/env python

"""
Test creation of a bolt object
@author pletzer@psu.edu
"""

from icqsol.tools.geometry.icqCone import Cone
from icqsol.tools.geometry.icqCylinder import Cylinder
from icqsol.tools.geometry.icqBox import Box

shaft = Cylinder(origin=[0., 0., 0.], lengths=[1., 0., 0.],
                 radius=0.1, n_theta=32)
head = Cone(origin=[-0.06, 0., 0.],
            lengths=[0.14, 0., 0.], radius=0.25)
notch1 = Box(origin=[-0.12, -0.015, -0.25],
             lengths=[0.03, 0.03, 0.25])
notch2 = notch1.rotate(axis=(1., 0., 0.), angleDeg=90.0)

geom = head + shaft - notch1 - notch2

geom.save('testBolt.vtk', file_format='vtk', file_type='ascii')
geom.show(filename='testBolt.png')
