#!/usr/bin/env python

"""
Test simple ring object
@author pletzer@psu.edu
"""

from icqsol.tools.geometry.icqGeometry import Geometry
from icqsol.tools.geometry.icqCylinder import Cylinder

geom = Geometry()
geom += Cylinder(radius=0.4, origin=(0.5, 0.5, 0.5), length=0.5)
geom -= Cylinder(radius=0.2, origin=(0.5, 0.5, 0.5), length=0.5)
geom.show()
