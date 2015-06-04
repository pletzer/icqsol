#!/usr/bin/env python

"""
Test simple pipe object
@author pletzer@psu.edu
"""

from icqsol.tools.geometry.sqGeometry import Geometry
from icqsol.tools.geometry.sqCylinder import Cylinder

bxLo = (0., 0., 0.)
bxHi = (1., 1., 1.)
geom = Geometry(bxLo=bxLo, bxHi=bxHi, full=False)
geom += Cylinder(radius=0.4, origin=(0.5, 0.5, 0.5), length=0.5)
geom -= Cylinder(radius=0.2, origin=(0.5, 0.5, 0.5), length=0.5)
geom.show()
