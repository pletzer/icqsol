from icqsol.tools.geometry.icqBox import Box
from icqsol.tools.geometry.icqSphere import Sphere

from csg.core import CSG

"""
Test conversion from a shape to a list of polygons
@author pletzer@psu.edu
"""

shp = Box(origin = [0., 0., 0.], lengths = [1., 1., 1.],)

#shp.debug()

# check whether one can convert to a list of polygons
polys = shp.toPolygons()

# check whether each polygon can be cloned
map(lambda p: p.clone(), polys)

# check that we can load the polygons
a = CSG.fromPolygons(polys)

shp2 = Sphere(radius=1.0, origin=(0., 0., 0.), n_theta=5, n_phi=2)
shp2.debug()
polys2 = shp2.toPolygons()
a2 = CSG.fromPolygons(polys2)
