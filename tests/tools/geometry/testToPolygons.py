from icqsol.tools.geometry.icqBox import Box

"""
Test conversion from a shape to a list of polygons
@author pletzer@psu.edu
"""

shp = Box(origin = [0., 0., 0.], lengths = [1., 1., 1.], n_x=1, n_y=1, n_z=1)

shp.debug()

# check whether one can convert to a list of polygons
polys = shp.toPolygons()

# check whether each polygon can be cloned
map(lambda p: p.clone(), polys)

