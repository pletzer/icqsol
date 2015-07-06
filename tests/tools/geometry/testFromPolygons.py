from icqsol.tools.geometry.icqShape import Shape
from csg.core import CSG

"""
Test construction of shape from a list of polygons
@author pletzer@psu.edu
"""

cube = CSG.cube()
shp = Shape()
shp.fromPolygons(cube.toPolygons())
