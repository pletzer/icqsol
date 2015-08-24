from icqShapeManager import ShapeManager
from csg.core import CSG

"""
Test construction of shape from a list of polygons
@author pletzer@psu.edu
"""
shape_mgr = ShapeManager()
cube = CSG.cube()
# Should we test if shp == cube?
shp = shape_mgr.shapeFromPolygons(shape_mgr.shapeToPolygons(cube))
