from icqsol.shapes.icqShapeManager import ShapeManager
from csg.core import CSG
from icqsol import util

"""
Test construction of shape from a list of polygons
@author alexander.net
"""
shape_mgr = ShapeManager(file_format=util.VTK_FORMAT, vtk_dataset_type=util.POLYDATA)
cube = CSG.cube()
# Should we test if shp == cube?
shp = shape_mgr.shapeFromPolygons(cube)
