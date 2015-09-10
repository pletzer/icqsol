from icqsol.shapes.icqShapeManager import ShapeManager
from csg.core import CSG

"""
Test construction of shape from a list of polygons
@author pletzer@psu.edu
"""
shape_mgr = ShapeManager(file_format='vtk', vtk_dataset_type='POLYDATA')
cube = CSG.cube()
# Should we test if shp == cube?
shp = shape_mgr.shapeFromPolygons(cube)
