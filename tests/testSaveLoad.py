"""
Test save/load methods
@author alexander.net
"""

from icqsol.shapes.icqShapeManager import ShapeManager
from icqsol import util

shape_mgr = ShapeManager(file_format=util.VTK_FORMAT, vtk_dataset_type=util.POLYDATA)
shp = shape_mgr.createShape('box', origin=[0., 0., 0.], lengths=[1., 1., 1.])
shape_mgr.saveShape(shape=shp, file_name='testSaveLoad.vtk', file_type=util.ASCII)

shp2 = shape_mgr.loadAsShape('testSaveLoad.vtk')
shape_mgr.saveShape(shape=shp2, file_name='testSaveLoad2.vtk', file_type=util.ASCII)
