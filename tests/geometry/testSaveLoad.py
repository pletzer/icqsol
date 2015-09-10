"""
Test save/load methods
@author pletzer@psu.edu
"""

from icqsol.shapes.icqShapeManager import ShapeManager

shape_mgr = ShapeManager(file_format='vtk', vtk_dataset_type='POLYDATA')
shp = shape_mgr.createShape('box', origin=[0., 0., 0.], lengths=[1., 1., 1.])
shape_mgr.saveShape(shape=shp, file_name='testSaveLoad.vtk', file_type='ascii')

shp2 = shape_mgr.loadAsShape('testSaveLoad.vtk')
shape_mgr.saveShape(shape=shp2, file_name='testSaveLoad2.vtk', file_type='ascii')
