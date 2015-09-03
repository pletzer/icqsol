"""
Test save/load methods
@author pletzer@psu.edu
"""

from icqsol.shapes.icqShapeManager import ShapeManager

shape_mgr = ShapeManager()
shp = shape_mgr.createShape('box', origin=[0., 0., 0.], lengths=[1., 1., 1.])
shape_mgr.save(file_name='testSaveLoad.vtk', file_format='vtk', file_type='ascii', shape=shp)

shp2 = shape_mgr.load('testSaveLoad.vtk')
shape_mgr.save(file_name='testSaveLoad2.vtk', file_format='vtk', file_type='ascii', shape=shp2)
