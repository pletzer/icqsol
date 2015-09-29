#!/usr/bin/env python

"""
Test colormaps
@author pletzer@psu.edu
"""

from icqsol.shapes.icqShapeManager import ShapeManager
from icqsol import util

shape_mgr = ShapeManager()
s = shape_mgr.createShape('box', origin=(0., 0., 0.), lengths=[10., 1., 1.])
s2 = shape_mgr.refineShape(s, refine=4)

# add a field
pdata = shape_mgr.addSurfaceFieldFromExpressionToShape(s2, 'myField', 'x', [0.0])

# color
pdataHot = shape_mgr.colorSurfaceField(pdata, 'hot', field_name='myField')
pdataCold = shape_mgr.colorSurfaceField(pdata, 'cold', field_name='myField')
pdataGnu = shape_mgr.colorSurfaceField(pdata, 'gnu', field_name='myField')
pdataBlackbody = shape_mgr.colorSurfaceField(pdata, 'blackbody', field_name='myField')

# write files
shape_mgr.setWriter(file_format='vtk', vtk_dataset_type='POLYDATA')
shape_mgr.saveVtkPolyData(pdataHot, file_name='hot.vtk', file_type='ascii')
shape_mgr.saveVtkPolyData(pdataCold, file_name='cold.vtk', file_type='ascii')
shape_mgr.saveVtkPolyData(pdataGnu, file_name='gnu.vtk', file_type='ascii')
shape_mgr.saveVtkPolyData(pdataBlackbody, file_name='blackbody.vtk', file_type='ascii')
