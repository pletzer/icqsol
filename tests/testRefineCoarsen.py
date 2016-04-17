#!/usr/bin/env python

"""
    Test colormaps
    @author alexander@gokliya.net
    """

from icqsol.shapes.icqShapeManager import ShapeManager
from icqsol import util

shape_mgr = ShapeManager()
s = shape_mgr.createShape('box', origin=(0., 0., 0.), lengths=[1., 1.1, 1.2])
pdata = shape_mgr.shapeToVTKPolyData(s)
shape_mgr.setWriter(file_format='vtk', vtk_dataset_type='POLYDATA')

# refine
pr = shape_mgr.refineVtkPolyData(pdata, max_edge_length=0.10)
shape_mgr.saveVtkPolyData(pr, file_name='boxRefined.vtk', file_type='ascii')

# coarsen
pc = shape_mgr.coarsenVtkPolyData(pr, min_cell_area=0.2)
shape_mgr.saveVtkPolyData(pc, file_name='boxRefinedCoarsened.vtk', file_type='ascii')
