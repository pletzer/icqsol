#!/usr/bin/env python

"""
    Test colormaps
    @author alexander@gokliya.net
    """

from icqsol.shapes.icqShapeManager import ShapeManager
from icqsol import util

shape_mgr = ShapeManager()
s = shape_mgr.createShape('box', origin=(0., 0., 0.), lengths=[10., 1., 1.])
pdata = s.getVtkPolydata()
shape_mgr.setWriter(file_format='vtk', vtk_dataset_type='POLYDATA')

# refine
sr = shape_mgr.refineVtkPolyData(pdata, min_edge_length=0.05)
shape_mgr.saveVtkPolyData(sr.getVtkPolyData(), file_name='boxRefined.vtk', file_type='ascii')

# coarsen
sc = shape_mgr.coarsenVtkPolyData(sr.getVtkPolyData(), max_cell_area=0.2)
shape_mgr.saveVtkPolyData(sc.getVtkPolyData(), file_name='boxRefinedCoarsened.vtk', file_type='ascii')
