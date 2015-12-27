#!/usr/bin/env python

from icqsol.shapes.icqShapeManager import ShapeManager
from icqsol.bem.icqLaplaceMatrices import LaplaceMatrices
from icqsol import util
import numpy

shape_mgr = ShapeManager(file_format='vtk', vtk_dataset_type='POLYDATA')
pdata = shape_mgr.loadAsVtkPolyData('data/tet.vtk')
response = LaplaceMatrices(pdata,
                           max_edge_length=float('inf'),
                           order=5)
g = response.getGreenMatrix()
gExact = numpy.array([[0.1915612707, 0.07688979401, 0.07688979401, 0.1312564214 + 1j*0.0416666667],
          [0.07688979401, 0.1915612707, 0.07688979401, 0.1312564214 + 1j*0.0416666667],
          [0.07688979401, 0.07688979401, 0.1915612707, 0.1312564214 + 1j*0.0416666667],
          [float('nan'), float('nan'), float('nan'), float('nan')]])

print g