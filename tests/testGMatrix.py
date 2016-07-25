#!/usr/bin/env python

from __future__ import print_function
from icqsol.shapes.icqShapeManager import ShapeManager
from icqsol.bem.icqLaplaceSolver import LaplaceSolver
from icqsol import util
import numpy

shape_mgr = ShapeManager(file_format='vtk', vtk_dataset_type='POLYDATA')
pdata = shape_mgr.loadAsVtkPolyData('data/tet.vtk')
slv = LaplaceSolver(pdata,
                    max_edge_length=float('inf'),
                    order=5)
g = slv.getGreenMatrix()

print('g:')
print g

gM1 = numpy.linalg.inv(g)
print('g^{-1}: ')
print(gM1)

print('g^{-1} * g: ')
print(gM1.dot(g))




