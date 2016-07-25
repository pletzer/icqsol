#!/usr/bin/env python

from __future__ import print_function
from icqsol.shapes.icqShapeManager import ShapeManager
from icqsol.bem.icqLaplaceSolver import LaplaceSolver
from icqsol import util
import numpy

shape_mgr = ShapeManager(file_format='vtk', vtk_dataset_type='POLYDATA')
pdata = shape_mgr.loadAsVtkPolyData('data/2triangles.vtk')
slv = LaplaceSolver(pdata,
                    max_edge_length=float('inf'),
                    order=5)
g = slv.getGreenMatrix()
print('g = ')
print g

gExact = numpy.array([[0.191561, 0.0378501], [0.0378501, 0.191561]])
print('g Mathematica = ')
print(gExact)

print('g error = ')
print(g - gExact)

print('g^{-1} = ', numpy.linalg.inv(g))

