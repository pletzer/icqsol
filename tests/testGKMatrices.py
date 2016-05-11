#!/usr/bin/env python

from __future__ import print_function
from icqsol.shapes.icqShapeManager import ShapeManager
from icqsol.bem.icqLaplaceMatrices import LaplaceMatrices
from icqsol import util
import numpy

shape_mgr = ShapeManager(file_format='vtk', vtk_dataset_type='POLYDATA')
pdata = shape_mgr.loadAsVtkPolyData('data/2triangles.vtk')
response = LaplaceMatrices(pdata,
                           max_edge_length=float('inf'),
                           order=5)
g = response.getGreenMatrix()
k = response.getNormalDerivativeGreenMatrix()
print('g = ', end="\n")
print g
print('k = ', end="\n")
print k

gExact = numpy.array([[0.191561, 0.0378501], [0.0378501, 0.191561]])
kExact = numpy.array([[0, -0.0344229], [-0.0344229, 0]])
print('g Mathematica = ', end="\n")
print(gExact, end="\n")
print('k Mathematica = ', end="\n")
print(kExact, end="\n")

print('g error = ', end="\n")
print(g - gExact, end="\n")
print('k error = ', end="")
print(k - kExact, end="\n")

print('g^{-1} = ', numpy.linalg.inv(g), end="\n")

# residue
k[0, 0] = -0.5
k[1, 1] = -0.5
