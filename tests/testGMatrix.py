#!/usr/bin/env python

from __future__ import print_function
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
k = response.getNormalDerivativeGreenMatrix()

# residue
n = k.shape[0]
for i in range(n):
    k[i, i] += 0.5

print('g:', end="\n")
print g

gM1 = numpy.linalg.inv(g)
print('g^{-1}: ', end="\n")
print(gM1, end="\n")

print('g^{-1} * g: ', end="\n")
print(gM1.dot(g), end="\n")

print('k + residue: ', end="\n")
print k

print('kres eigenvalues', end="\n")
print(numpy.linalg.eig(k), end="\n")

print('g^{-1} * kres: ', end="\n")
response = gM1.dot(k)
print(gM1.dot(k), end="\n")

print('g^{-1} * kres * 1: ', end="\n")
print(response.dot([1.,]*n), end="\n")




