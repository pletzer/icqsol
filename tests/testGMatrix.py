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
k = response.getNormalDerivativeGreenMatrix()

# residue
n = k.shape[0]
for i in range(n):
    k[i, i] += 0.5

print 'g:'
print g

gM1 = numpy.linalg.inv(g)
print 'g^{-1}: '
print gM1

print 'g^{-1} * g: '
print gM1.dot(g)

print 'k + residue: '
print k

print 'kres eigenvalues'
print numpy.linalg.eig(k)

print 'g^{-1} * kres: '
response = gM1.dot(k)
print gM1.dot(k)

print 'g^{-1} * kres * 1: '
print response.dot([1.,]*n)




