#!/usr/bin/env python

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
print 'g = '
print g
print 'k = '
print k

print 'g Mathematica = '
print numpy.array([[0.191561, 0.0378501], [0.0378501, 0.191561]])
print 'k Mathematica = '
print numpy.array([[0, 0.0344229], [0.0344229, 0]])

print 'g^{-1} = ', numpy.linalg.inv(g)

# residue
k[0, 0] = -0.5
k[1, 1] = -0.5

normalDerivResponseMat = numpy.linalg.inv(g).dot(k)
print 'normal deriv response = ', normalDerivResponseMat

print 'normal deriv to [0, 1] = ', normalDerivResponseMat.dot([0, 1.])
