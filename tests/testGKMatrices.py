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

gExact = numpy.array([[0.191561, 0.0378501], [0.0378501, 0.191561]])
kExact = numpy.array([[0, -0.0344229], [-0.0344229, 0]])
print 'g Mathematica = '
print gExact
print 'k Mathematica = '
print kExact

print 'g error = '
print g - gExact
print 'k error = ',
print k - kExact

print 'g^{-1} = ', numpy.linalg.inv(g)

# residue
k[0, 0] = -0.5
k[1, 1] = -0.5
