#!/usr/bin/env python

"""
Solving Laplace equation
@author alexander@gokliya.net
"""

import argparse
from icqsol.shapes.icqShapeManager import ShapeManager
from icqsol.bem.icqLaplaceMatrices2 import LaplaceMatrices2
from icqsol import util

description = 'Apply a surface field to a shape'
parser = argparse.ArgumentParser(description=description)

parser.add_argument('--n_theta', dest='n_theta', default=4, type=int,
                    help='Number of longitude sections')

parser.add_argument('--n_phi', dest='n_phi', default=2, type=int,
                    help='Number of latitude sections')

parser.add_argument('--order', dest='order', default=5, type=int,
                    help='Order of the Gauss quadrature scheme')

args = parser.parse_args()
assert(args.n_theta >= 4)
assert(args.n_phi >= 2)
assert(args.order >= 1)
assert(args.order <= 5)

shape_mgr = ShapeManager(file_format=util.VTK_FORMAT,
                         vtk_dataset_type='POLYDATA')
s = shape_mgr.createShape('sphere',
                          origin=(0., 0., 0.),
                          radius=1.0,
                          n_theta=args.n_theta,
                          n_phi=args.n_phi)

pdata = shape_mgr.shapeToVTKPolyData(s)
response = LaplaceMatrices2(pdata,
                            max_edge_length=float('inf'),
                            order=args.order)

# compute the normal electric field from the potential
expr = '1.0/(4.*pi*sqrt(x**2 + y**2 + z**2))'
en = response.computeNeumannFromDirichlet(expr)
print 'normal electric field: ', en