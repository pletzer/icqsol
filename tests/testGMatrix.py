#!/usr/bin/env python

import argparse
from icqsol.shapes.icqShapeManager import ShapeManager
from icqsol.bem.icqLaplaceMatrices import LaplaceMatrices
from icqsol import util

parser = argparse.ArgumentParser(description='Compute response matrices.')

parser.add_argument('--input', dest='input', default='',
                    help='Input file (VTK)')
args = parser.parse_args()

shape_mgr = ShapeManager(file_format='vtk', vtk_dataset_type='POLYDATA')
pdata = shape_mgr.loadAsVtkPolyData(args.input)
response = LaplaceMatrices(pdata,
                           max_edge_length=float('inf'),
                           order=8,
                           maxN=10)
print 'G = ', response.getGreenMatrix()