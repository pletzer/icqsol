#!/usr/bin/env python

"""
Solve Laplace equation given Dirichlet boundary conditions
"""
import argparse
import time
import sys
import re
from numpy import linspace

from icqsol.shapes.icqShapeManager import ShapeManager
from icqsol.bem.icqLaplaceMatrices import LaplaceMatrices
from icqsol import util

# time stamp
tid = re.sub(r'\.', '', str(time.time()))

description = 'Solve Laplace equation given Dirichlet boundary conditions'
parser = argparse.ArgumentParser(description=description)

parser.add_argument('--input', dest='input', default='',
                    help='Input file (PLY or VTK).')

parser.add_argument('--dirichlet', dest='dirichlet', default='sin(pi*x)*cos(pi*y)*z',
                    help='Dirichlet boundary conditions, expression of x, y, and z.')

parser.add_argument('--refine', dest='refine', default=0.0, type=float,
                    help='Maximum edge length (use 0 if no refinement).')
                    
parser.add_argument('--order', dest='order', default=5, type=int,
                    help='Quadrature order (1 <= order <= 5).')

parser.add_argument('--input_name', dest='input_name', default='voltage',
                    help='Set the name of the input field.')

parser.add_argument('--output_name', dest='output_name', default='normal_electric_field',
                    help='Set the name of the output field.')

parser.add_argument('--ascii', dest='ascii', action='store_true',
                    help='Save data in ASCII format (default is binary).')

parser.add_argument('--output', dest='output',
                    default='solveLaplace-{0}.vtk'.format(tid),
                    help='VTK Output file.')
                    
parser.add_argument('--verbose', dest='verbose', action='store_true',
                    help='VTK Output file.')

args = parser.parse_args()

assert(args.order >= 1)
assert(args.order <= 5)

if not args.dirichlet:
    print 'ERROR: must specify --dirichlet <expression>'
    sys.exit(2)

if not args.input:
    print 'ERROR: must specify input file: --input <file>'
    sys.exit(3)

# make sure the field names contain no spaces
args.input_name = re.sub('\s', '_', args.input_name)
args.output_name = re.sub('\s', '_', args.output_name)

# Get the format of the input - either vtk or ply.
file_format = util.getFileFormat(args.input)

if file_format == util.PLY_FORMAT:
    shape_mgr = ShapeManager(file_format=util.PLY_FORMAT)
else:
    # We have a VTK file, so Get the dataset type.
    vtk_dataset_type = util.getVtkDatasetType(args.input)
    shape_mgr = ShapeManager(file_format=util.VTK_FORMAT,
                             vtk_dataset_type=vtk_dataset_type)

pdata = shape_mgr.loadAsVtkPolyData(args.input)

maxEdgeLength = float('inf')
if args.refine > 0:
    maxEdgeLength = args.refine

solver = LaplaceMatrices(pdata, maxEdgeLength, order=args.order)

# set the output field names
solver.setPotentialName(args.input_name)
solver.setNormalDerivativeName(args.output_name)

# in place operation, pdata will be modified
normalDeriv = solver.computeNeumannFromDirichlet(args.dirichlet)

if args.verbose:
    minNormDeriv = min(normalDeriv)
    maxNormDeriv = max(normalDeriv)
    avgNormDeriv = normalDeriv.sum()/len(normalDeriv)
    print 'normal derivative min/avg/max: {0}/{1}/{2}'.format(minNormDeriv, 
                                                              avgNormDeriv,
                                                              maxNormDeriv)
    # compute the response matrices
    gMat = solver.getGreenMatrix()
    kMat = solver.getNormalDerivativeGreenMatrix()

    # compute the potential
    import numpy
    from math import pi, sin, cos, log, exp, sqrt
    potential = numpy.zeros(kMat.shape[0], numpy.float64)
    pointArray = solver.getPoints()
    cellArray = solver.getCells()
    numCells = cellArray.shape[0]
    for i in range(numCells):
        ia, ib, ic = cellArray[i, :]
        x, y, z = (pointArray[ia, :] + pointArray[ib, :] + pointArray[ic, :]) / 3.
        potential[i] = eval(args.dirichlet)
    error = gMat.dot(normalDeriv) - kMat.dot(potential) - 0.5*potential
    print 'total error: ', error.sum()

if args.output:
    # Always produce VTK POLYDATA.
    shape_mgr.setWriter(file_format=util.VTK_FORMAT,
                        vtk_dataset_type=util.POLYDATA)
    if args.ascii:
        file_type = util.ASCII
    else:
        file_type = util.BINARY
    shape_mgr.saveVtkPolyData(vtk_poly_data=solver.getVtkPolyData(),
                              file_name=args.output,
                              file_type=file_type)
