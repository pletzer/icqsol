from __future__ import print_function
from ctypes import cdll, POINTER, byref, c_void_p, c_double, c_long
import numpy
import pkg_resources
from icqsol.shapes.icqShapeManager import ShapeManager
from icqsol import util

# Extract the shared library from the egg
libName = pkg_resources.resource_filename('icqsol', 'icqInsideLocatorCpp.so')

# Open the shared library 
lib = cdll.LoadLibrary(libName)

# Opaque handle
handle = c_void_p(0)

# load vtkPolyData from file
shape_mgr = ShapeManager(file_format=util.VTK_FORMAT, vtk_dataset_type='POLYDATA')
s = shape_mgr.createShape('box', origin=(-0.5, -0.5, -0.5),)
pdata = shape_mgr.shapeToVTKPolyData(s)

# Constructor
addr = int(pdata.GetAddressAsString('vtkPolyData')[5:], 0)
lib.icqInsideLocatorInit(byref(handle), c_long(addr))

# Tests

tests = [((0., 0., 0.), 1), # center
         ((10, 0., 0.), 0), # well outside
         ((0.1, 0.2, 0.3), 1), # somewhere inside
         ((0.499, 0., 0.), 1), # close to a face
         ((0.499, 0.499, 0.), 1), # close to an edge
         ((0.499, 0.499, 0.499), 1), # close to a vertex
         ((0.501, 0.499, 0.499), 0), # just outside
         ((0.5, 0.5, 0.5), 0), # on face and on vertex
         ((0.5, 0.2, 0.5), 0), # on face and on edge
         ((0.5, 0.2, 0.3), 1), # right on a face but well inside triangle
         ]
for test in tests:
	point = numpy.array(test[0])
	result = test[1]
	print('point = ', point)
	inside = lib.icqInsideLocatorIsPointInside(byref(handle), point.ctypes.data_as(POINTER(c_double)))
	print('inside = ', inside)
	assert(inside - result == 0)

# Destructor
lib.icqInsideLocatorDel(byref(handle))
