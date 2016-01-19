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
s = shape_mgr.createShape('sphere', origin=(0., 0., 0.), radius=1.0)
pdata = shape_mgr.shapeToVTKPolyData(s)

# Constructor
addr = int(pdata.GetAddressAsString('vtkPolyData')[5:], 0)
lib.icqInsideLocatorInit(byref(handle), c_long(addr))

# Tests

tests = [((0.,0.,0.), 1)]
for test in tests:
	point = numpy.array(test[0])
	result = test[1]
	print 'point = ', point
	inside = lib.icqInsideLocatorIsPointInside(byref(handle), point.ctypes.data_as(POINTER(c_double)))
	print 'inside = ', inside
	assert(inside - result == 0)

# Destructor
lib.icqInsideLocatorDel(byref(handle))