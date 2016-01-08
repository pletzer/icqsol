from ctypes import cdll, POINTER, byref, c_void_p, c_double
import numpy
import pkg_resources
from icqsol.bem.icqQuadrature import triangleQuadrature

# Extract the shared library from the egg
libName = pkg_resources.resource_filename('icqsol', 'icqQuadratureCpp.so')

# Open the shared library 
lib = cdll.LoadLibrary(libName)

# Opaque handle
handle = c_void_p(0)

# Constructor
lib.icqQuadratureInit(byref(handle))

# Tests

def func(p):
	pNorm = numpy.linalg.norm(p)
	return 1.0/pNorm

maxOrder = lib.icqQuadratureGetMaxOrder(byref(handle))

pa = numpy.array([0., 0., 0.])
pb = numpy.array([1., 0., 0.])
pc = numpy.array([0., 1., 0.])

paPtr = pa.ctypes.data_as(POINTER(c_double))
pbPtr = pb.ctypes.data_as(POINTER(c_double))
pcPtr = pc.ctypes.data_as(POINTER(c_double))

lib.icqQuadratureEvaluate.restype = c_double
for order in range(1, maxOrder + 1):
	integral = lib.icqQuadratureEvaluate(byref(handle), order, paPtr, pbPtr, pcPtr)
	integral2 = triangleQuadrature(order=order, pa=pa, pb=pb, pc=pc, func=func)
	print 'order = ', order, ' integral = ', integral, ' integral2 = ', integral2
	assert(abs(integral - integral2) < 1.e-10)



# Destructor
lib.icqQuadratureDel(byref(handle))