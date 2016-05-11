#!/usr/bin/env python

from __future__ import print_function
import numpy

class FiniteElementIntegrals:

	def __init__(self, vertices):
		"""
		Constructor
		@param vertices list of vertices
		"""

		# distance along each edge of the tettrahedron's vertices
		rEdge = [numpy.array(vertices[1]) - numpy.array(vertices[0]), 
		         numpy.array(vertices[2]) - numpy.array(vertices[0]), 
		         numpy.array(vertices[3]) - numpy.array(vertices[0])]

		# the parallelepiped volume of the cell
		self.jac = numpy.dot(numpy.cross(rEdge[0], rEdge[1]), rEdge[2])

		# an approximation of the co-variant metric tensor
		dxMat = numpy.zeros( (3, 3), numpy.float64 )
		for i in range(3):
			for j in range(3):
				dxMat[i, j] = numpy.dot(rEdge[i], rEdge[j])

		# contravariant approximation of the metric tensor
		self.invDxMat = numpy.linalg.inv(dxMat)

		self.dBasisOverDXi = numpy.array( [[-1., -1., -1.], 
			                               [1., 0., 0.], 
			                               [0., 1., 0.], 
			                               [0., 0., 1.]])

		# derivative of basis function with respect to cartesian coordinates
		self.grads = numpy.dot(self.dBasisOverDXi, self.invDxMat)

	def integrateBasisBasis(self):
		"""
		Compute the basis function times basis function cell integral
		@return a 4x4, contribution matrix containing the coupling between nodes
		"""
		res = numpy.zeros( (4, 4), numpy.float64)
		for i in range(4):
			res[i, i] = self.jac / 60.0
			res[i, i+1:] = self.jac / 120.0
			res[i+1:, i] = self.jac / 120.0
		return res

	def integrateGradientDotGradient(self):
		"""
		Compute the cell integral of the gradient of the basis function 
		dotted with the gradient of another gradient basis function
		@return a 4x4, contribution matrix containing the coupling between nodes
		"""
		integral = 1. / 6.
		res = numpy.zeros( (4, 4), numpy.float64 )
		# iterate over basis functions
		for iBasis1 in range(4):
			for iBasis2 in range(iBasis1, 4):
				dotProd = numpy.dot(self.grads[iBasis1, :], 
					                self.grads[iBasis2, :])
				res[iBasis1, iBasis2] = self.jac * integral * dotProd
				res[iBasis2, iBasis1] = res[iBasis1, iBasis2]
		return res

##############################################################################
def test():
	v0 = (0., 0., 0.)
	v1 = (1., 0., 0.)
	v2 = (0., 1., 0.)
	v3 = (0., 0., 1.)
	verts = (v0, v1, v2, v3)
	fei = FiniteElementIntegrals(verts)
	print('basis basis = ', fei.integrateBasisBasis())
	print('grad dot grad = ', fei.integrateGradientDotGradient())

if __name__ == '__main__':
	test()
