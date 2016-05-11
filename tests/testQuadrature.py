#!/usr/bin/env python

"""
Test Gauss quadrature
@author alexander@gokliya.net
"""

from __future__ import print_function
import unittest
import numpy
from math import exp, e
from icqsol.bem.icqQuadrature import triangleQuadrature

class TestQuadrature(unittest.TestCase):

    def testTrianglePolynomial(self):

        pa = numpy.array([0., 0., 0.])
        pb = numpy.array([1., 0., 0.])
        pc = numpy.array([0., 1., 0.])

        def func(point):
            x, y, z = point
            return 1. + 2*x + 3*y**2

        exact = 13./12.
        
        order = 1
        value = triangleQuadrature(order=order,
                                   pa=pa, pb=pb, pc=pc,
                                   func=func)
        assert(abs(value - exact) < 0.084)

        for order in range(2, 9):
            value = triangleQuadrature(order=order,
                                       pa=pa, pb=pb, pc=pc,
                                       func=func)
            assert(abs(value - exact) < 1.e-10)

    def testTrianglePolynomial2(self):

        pa = numpy.array([1., 2., 3.])
        pb = numpy.array([1.1, 2.3, 3.4])
        pc = numpy.array([0.8, 1., 0.])
        jacVec = numpy.cross(pb - pa, pc - pa)
        jac = numpy.sqrt(jacVec.dot(jacVec))

        def func(point):
            x, y, z = point
            return 1. + 2*x + 3*y**2 + 3*z**3

        exact = jac*780053./30000.
    
        order = 1
        value = triangleQuadrature(order=order,
                                   pa=pa, pb=pb, pc=pc,
                                   func=func)
        error = value - exact
        assert(abs(error) < 2.9)
                               
        order = 2
        value = triangleQuadrature(order=order,
                                   pa=pa, pb=pb, pc=pc,
                                   func=func)
        error = value - exact
        assert(abs(error) < 0.05)
        
        for order in range(3, 9):
            value = triangleQuadrature(order=order,
                                       pa=pa, pb=pb, pc=pc,
                                       func=func)
            error = value - exact
            assert(abs(value - exact) < 1.e-12)

    def testTriangleExp(self):

        pa = numpy.array([1., 2., 3.])
        pb = numpy.array([1.1, 2.3, 3.4])
        pc = numpy.array([0.8, 1., 0.])
        jacVec = numpy.cross(pb - pa, pc - pa)
        jac = numpy.sqrt(jacVec.dot(jacVec))

        def func(point):
            x, y, z = point
            return exp(x)

        exact = jac*50/3*(e**(4./5.)-3*e + 2*e**(11./10.))
            
        for order in range(1, 9):
            value = triangleQuadrature(order=order,
                                       pa=pa, pb=pb, pc=pc,
                                       func=func)
            error = value - exact
            print('order = {0} value = {1} error = {2}'.format(order, value, error))
            assert(abs(value - exact) < 0.002)


if __name__ == '__main__':
    unittest.main()
