#!/usr/bin/env python

from __future__ import print_function
import operator
from functools import reduce

# https://pomax.github.io/bezierinfo/legendre-gauss.html
# u, weight
gaussPtsAndWeights = {
    1: [(0.5, 1.)],
    2: [(0.21132486540518711775, 0.50000000000000000000),
        (0.78867513459481288225, 0.50000000000000000000)],
    3: [(0.50000000000000000000, 0.44444444444444444444),
        (0.11270166537925831148, 0.27777777777777777778),
        (0.88729833462074168852, 0.27777777777777777778)],
    4: [(0.33000947820757186760, 0.32607257743127307131),
        (0.66999052179242813240, 0.32607257743127307131),
        (0.069431844202973712388, 0.17392742256872692869),
        (0.93056815579702628761, 0.17392742256872692869)],
    5: [(0.50000000000000000000, 0.28444444444444444444),
        (0.23076534494715845448, 0.23931433524968323402),
        (0.76923465505284154552, 0.23931433524968323402),
        (0.046910077030668003601, 0.11846344252809454376),
        (0.95308992296933199640, 0.11846344252809454376)],
    }


def lineQuadrature(order, pa, pb, func):
    """
    Compute the integral of a function over a line
    @param order integration order (1 <= order <= 5)
    @param pa starting point
    @param pb end point
    @param func function of position -> real
    @return integral
    """
    pba = pb - pa
    gpws = gaussPtsAndWeights[order]
    return pba * reduce(operator.add, [gpw[1]*func(pa + gpw[0]*pba) \
        for gpw in gpws])

##############################################################################


def testConstant():

    pa = 1.2
    pb = 2.0

    def f(x):
        return 1.0
    for order in range(1, 6):
        res = lineQuadrature(order, pa, pb, f)
        print('testConstant: order = ', order, ' integral = ', res)
        assert(abs(res - 0.8) < 1.e-12)


def testLinear():

    pa = 1.2
    pb = 2.0

    def f(x):
        return 2.0 - x
    for order in range(1, 6):
        res = lineQuadrature(order, pa, pb, f)
        print('testLinear: order = ', order, ' integral = ', res)
        assert(abs(res - 0.32) < 1.e-12)


def testQuadratic():

    pa = 1.2
    pb = 2.0

    def f(x):
        return 2.0 - x + 3*x**2
    for order in range(2, 6):
        res = lineQuadrature(order, pa, pb, f)
        print('testQuadratic: order = ', order, ' integral = ', res)
        assert(abs(res - 6.592) < 1.e-12)

if __name__ == '__main__':
    testConstant()
    testLinear()
    testQuadratic()
