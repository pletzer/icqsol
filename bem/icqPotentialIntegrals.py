#!/usr/bin/env python
"""
Analytic expressions for Laplacian singular kernel integrals
http://arxiv.org/pdf/1201.4938.pdf
"""

from __future__ import print_function
from math import cos, sin
from icqsol.bem.icqReferenceTriangle import ReferenceTriangle
from icqsol.bem.icqQuadrature1D import lineQuadrature


class PotentialIntegrals:

    def __init__(self, xa, xb, xc, order=5):
        """
        Constructor
        @param xa first triangle vertex
        @param xb second triangle vertex
        @param xc third triangle vertex
        @param order Gauss integration order
        """
        rt = ReferenceTriangle(xa, xb, xc)
        self.r1 = rt.getR1()
        self.r2 = rt.getR2()
        self.a = rt.getA()
        self.bigTheta = rt.getBigTheta()
        self.sinBigTheta = sin(self.bigTheta)
        self.tanBigTheta = self.sinBigTheta / cos(self.bigTheta)
        self.xbdiff = rt.getPointBDiff()
        self.xcdiff = rt.getPointCDiff()
        self.order = order

    def bigR(self, t):
        return self.r1/(cos(t) - self.a*sin(t))

    def getIntegralOneOverR(self):
        def integrand(t):
            return self.bigR(t)
        return lineQuadrature(self.order, 0.0, self.bigTheta, integrand)

##########################################################################


def testRightTriangle():
    for order in range(1, 6):
        potInt = PotentialIntegrals([0., 0., 0.],
                                    [1., 0., 0.], [0., 1., 0.], order)
        print('order = ', order,
            ' integral of 1/R is ', potInt.getIntegralOneOverR())

if __name__ == '__main__':
    testRightTriangle()
