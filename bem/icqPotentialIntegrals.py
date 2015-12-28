#!/usr/bin/env python
"""
Analytic expressions for Laplacian singular kernel integrals
http://arxiv.org/pdf/1201.4938.pdf
"""

import numpy
from math import cos, sin, sqrt
from icqsol.bem.icqReferenceTriangle import ReferenceTriangle
from icqsol.bem.icqQuadrature1D import lineQuadrature    

class PotentialIntegrals:

    def __init__(self, xa, xb, xc, order=5):
        """
        Constructor
        @param xa first triangle vertex
        @param xb second triangle vertex
        @param xc third triangle vertex
        """
        rt = ReferenceTriangle(xa, xb, xc)
        self.r1 = rt.getR1()
        self.r2 = rt.getR2()
        self.a = rt.getA()
        self.bigTheta = rt.getBigTheta()
        self.order = order
        
    def bigR(self, t):
        return self.r1/(cos(t) - self.a*sin(t))

    def getIntegralOneOverR(self, elev):
    	def integrand(t):
    	    return sqrt(self.bigR(t)**2 + elev**2)
    	return lineQuadrature(self.order, 0.0, self.bigTheta, integrand) \
    	       - self.bigTheta * abs(elev)

##########################################################################

def testRightTriangle():
    for order in range(1, 6):
        potInt = PotentialIntegrals([0., 0., 0.], 
            [1., 0., 0.], [0., 1., 0.], order)
        print 'order = ', order
        for elev in (-1., -0.1, -0.01, 0.0, 0.01, 0.1, 1.):
        	print '\telev = ', elev, \
            ' integral of 1/R is ', potInt.getIntegralOneOverR(elev)

if __name__ == '__main__':
    testRightTriangle()
