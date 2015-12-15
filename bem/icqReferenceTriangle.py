#!/usr/bin/env python
"""
Map an arbitrary triangle to a reference triangle
http://arxiv.org/pdf/1201.4938.pdf
"""

import numpy
import math


class ReferenceTriangle:

    def __init__(self, xa, xb, xc):
        """
        Constructor
        @param xa first triangle vertex
        @param xb second triangle vertex
        @param xc third triangle vertex
        """
        self.xa = numpy.array(xa)
        self.xb = numpy.array(xb)
        self.xc = numpy.array(xc)
        xbPrime = self.xb - self.xa
        xcPrime = self.xc - self.xa
        xbPrimeCrossxaPrime = numpy.cross(xbPrime, xcPrime)
        self.normal = xbPrimeCrossxaPrime / numpy.linalg.norm(xbPrimeCrossxaPrime)
        self.r1 = numpy.linalg.norm(xbPrime)
        self.r2 = numpy.linalg.norm(xcPrime)
        self.bigTheta = math.atan2(
           numpy.linalg.norm(xbPrimeCrossxaPrime), 
           math.sqrt(xbPrime.dot(xcPrime)))
        self.a = (self.r2*math.cos(self.bigTheta) - self.r1) / \
                 self.r2/math.sin(self.bigTheta)

    def getR1(self):
        """
        Get the triangle length "r1" along theta = 0
        @return length
        """
        return self.r1

    def getR2(self):
        """
        Get the triangle length "r2" along theta = bigTheta
        @return length
        """
        return self.r2

    def getA(self):
        """
        Get the coefficient "a" in the opposite edge equation x' = r1 + a*y'
        @return coefficient
        """
        return self.a

    def getBigTheta(self):
        """
        Get the max theta angle
        @return angle in radians
        """
        return self.bigTheta

    def getElevation(self, point):
        """
        Get the elevation of "point" above (or below) the triangle
        """
        return (numpy.array(point) - self.xa).dot(self.normal)

##########################################################################

def test1():
    rt = ReferenceTriangle([0., 0., 0.], [1., 0., 0.], [0., 2., 0.])
    assert(abs(rt.getR1() - 1.0) < 1.e-10)
    assert(abs(rt.getR2() - 2.0) < 1.e-10)
    assert(abs(rt.getA() + 1.0/2.0) < 1.e-10)
    assert(abs(rt.getBigTheta() - math.pi/2.0) < 1.e-10)

if __name__ == '__main__':
    test1()
