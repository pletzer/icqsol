#!/usr/bin/env python
"""
Map an arbitrary triangle to a reference triangle
http://arxiv.org/pdf/1201.4938.pdf
"""

from __future__ import print_function
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

        # triangle edge lengths
        x0 = numpy.array(xa)
        self.xbdiff = numpy.array(xb) - x0
        self.xcdiff = numpy.array(xc) - x0

        xbdiffCrossXcdiff = numpy.cross(self.xbdiff, self.xcdiff)

        # the normal vector at the reference position
        self.normal = xbdiffCrossXcdiff / numpy.linalg.norm(xbdiffCrossXcdiff)

        # the edge length at theta = 0
        self.r1 = numpy.linalg.norm(self.xbdiff)

        # the edge length at theta = bigTheta
        self.r2 = numpy.linalg.norm(self.xcdiff)

        # the max integration angle
        self.bigTheta = math.atan2(
            numpy.linalg.norm(xbdiffCrossXcdiff), self.xbdiff.dot(self.xcdiff)
            )

        # the slope of the b-c triangle edge, x = r1 + a*y
        self.a = (self.r2*math.cos(self.bigTheta) - self.r1) / \
            self.r2/math.sin(self.bigTheta)

    def getPointBDiff(self):
        """
        Get the difference between the second and the first vertex points
        @return vector
        """
        return self.xbdiff

    def getPointCDiff(self):
        """
        Get the difference between the third and the first vertex points
        @return vector
        """
        return self.xcdiff

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
        Get the max theta angle of integration
        @return angle in radians
        """
        return self.bigTheta

##########################################################################


def test1():
    # simple, right-angle triangle
    rt = ReferenceTriangle([0., 0., 0.], [1., 0., 0.], [0., 2., 0.])
    assert(abs(rt.getR1() - 1.0) < 1.e-10)
    assert(abs(rt.getR2() - 2.0) < 1.e-10)
    assert(abs(rt.getA() + 1.0/2.0) < 1.e-10)
    assert(abs(rt.getBigTheta() - math.pi/2.0) < 1.e-10)


def test2():
    # triangle with offset
    rt = ReferenceTriangle([10., 20., 30.], [11., 20., 30.], [10., 22., 30.])
    assert(abs(rt.getR1() - 1.0) < 1.e-10)
    assert(abs(rt.getR2() - 2.0) < 1.e-10)
    assert(abs(rt.getA() + 1.0/2.0) < 1.e-10)
    assert(abs(rt.getBigTheta() - math.pi/2.0) < 1.e-10)


def test3():
    # big angle
    rt = ReferenceTriangle([10., 20., 30.], [11., 20., 30.], [9.0, 22., 30.])
    assert(abs(rt.getR1() - 1.0) < 1.e-10)
    assert(abs(rt.getR2() - math.sqrt(5.)) < 1.e-10)
    assert(abs(rt.getA() + 1.0) < 1.e-10)

if __name__ == '__main__':
    test1()
    test2()
    test3()
