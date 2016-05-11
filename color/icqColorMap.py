#!/usr/bin/env python

from __future__ import print_function
import math


class ColorMap:


    def __init__(self, fmin, fmax):
        """
        Constructor 
        @param fmin minimum field value
        @param fmax maximum field value
        """

        self.fmin = fmin
        self.fmax = fmax

    @classmethod
    def hatFunc(self, x, x0, x1, x2, x3):
        """
        Hat function
        @param x normalized abscissa, between 0 and 1
        @param x0 x value where function rises from zero
        @param x1 x value where function reaches 1
        @param x2 x value where function drops from 1
        @param x3 x value where function reaches 0
        @return value
        """
        slope01 = 1./(x1 - x0)
        slope23 = 1./(x3 - x2)
        return min(1., max(0., slope01*(x - x0))) \
             - min(1., max(0., slope23*(x - x2)))

    def hot(self, f):
        """
        Get hot color 
        @param f field value
        @return red, green, blue components in range 0 to 255
        """
        x = (f - self.fmin)/(self.fmax - self.fmin)
        r = int(255*math.cos((1. - x)*math.pi/2.0)**2 + 0.5)
        g = int(255*math.sin(x*math.pi)**2 + 0.5)
        b = int(255*math.cos(x*math.pi/2.0)**2 + 0.5)
        return r, g, b

    def cold(self, f):
        """
        Get cold color 
        @param f field value
        @return red, green, blue components in range 0 to 255
        """
        r, g, b = self.hot(f)
        # reverse order
        return b, g, r

    def gnu(self, f):
        """
        Get gnu color
        @param f field value
        @return red, green, blue components in range 0 to 255
        """
        x = (f - self.fmin)/(self.fmax - self.fmin)
        r = int(255*self.hatFunc(x, 3./8., 5./8., 7./8., 9./8.) + 0.5)
        g = int(255*self.hatFunc(x, 1./8., 3./8., 5./8., 7./8.) + 0.5)
        b = int(255*self.hatFunc(x, -1./8., 1./8., 3./8., 5./8.) + 0.5)
        return r, g, b

    def blackbody(self, f):
        """
        Get black body color 
        @param f field value
        @return red, green, blue components in range 0 to 255
        """
        x = (f - self.fmin)/(self.fmax - self.fmin)
        r = int(255*min(1., 2*math.cos(x*math.pi/2.)**2) + 0.5)
        g = int(255*math.sin(x*math.pi)**2 + 0.5)
        b = int(255*min(1., 2*math.cos((1. - x)*math.pi/2.)**2) + 0.5)
        return r, g, b

##############################################################################


def testHat():
    n = 10
    dx = 1.0/float(10)
    x0, x1, x2, x3 = 0.1, 0.3, 0.4, 0.8
    for i in range(n + 1):
        x = i * dx
        print('x = {} y = {}'.format(x, ColorMap.hatFunc(x, x0, x1, x2, x3)))

def testGnu():
    fmin, fmax = 0., 1.
    cm = ColorMap(fmin, fmax)
    n = 10
    df = (fmax - fmin)/float(10)
    for i in range(n + 1):
        f = fmin + i * df
        print('f = {} rgb = {}'.format(f, cm.gnu(f)))

if __name__ == '__main__': 
    testHat()
    testGnu()
