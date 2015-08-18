#!/usr/bin/env python

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

    def hot(self, f):
        """
        Get hot color 
        @param f field value
        @return red, green, blue components in range 0 to 255
        """
        x = (f - self.fmin)/(self.fmax - self.fmin)
        r = int(255*math.cos((1. - x)*math.pi)**2 + 0.5)
        g = int(255*math.sin(x*math.pi)**2 + 0.5)
        b = int(255*math.cos(x*math.pi)**2 + 0.5)
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
