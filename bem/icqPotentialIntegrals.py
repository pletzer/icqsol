#!/usr/bin/env python
"""
Analytic expressions for Laplacian singular kernel integrals
http://arxiv.org/pdf/1201.4938.pdf
"""

import numpy
from math import cos, sin, log, sqrt, pi, asin

class PotentialIntegrals:

    def __init__(self, z, r1, r2, bigTheta):
        """
        Constructor
        @param z distance of target normal to the triangle
        @param r1 distance from a -> b 
        @param r2 distance from a -> c
        @param bigTheta angle at reference point a
        """
        self.r1 = r1
        self.r2 = r2
        self.z = z
        zSquare = z**2
        self.bigTheta = bigTheta
        
        r1Square = r1**2
        r2Square = r2**2
        cosBigTheta = cos(bigTheta)
        sinBigTheta = sin(bigTheta)
        onePlusASquare = 1. + self.a**2

        self.a = (r2*cosbigTheta - r1)/(r2*sinBigTheta)
        betaSquare = (r1Square + zSquare*onePlusASquare)/onePlusASquare

        self.beta = sqrt(betaSquare)
        self.alphaSquare = zSquare/betaSquare
        self.alphaPrimeSquare = 1. - self.alphaSquare)
        self.alphaPrime = sqrt(self.alphaPrimeSquare)

        # m, n
        self.bigJ = {
        (0, 0): lambda t: t,
        (0, -1): lambda t: log((1. + sin(t))/(1. - sin(t))),
        (0, 1): lambda t: sin(t),
        (1, 0): lambda t: -cos(t),
        (0, 3): lambda t: sin(t) - sin(t)**3/3.,
        (1, -2): lambda t: 1./cos(t),
        (1, 2): lambda t: -cos(t)**3/3.,
        (2, -1): lambda t: -sin(t) + log(pi/4. + t/2.),
        }

        def bigDelta(t):
            return 1. - self.alphaSquare*sin(t)**2

        # p, m, n
        self.bigI = {
        (-3, 0, 1): lambda t: sin(t)/bigDelta(t),
        (-3, 1, 0): lambda t: -cos(t)/bigDelta(t)/self.alphaPrimeSquare,
        (-3, 0, 3): lambda t: -self.alphaPrimeSquare*sin(t)/self.alphaSquare/bigDelta(t) + asin(self.alpha*sin(t))/self.alpha**3,
        (-1, 0, 1): lambda t: asin(self.alpha*sin(t))/self.alpha,
        (-1, 1, 0): lambda t: -log(self.alpha*cos(t) + bigDelta(t))/self.alpha,
        (-1, 2, -1): lambda t: log((bigDelta(t) + self.alphaPrime*sin(t))/(bigDelta(t) - self.alphaPrime*sin(t)))/2./self.alphaPrime - asin(self.alpha*sin(t))/self.alpha,
        (1, 0, -1): lambda t: self.alphaPrime*log((bigDelta(t) + self.alphaPrime*sin(t))/(bigDelta(t) - self.alphaPrime*sin(t)))/2. - self.alpha*asin(self.alpha*sin(t)),
        (1, 2, -1): lambda t: -bigDelta(t)*sin(t)/2. + (2.*self.alphaSquare - 1.)*asin(self.alpha*sin(t))/2./self.alpha + self.alphaPrime*log((bigDelta(t) + self.alphaPrime*sin(t))/(bigDelta(t) - self.alphaPrime*sin(t)))/2.,
        (1, 1, -2): lambda t: bigDelta(t)/cos(t) - self.alpha*log(self.alpha*cos(t) + bigDelta(t)),
        (1, 1, 0): lambda t: -bigDelta(t)*cos(t)/2. - self.alphaPrimeSquare*log(self.alpha*cos(t) + bigDelta(t))/2./self.alpha,
        }

    def getOneOverR(self):
        if self.z != 0:
            return self.beta*self.bigI[1, 0, -1] - abs(self.z)*self.bigJ[0, 0]
        else:
            return self.beta*self.bigJ[0, -1]


    def getOneOverRCube(self):
        if self.z != 0:
            return -self.bigI[-1, 0, 1]/self.beta + self.bigJ[0, 0]/abs(self.z)
        else:
            return -self.bigJ[0, 1]/self.beta

if __name__ == '__main__':
    test()
