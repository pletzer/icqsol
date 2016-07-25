#!/usr/bin/env python

from __future__ import print_function
import vtk
import numpy
from icqsol.bem.icqBaseSolver import BaseSolver
from icqsol.bem.icqPotentialIntegrals import PotentialIntegrals
from icqsol.bem.icqQuadrature import gaussPtsAndWeights
from icqsol.util.icqSharedLibraryUtils import getSharedLibraryName
from ctypes import cdll, POINTER, byref, c_void_p, c_double, c_long

FOUR_PI = 4. * numpy.pi

class BaseLaplaceSolver(BaseSolver):

    def __init__(self, pdata, max_edge_length, order=5):
        """
        Constructor
        @param pdata instance of vtkPolyData
        @param max_edge_length maximum edge length, used to turn
                               polygons into triangles
        """

        BaseSolver.__init__(self, pdata, max_edge_length, order)

        libName = getSharedLibraryName('icqLaplaceMatricesCpp')
        self.lib = cdll.LoadLibrary(libName)

        shp = (self.numTriangles, self.numTriangles)
        self.gMat = numpy.zeros(shp, numpy.float64)

        self.__computeResponseMatrix()

    def __computeDiagonalTerms(self):

        # Gauss points and weights
        gpws = gaussPtsAndWeights[self.order]
        npts = gpws.shape[1]
        xsis, etas, weights = gpws[0, :], gpws[1, :], gpws[2, :]

        # iterate over the source triangles
        for jSrc in range(self.numTriangles):

            ia, ib, ic = self.ptIdList[jSrc]

            # The triangle vertex positions
            paSrc = numpy.array(self.points.GetPoint(ia))
            pbSrc = numpy.array(self.points.GetPoint(ib))
            pcSrc = numpy.array(self.points.GetPoint(ic))
            dbSrc = pbSrc - paSrc
            dcSrc = pcSrc - paSrc

            # Iterate over the observer points
            g = 0
            for ipt in range(npts):
                # Observer point
                xObs = paSrc + xsis[ipt]*dbSrc + etas[ipt]*dcSrc
                # Three triangles having observer point as one corner
                pot0ab = PotentialIntegrals(xObs, paSrc, pbSrc, self.order)
                pot0bc = PotentialIntegrals(xObs, pbSrc, pcSrc, self.order)
                pot0ca = PotentialIntegrals(xObs, pcSrc, paSrc, self.order)
                g += weights[ipt] * (pot0ab.getIntegralOneOverR() + \
                                     pot0bc.getIntegralOneOverR() + \
                                     pot0ca.getIntegralOneOverR())

            self.gMat[jSrc, jSrc] = g / (-FOUR_PI)

    def __computeOffDiagonalTerms(self):

        addr = int(self.pdata.GetAddressAsString('vtkPolyData')[5:], 0)
        self.lib.computeOffDiagonalTerms(c_long(addr),
                                         self.gMat.ctypes.data_as(POINTER(c_double)))

    def __computeResponseMatrix(self):

        self.__computeDiagonalTerms()
        self.__computeOffDiagonalTerms()

    def getGreenMatrix(self):
        """
        Return the Green function matrix
        @return matrix
        """
        return self.gMat
