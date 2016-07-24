#!/usr/bin/env python

from __future__ import print_function
import vtk
import numpy
from icqsol.shapes.icqRefineSurface import RefineSurface
from icqsol.bem.icqPotentialIntegrals import PotentialIntegrals
from icqsol.bem.icqQuadrature import gaussPtsAndWeights
from icqsol.util.icqSharedLibraryUtils import getSharedLibraryName
from ctypes import cdll, POINTER, byref, c_void_p, c_double, c_long

FOUR_PI = 4. * numpy.pi

class BaseSolver:

    def __init__(self, pdata, max_edge_length, order=5):
        """
        Constructor
        @param pdata instance of vtkPolyData
        @param max_edge_length maximum edge length, used to turn
                               polygons into triangles
        """
        libName = getSharedLibraryName('icqLaplaceMatricesCpp')
        self.lib = cdll.LoadLibrary(libName)

        # triangulate
        rs = RefineSurface(pdata)
        rs.refine(max_edge_length=max_edge_length)
        self.pdata = rs.getVtkPolyData()

        # store the point indices for each cell
        self.ptIdList = []
        ptIds = vtk.vtkIdList()
        polys = self.pdata.GetPolys()
        polys.InitTraversal()
        for i in range(polys.GetNumberOfCells()):
            polys.GetNextCell(ptIds)
            assert(ptIds.GetNumberOfIds() == 3)
            self.ptIdList.append([ptIds.GetId(0),
                                  ptIds.GetId(1),
                                  ptIds.GetId(2)])

        self.points = self.pdata.GetPoints()
        self.polys = self.pdata.GetPolys()
        self.numTriangles = self.polys.GetNumberOfCells()

        shp = (self.numTriangles, self.numTriangles)
        self.gMat = numpy.zeros(shp, numpy.float64)

        self.order = order
        self.__computeMatrices()

        self.normalEJumpName = 'normal_electric_field_jump'

    def getVtkPolyData(self):
        """
        Get the (modified) vtkPolyData object
        @return object
        """
        return self.pdata

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

    def __computeMatrices(self):

        self.__computeDiagonalTerms()
        self.__computeOffDiagonalTerms()

    def getGreenMatrix(self):
        """
        Return the Green function matrix
        @return matrix
        """
        return self.gMat

    def getPoints(self):
        """
        Get grid points
        @return points
        """
        points = self.pdata.GetPoints()
        numPoints = points.GetNumberOfPoints()
        res = numpy.zeros((numPoints, 3), numpy.float64)
        for i in range(numPoints):
            res[i, :] = points.GetPoint(i)
        return res

    def getCells(self):
        """
        Get cell connectivity
        @return cells
        """
        polys = self.pdata.GetPolys()
        numCells = polys.GetNumberOfCells()
        res = numpy.zeros((numCells, 3), numpy.int)
        polys.InitTraversal()
        ptIds = vtk.vtkIdList()
        for i in range(numCells):
            polys.GetNextCell(ptIds)
            res[i, :] = ptIds.GetId(0), ptIds.GetId(1), ptIds.GetId(2)
        return res
