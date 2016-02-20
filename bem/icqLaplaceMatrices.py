#!/usr/bin/env python

import vtk
import numpy
from icqsol.shapes.icqRefineSurface import RefineSurface
from icqsol.bem.icqPotentialIntegrals import PotentialIntegrals
from icqQuadrature import gaussPtsAndWeights
import pkg_resources
from ctypes import cdll, POINTER, byref, c_void_p, c_double, c_long

FOUR_PI = 4. * numpy.pi


class LaplaceMatrices:

    def __init__(self, pdata, max_edge_length, order=5):
        """
        Constructor
        @param pdata instance of vtkPolyData
        @param max_edge_length maximum edge length, used to turn
                               polygons into triangles
        """

        libName = pkg_resources.resource_filename('icqsol', 'icqLaplaceMatricesCpp.so')
        self.lib = cdll.LoadLibrary(libName)

        self.normalEJumpName = 'normal_electric_field_jump'

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
        
    def getArrayIndexFromName(self, data, name):
        """
        Get the array index from its name
        @param data either a vtkCellData or a vtkPointData object
        @param name name
        @return index >= 0 if name exists or -1 if it does not
        """
        numArrays = data.GetNumberOfArrays()
        index = -1
        for i in range(numArrays):
            arr = data.GetArray(i)
            if arr.GetName() == name:
                index = i
                break
        return index

    def getPotentialArrayIndexFromName(self, name):
        """
        Get the potential array index from its name, if the array 
        is a point data type then project onto cells
        @param name name
        @return index >= 0 if name exists or -1 if it does not
        """
        cellData = self.pdata.GetCellData()
        pointData = self.pdata.GetPointData()
        index = self.getArrayIndexFromName(cellData, name)
        if index < 0:
            # Maybe a point array?
            index2 = self.getArrayIndexFromName(pointData, name)

            if index2 >= 0:
                # Project from points to cells
                pointArr = pointData.GetArray(index2)
                numCells = self.pdata.GetNumberOfPolys()
                cellArr = vtk.vtkDoubleArray()
                cellArr.SetName(name) # same name as the point array
                cellArr.SetNumberOfComponents(1)
                cellArr.SetNumberOfTuples(numCells)
                cells = self.pdata.GetPolys()
                ptIds = vtk.vtkIdList()
                cells.InitTraversal()
                for i in range(numCells):
                    cells.GetNextCell(ptIds)
                    # We know for sure that all cells are triangles
                    ia, ib, ic = ptIds.GetId(0), \
                        ptIds.GetId(1), ptIds.GetId(2)
                    va = pointArr.GetComponent(ia, 0)
                    vb = pointArr.GetComponent(ib, 0)
                    vc = pointArr.GetComponent(ic, 0)
                    cellArr.SetComponent(i, 0, (va + vb + vc)/3.)
                # Add the cell array
                cellData.AddArray(cellArr)
                return self.getPotentialArrayIndexFromName(name)

        return index

    def setPotentialFromExpression(self, expression, potName='v'):
        """
        Set the potential from expression
        @param potName name of the potential field saved in vtkPolyData 
        @param expression expression of x, y, and z
        """
        from math import sqrt, pi, sin, cos, tan, log, exp
        
        n = self.pdata.GetNumberOfPolys()
        potentialData = vtk.vtkDoubleArray()
        potentialData.SetNumberOfComponents(1)
        potentialData.SetNumberOfTuples(n)
        potentialData.SetName(potName)
        midPoint = numpy.zeros((3,), numpy.float64)
        ptIds = vtk.vtkIdList()
        cells = self.pdata.GetPolys()
        cells.InitTraversal()
        for i in range(n):
            cell = cells.GetNextCell(ptIds)
            npts = ptIds.GetNumberOfIds()
            midPoint *= 0 # reset
            for j in range(npts):
                midPoint += self.points.GetPoint(ptIds.GetId(j))
            midPoint /= float(npts)
            x, y, z = midPoint
            v = eval(expression)
            potentialData.SetTuple(i, [v])
        self.pdata.GetCellData().AddArray(potentialData)

    def setNormalElectricFieldJumpName(self, name):
        """
        Set the name of the normal electric field jump
        @param name name
        """
        self.normalEJumpName = name

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
                # Three triangles having oberver point as one corner
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

    def computeNormalElectricFieldJump(self, potName='v'):
        """
        Get the jump of the normal electric field - dv/dn
        from potential v
        @param potName name of the potential field in the vtkPolyData object
        @return response
        """
        # Has the potential been set?
        potIndex = self.getPotentialArrayIndexFromName(potName)
        if potIndex < 0:
            raise RuntimeError, \
                'ERROR: could not find any cell field named {0}!'.format(potName)
        potArray = self.pdata.GetCellData().GetArray(potIndex)

        # Set the potential.
        n = self.pdata.GetNumberOfPolys()
        v = numpy.zeros((n,), numpy.float64)
        for i in range(n):
            v[i] = potArray.GetComponent(i, 0)

        gMat = self.getGreenMatrix()

        normalEJump = -numpy.linalg.solve(gMat, v)

        # Add normal electric field jump
        normalEJumpData = vtk.vtkDoubleArray()
        normalEJumpData.SetNumberOfComponents(1)
        normalEJumpData.SetNumberOfTuples(n)
        normalEJumpData.SetName(self.normalEJumpName)
        for i in range(n):
            normalEJumpData.SetTuple(i, [normalEJump[i]])
        self.pdata.GetCellData().AddArray(normalEJumpData)

        return normalEJump

###############################################################################


def testSingleTriangle():

    "Single triangle"

    h = 0.1
    # create set of points
    points = vtk.vtkPoints()
    points.SetNumberOfPoints(3)
    points.SetPoint(0, [1., -1.*h/3., -1.*h/3.])
    points.SetPoint(1, [1., 2.*h/3., -1.*h/3.])
    points.SetPoint(2, [1., -1.*h/3., 2.*h/3.])

    # create vtkPolyData object
    pdata = vtk.vtkPolyData()
    pdata.SetPoints(points)
    ptIds = vtk.vtkIdList()
    ptIds.SetNumberOfIds(3)
    ptIds.SetId(0, 0)
    ptIds.SetId(1, 1)
    ptIds.SetId(2, 2)
    pdata.Allocate(1, 1)
    pdata.InsertNextCell(vtk.VTK_POLYGON, ptIds)

    for order in range(1, 6):
        lslm = LaplaceMatrices(pdata, max_edge_length=1000.)
        print 'order = ', order
        print 'g matrix: ', lslm.getGreenMatrix()


def testTwoTrianglesCoplanar():

    "Two triangles"

    # create set of points
    points = vtk.vtkPoints()
    points.SetNumberOfPoints(4)
    points.SetPoint(0, [0., 0., 0.])
    points.SetPoint(1, [1., 0., 0.])
    points.SetPoint(2, [0., 1., 0.])
    points.SetPoint(3, [1., 1., 0.])

    # create vtkPolyData object
    pdata = vtk.vtkPolyData()
    pdata.SetPoints(points)

    pdata.Allocate(2, 1)
    ptIds = vtk.vtkIdList()
    ptIds.SetNumberOfIds(3)

    ptIds.SetId(0, 0)
    ptIds.SetId(1, 1)
    ptIds.SetId(2, 2)
    pdata.InsertNextCell(vtk.VTK_POLYGON, ptIds)
    ptIds.SetId(0, 1)
    ptIds.SetId(1, 3)
    ptIds.SetId(2, 2)
    pdata.InsertNextCell(vtk.VTK_POLYGON, ptIds)

    for order in range(1, 6):
        lslm = LaplaceMatrices(pdata,
                               max_edge_length=1000.,
                               order=order)
        print 'order = ', order
        print 'g matrix: ', lslm.getGreenMatrix()


def testTwoTriangles():

    "Two triangles"

    # create set of points
    points = vtk.vtkPoints()
    points.SetNumberOfPoints(4)
    points.SetPoint(0, [0., 0., 0.])
    points.SetPoint(1, [1., 0., 0.])
    points.SetPoint(2, [0., 1., 0.])
    points.SetPoint(3, [0., 0., 1.])

    # create vtkPolyData object
    pdata = vtk.vtkPolyData()
    pdata.SetPoints(points)

    pdata.Allocate(2, 1)
    ptIds = vtk.vtkIdList()
    ptIds.SetNumberOfIds(3)

    ptIds.SetId(0, 0)
    ptIds.SetId(1, 1)
    ptIds.SetId(2, 3)
    pdata.InsertNextCell(vtk.VTK_POLYGON, ptIds)
    ptIds.SetId(0, 0)
    ptIds.SetId(1, 3)
    ptIds.SetId(2, 2)
    pdata.InsertNextCell(vtk.VTK_POLYGON, ptIds)

    for order in range(1, 6):
        lslm = LaplaceMatrices(pdata,
                               max_edge_length=1000.,
                               order=order)
        print 'order = ', order
        print 'g matrix: ', lslm.getGreenMatrix()

if __name__ == '__main__':
    testSingleTriangle()
    testTwoTrianglesCoplanar()
    testTwoTriangles()
