#!/usr/bin/env python

from __future__ import print_function
import vtk
import numpy
from icqsol.bem.icqBaseSolver import BaseSolver
from icqsol.util.icqDataFetcher import getArrayIndexFromNameAndProjectOntoCells

class LaplaceSolver(BaseSolver):

    def __init__(self, pdata, max_edge_length, order=5):
        """
        Constructor
        @param pdata instance of vtkPolyData
        @param max_edge_length maximum edge length, used to turn
                               polygons into triangles
        """
        super(LaplaceSolver, self).__init__(pdata, max_edge_length, order)
        self.normalEJumpName = 'normal_electric_field_jump'

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


    def computeNormalElectricFieldJump(self, potName='v'):
        """
        Get the jump of the normal electric field - dv/dn
        from potential v
        @param potName name of the potential field in the vtkPolyData object
        @return response
        """
        # Has the potential been set?
        potIndex = getArrayIndexFromNameAndProjectOntoCells(self.pdata, potName)
        if potIndex < 0:
            msg = 'ERROR: could not find any cell field named {0}!'.format(potName)
            raise RuntimeError(msg)
        
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
        print('order = ', order)
        print('g matrix: ', lslm.getGreenMatrix())


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
        print('order = ', order)
        print('g matrix: ', lslm.getGreenMatrix())


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
        print('order = ', order)
        print('g matrix: ', lslm.getGreenMatrix())

if __name__ == '__main__':
    testSingleTriangle()
    testTwoTrianglesCoplanar()
    testTwoTriangles()
