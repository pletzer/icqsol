#!/usr/bin/env python

from __future__ import print_function
import vtk
import numpy
from icqsol.bem.icqBaseLaplaceSolver import BaseLaplaceSolver


class PoissonSolver(BaseLaplaceSolver):

    def __init__(self, pdata, max_edge_length, order=5):
        """
        Constructor
        @param pdata instance of vtkPolyData
        @param max_edge_length maximum edge length, used to turn
                               polygons into triangles
        """
        #super(BaseLaplaceSolver, self).__init__(pdata, max_edge_length, order)
        BaseLaplaceSolver.__init__(self, pdata, max_edge_length, order)
        self.responseName = 'v'
        self.sourceName = 'charge'

    def computeResponseField(self):
        """
        Compute the response field, in this case the potential due to a charge source
        @return response
        """
        
        srcIndex = self.getSourceArrayIndex()
        src = self.getSourceArray(srcIndex)

        # Get the response matrix.
        gMat = self.getGreenMatrix()

        # Compute the response.
        rsp = numpy.dot(gMat, src)

        self.addResponseField(rsp)

        return rsp

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
        lslm = PoissonSolver(pdata, max_edge_length=1000.)
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
        lslm = PoissonSolver(pdata,
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
        lslm = PoissonSolver(pdata,
                             max_edge_length=1000.,
                             order=order)
        print('order = ', order)
        print('g matrix: ', lslm.getGreenMatrix())

if __name__ == '__main__':
    testSingleTriangle()
    testTwoTrianglesCoplanar()
    testTwoTriangles()
