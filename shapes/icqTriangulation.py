#!/usr/bin/env python

from __future__ import print_function
import numpy
import vtk


class Triangulation:

    def __init__(self):
        """
        Constructor
        """
        self.points = vtk.vtkPoints()
        self.polydata = vtk.vtkPolyData()
        self.delny = vtk.vtkDelaunay3D()

        self.polydata.SetPoints(self.points)
        if vtk.VTK_MAJOR_VERSION >= 6:
            self.delny.SetInputData(self.polydata)
        else:
            self.delny.SetInput(self.polydata)

    def setInputPoints(self, points,):
        """
        Set input points
        @param points points
        """
        numPoints = points.shape[0]
        for i in range(numPoints):
            self.points.InsertPoint(i, points[i, 0],
                                    points[i, 1], points[i, 2])

    def triangulate(self):
        """
        Triangulate
        """
        self.delny.Update()

    def getVTKUnstructuredGrid(self):
        """
        Get the VTK unstructured grid
        @return VTK object
        """
        return self.delny.GetOutput()

    def getCells(self):
        """
        Get cells
        @return list of vertex indices
        """
        cells = self.delny.GetOutput().GetCells()
        numCells = cells.GetNumberOfCells()

        npts = 4  # number of points per cell
        cellArr = numpy.zeros((numCells, npts), numpy.int)
        ptIds = vtk.vtkIdList()
        for i in range(numCells):
            cells.GetCell(i, ptIds)
            for j in range(npts):
                cellArr[i, j] = int(ptIds.GetId(j))

        return cellArr

##############################################################################


def test():

    points = numpy.array([(0, 0, 0), (0, 0, 1), (0, 1, 0),
                          (0, 1, 1), (1, 0, 0), (1, 0, 1),
                          (1, 1, 0), (1, 1, 1)], numpy.float64)
    tri = Triangulation()
    tri.setInputPoints(points)
    tri.triangulate()
    print(tri.getVTKUnstructuredGrid())
    cells = tri.getCells()
    print('cells = ', cells)

if __name__ == '__main__':
    test()
