#!/usr/bin/env python

import math
import numpy
import vtk

class CoarsenSurface:

    def __init__(self, pdata):
        """
        Constructor
        @param pdata vtkPolyData instance
        """
        self.polydata = vtk.vtkPolyData()
        self.polydata.DeepCopy(pdata)
        self.points = self.polydata.GetPoints()

        self.polydata.BuildLinks()

        self.pointData = self.polydata.GetPointData()
        self.numPointData = self.pointData.GetNumberOfArrays()

        # compute the area of each cell
        polys = self.polydata.GetPolys()
        numPolys = polys.GetNumberOfCells()
        self.polyAreas = numpy.zeros((numPolys,), numpy.float64)
        polys.InitTraversal()
        ptIds = vtk.vtkIdList()
        for polyId in range(numPolys):
            polys.GetNextCell(ptIds)
            self.polyAreas[polyId] = self.getPolygonArea(ptIds)

        # sort the cell by deacreasing cell areas
        self.sortedPolyIndices = numpy.fliplr([numpy.argsort(self.polyAreas)])[0]


    def getVtkPolyData(self):
        """
        Get the vtkPolyData instance associated with the grid
        @return vtkPolyData instance
        """
        return self.polydata

    def coarsen(self, min_cell_area = 1.e-10):
        """
        Coarsen surface by removing very small cells
        @param min_cell_area cell area tolerance
        @note operation is in place
        """

        points = self.polydata.GetPoints()

        polys = self.polydata.GetPolys()
        numPolys = polys.GetNumberOfCells()
        if numPolys <= 1:
            # must have at least one polygon
            return

        ptIds = vtk.vtkIdList()

        # holds the cell Ids to delete
        cellIdsToRemove = []
        cellIds = vtk.vtkIdList()

        # iterate over cells
        polys.InitTraversal()
        for polyId in range(numPolys):

            polys.GetNextCell(ptIds)

            polygonArea = self.getPolygonArea(ptIds)
            if polygonArea < min_cell_area:

                # number of points spanning the cell 
                # (also equal to number of edges)
                n = ptIds.GetNumberOfIds()

                # iterate over edges
                neighCellIds = set() # unique entries
                for iEdge in range(n):
                    ptId1, ptId2 = ptIds.GetId(iEdge), ptIds.GetId((iEdge + 1) % n)
                    self.polydata.GetCellEdgeNeighbors(polyId, ptId1, ptId2, cellIds)
                    for j in range(cellIds.GetNumberOfIds()):
                        neighCellIds.add(cellIds.GetId(j))

                if len(neighCellIds) < n:
                    # must be a boundary cell, skip for the time being...
                    # might need to do something special here. Just want to 
                    # be sure the boundary does not move...
                    continue

                # move the points spanning the polygon to the 
                # center of the polygon. This will essentially 
                # reduce the cell to a set of points which are 
                # on top of each other.
                self.movePointsToBaryCenter(ptIds)

                # tag the cell for removal since it now has zero 
                # area
                cellIdsToRemove.append(polyId)
        
        # remove the tagged, zero-area cells
        for cellId in cellIdsToRemove:
            self.polydata.DeleteCell(cellId)

        # now remove
        self.polydata.RemoveDeletedCells()

    def movePointsToBaryCenter(self, ptIds):
        """
        Move points to barycenter location
        @param ptIds point Ids
        """

        p = numpy.zeros((3,), numpy.float64)
        n = ptIds.GetNumberOfIds()
        # compute the cell's gravity center
        barycenter = numpy.zeros((3,), numpy.float64)
        for j in range(n):
            p[:] = self.points.GetPoint(ptIds.GetId(j))
            barycenter += p
        barycenter /= float(n)

        # interpolate the nodal fields to the barycenter locations.
        # This simply amounts to averaging the nodal fields
        for el in range(self.numPointData):
            arr = self.pointData.GetArray(el)
            numComps = arr.GetNumberOfComponents()
            vals = numpy.zeros((numComps,), numpy.float64)
            baryVals = numpy.zeros((numComps,), numpy.float64)
            # mid cell values
            for j in range(n):
                vals[:] = arr.GetTuple(ptIds.GetId(j))
                baryVals += vals
            baryVals /= float(n)
            # set the field values to the mid cell values
            for j in range(n):
                arr.SetTuple(ptIds.GetId(j), baryVals)


    def getPolygonArea(self, polyPtIds):
        """
        Compute the (scalar) area of a polygon
        @param polyPtIds list of point indices
        @return area
        """
        area = numpy.zeros((3,), numpy.float64)
        p0 = numpy.array(self.points.GetPoint(polyPtIds.GetId(0)))
        numPolyPts = polyPtIds.GetNumberOfIds()
        for i in range(1, numPolyPts - 1):
            p1 = numpy.array(self.points.GetPoint(polyPtIds.GetId(i    )))
            p2 = numpy.array(self.points.GetPoint(polyPtIds.GetId(i + 1)))
            area += numpy.cross(p1 - p0, p2 - p0)
        return numpy.linalg.norm(area)

    def removeDegeneratePoints(self, polyPtIds):
        """
        Remove degenerate points
        @param polyPtIds list of point indices (modified on output)
        """
        indicesToDelete = []
        numPolyPts = len(polyPtIds)
        for i in range(numPolyPts):
            p0 = numpy.array(self.points.GetPoint(polyPtIds[i]))
            p1 = numpy.array(self.points.GetPoint(polyPtIds[(i + 1)%numPolyPts]))
            if numpy.linalg.norm(p1 - p0) < 1.e-15:
                indicesToDelete.append(i)
        for j in range(len(indicesToDelete) - 1, -1, -1):
            i = indicesToDelete[j]
            del polyPtIds[i]

##############################################################################


def printVtkPolyData(pdata):

    points = pdata.GetPoints()
    polys = pdata.GetPolys()
    ptIds = vtk.vtkIdList()
    numPolys = polys.GetNumberOfCells()
    print 'Number of polygons: {0}'.format(numPolys)
    polys.InitTraversal()
    for i in range(numPolys):
        numPts = ptIds.GetNumberOfIds()
        print '\tCell {0} has {1} points: '.format(i, numPts)
        for j in range(numPts):
            ptId = ptIds.GetId(j)
            pt = points.GetPoint(ptId)
            print '\t\t{0} -> {1}'.format(ptId, pt)


def testNoCoarsening():

    points = vtk.vtkPoints()
    points.InsertNextPoint((0., 0., 0.))
    points.InsertNextPoint((1., 0., 0.))
    points.InsertNextPoint((1., 1., 0.))

    pdata = vtk.vtkPolyData()
    pdata.SetPoints(points)

    ptIds = vtk.vtkIdList()
    ptIds.InsertNextId(0)
    ptIds.InsertNextId(1)
    ptIds.InsertNextId(2)
    pdata.Allocate(1, 1)
    pdata.InsertNextCell(vtk.VTK_POLYGON, ptIds)

    rs = CoarsenSurface(pdata)
    rs.coarsen(min_cell_area=1.e-5)
    pdata2 = rs.getVtkPolyData()
    assert(pdata2.GetNumberOfPolys() == 1)


def testAddingThreePointsThenMore():

    points = vtk.vtkPoints()
    points.InsertNextPoint((0., 0., 0.))
    points.InsertNextPoint((2., 0., 0.))
    points.InsertNextPoint((2., 1., 0.))

    pdata = vtk.vtkPolyData()
    pdata.SetPoints(points)

    ptIds = vtk.vtkIdList()
    ptIds.InsertNextId(0)
    ptIds.InsertNextId(1)
    ptIds.InsertNextId(2)
    pdata.Allocate(1, 1)
    pdata.InsertNextCell(vtk.VTK_POLYGON, ptIds)

    rs = CoarsenSurface(pdata)
    rs.coarsen(min_cell_area=1.e-5)
    assert(rs.getVtkPolyData().GetNumberOfPolys() == 1)


def testStartingWithTwoCells():

    points = vtk.vtkPoints()
    points.InsertNextPoint((0., 0., 0.))
    points.InsertNextPoint((2., 0., 0.))
    points.InsertNextPoint((2., 1., 0.))
    points.InsertNextPoint((0., 1., 0.))

    pdata = vtk.vtkPolyData()
    pdata.SetPoints(points)

    pdata.Allocate(2, 1)
    ptIds = vtk.vtkIdList()
    ptIds.SetNumberOfIds(3)
    ptIds.SetId(0, 0)
    ptIds.SetId(1, 1)
    ptIds.SetId(2, 2)
    pdata.InsertNextCell(vtk.VTK_POLYGON, ptIds)
    ptIds.SetId(0, 2)
    ptIds.SetId(1, 3)
    ptIds.SetId(2, 0)
    pdata.InsertNextCell(vtk.VTK_POLYGON, ptIds)

    rs = CoarsenSurface(pdata)
    rs.coarsen(min_cell_area=1.e-5)
    assert(rs.getVtkPolyData().GetNumberOfPolys() == 2)


if __name__ == '__main__':
    testNoCoarsening()
    testAddingThreePointsThenMore()
    testStartingWithTwoCells()
