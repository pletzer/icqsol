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
        self.polydata = pdata

        # save the points and point data as separate class members
        # so as not to pollute pdata
        self.points = vtk.vtkPoints()
        self.points.DeepCopy(pdata.GetPoints())

        self.pointData = {}
        pd = pdata.GetPointData()
        for i in range(pd.GetNumberOfArrays()):
            arr = pd.GetArray(i)
            name = arr.GetName()
            self.pointData[name] = vtk.vtkDoubleArray()
            self.pointData[name].DeepCopy(arr)

        self.cellData = {}
        cd = pdata.GetCellData()
        for i in range(cd.GetNumberOfArrays()):
            name = cd.GetArray(i).GetName()
            self.cellData[name] = vtk.vtkDoubleArray()
            self.cellData[name].SetName(name)

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

        pd = self.polydata.GetPointData()
        numArrays = pd.GetNumberOfArrays()

        points = self.polydata.GetPoints()

        polys = self.polydata.GetPolys()
        p = numpy.zeros((3,), numpy.float64)
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

            if self.getPolygonArea(ptIds) < min_cell_area:

                # number of points spanning the cell (or number of edges)
                n = ptIds.GetNumberOfIds()

                # iterate over edges
                neighCellIds = set() # unique entries
                for iEdge in range(n):
                    ptId1, ptId2 = ptIds.GetId(iEdge), ptIds.GetId((iEdge + 1) % n)
                    self.polydata.GetCellEdgeNeighbors(polyId, ptId1, ptId2, cellIds)
                    for j in range(cellIds.GetNumberOfIds()):
                        neighCellIds.add(cellIds.getId(j))

                if len(neighCellIds) < n:
                    # must be a boundary cell, skip for the time being...
                    # might need to do something special here. Just want to 
                    # be sure the boundary does not move...
                    continue

                # compute the cell's gravity center
                barycenter = numpy.zeros((3,), numpy.float64)
                for j in range(n):
                    p[:] = points.GetPoint(ptId)
                    barycenter += p
                barycenter /= float(n)

                # interpolate the nodal fields to the barycenter locations
                for el in range(numArrays):
                    arr = pd.GetArray(el)
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

                # tag the cell for removal
                cellIdsToRemove.append(cellId)
        
        # remove the tagged, zero-area cells
        for cellId in cellIdsToRemove:
            self.polydata.DeleteCell(cellId)

        # now remove
        self.polydata.RemoveDeletedCells()

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
