#!/usr/bin/env python

import math
import numpy
import vtk

class CoarsenSurface:

    EPS = 1.23456789e-10

    def __init__(self, pdata):
        """
        Constructor
        @param pdata vtkPolyData instance
        """
        self.polydata = vtk.vtkPolyData()
        self.polydata.DeepCopy(pdata)
        self.points = self.polydata.GetPoints()

        # reuired in order to get the cell Ids sharing an edge
        self.polydata.BuildLinks()

        # will need to be able to interpolate the nodal data to the
        # new vertices
        self.pointData = self.polydata.GetPointData()
        self.numPointData = self.pointData.GetNumberOfArrays()

        # polygons
        self.polys = self.polydata.GetPolys()

        # number of polygons
        self.numPolys = self.polys.GetNumberOfCells()

        # polygon areas -- polygons will be sorted according to 
        # their polygon areas
        self.polyAreas = numpy.zeros((self.numPolys,), numpy.float64)
        self.polys.InitTraversal()
        ptIds = vtk.vtkIdList()
        for polyId in range(self.numPolys):
            self.polys.GetNextCell(ptIds)
            self.polyAreas[polyId] = self.getPolygonArea(ptIds)

        # cell indices sorted by increasing polygon areas
        self.sortedPolyIndices = numpy.argsort(self.polyAreas)

        # required so we can get the connectivity between points and 
        # cells
        self.polydata.BuildLinks()

    def getVtkPolyData(self):
        """
        Get the vtkPolyData instance associated with the grid
        @return vtkPolyData instance
        """
        return self.polydata

    def getPolygonArea(self, ptIds):
        """
        Compute the (scalar) area of a polygon
        @param ptIds list of point indices
        @return area
        """
        area = numpy.zeros((3,), numpy.float64)
        p0 = numpy.array(self.points.GetPoint(ptIds.GetId(0)))
        numPts = ptIds.GetNumberOfIds()
        for i in range(1, numPts - 1):
            p1 = numpy.array(self.points.GetPoint(ptIds.GetId(i    )))
            p2 = numpy.array(self.points.GetPoint(ptIds.GetId(i + 1)))
            area += numpy.cross(p1 - p0, p2 - p0)
        return numpy.linalg.norm(area)

    def colapsePolygon(self, polyId):
        """
        Colapse the vertices of the cell and adjust the 
        neighboring polygons' area
        @param polyId Id of the polygon
        """
        ptIds = vtk.vtkIdList()
        neighPolyIds = vtk.vtkIdList()
        self.polydata.GetCellPoints(polyId, ptIds)

        # iterate over the points to compute the barycenter
        barycenter = numpy.zeros((3,), numpy.float64)
        numPts = ptIds.GetNumberOfIds()
        for i in range(numPts):
            barycenter += self.points.GetPoint(ptIds.GetId(i))
        barycenter /= float(numPts)

        pId = vtk.vtkIdList()
        pId.SetNumberOfIds(1)

        # move each vertex of polyId to the barycenter position
        for i in range(numPts):

            # Id of this point
            pI = ptIds.GetId(i)

            # move this point
            self.points.SetPoint(pI, barycenter)



            # get a list of the polys that have this vertex and 
            # correct their area
            pId.SetId(pI)
            self.polydata.GetCellNeighbors(polyId, pId, neighPolyIds)
            numNeigh = neighPolyIds.GetNumberOfIds()
            for j in range(numNeigh):
                neighPolyId = neighPolyIds.GetId(j)
                # correct the area
                area = self.getPolygonArea(neighPolyId)
                self.polyAreas[neighPolyId] = area

        # this poly has now zero area
        self.polyAreas[polyId] = 0.

        # reset the nodal field values to account for the vertices 
        # having moved
        self.averagePointData(ptIds)

    def coarsen(self, min_cell_area = 1.e-10):
        """
        Coarsen surface by colapsing small polygons
        @param min_cell_area cell area tolerance
        @note operation is in place
        """
        validPolyIndices = numpy.where(
            self.polyAreas[self.sortedPolyIndices] > self.EPS)
        if len(validPolyIndices) == 0:
            # no valid polygons with area > 0
            # nothing to do 
            return 

        polyId = validPolyIndices[0]
        go = True

        while go and polyId < self.numPolys:

            polyArea = self.polyAreas[polyId]
            if polyArea > min_cell_area:
                # done
                go = False
                break

            else:

                # colapse polygon. WILL NEED TO DO SOMETHING 
                # ABOUT VERTICES THAT ARE AT THE BOUNDARY
                self.colapsePolygon(polyId)

                # re-order the polygons since some 
                # areas have changed
                self.sortedPolyIndices = numpy.argsort(self.polyAreas)
                validPolyIndices = numpy.where(
                    self.polyAreas[self.sortedPolyIndices] > self.EPS)
                if len(validPolyIndices) == 0:
                    go = False
                    break

                # new poly index
                polyId = validPolyIndices[0]
 
    def averagePointData(self, ptIds):
        """
        Average the field at the point locations and set the nodal field values 
        to the average value
        @param ptIds set of point Ids
        """
        for el in range(self.numPointData):
            arr = self.pointData.GetArray(el)
            numComps = arr.GetNumberOfComponents()
            vals = numpy.zeros((numComps,), numpy.float64)
            baryVals = numpy.zeros((numComps,), numpy.float64)
            # mid cell values
            for ptId in ptIds:
                vals[:] = arr.GetTuple(ptId)
                baryVals += vals
            n = len(ptIds)
            baryVals /= float(n)
            # set the field values to the mid cell values
            for j in range(n):
                arr.SetTuple(ptIds.GetId(j), baryVals)

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
