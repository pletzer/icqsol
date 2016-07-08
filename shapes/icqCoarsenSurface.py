#!/usr/bin/env python

from __future__ import print_function
from operator import itemgetter
import math
import numpy
import vtk
import math

class CoarsenSurface:

    EPS = 1.23456789e-10
    TWOPI = 2.0 * math.pi

    def __init__(self, pdata):
        """
        Constructor
        @param pdata vtkPolyData instance
        """
        self.polydata = vtk.vtkPolyData()
        self.polydata.DeepCopy(pdata)

        # will need to be able to interpolate the nodal data to the
        # new vertices
        self.pointData = self.polydata.GetPointData()
        self.numPointData = self.pointData.GetNumberOfArrays()

        # required so we can get the connectivity between points and 
        # cells
        self.polydata.BuildLinks()

    def getVtkPolyData(self):
        """
        Get the vtkPolyData instance associated with the grid
        @return vtkPolyData instance
        """
        return self.polydata

    def getPolygonVectorArea(self, ptIds):
        """
        Get the polygon vector area 
        @param ptIds boundary points of the cell
        @return vector
        """
        points = self.polydata.GetPoints()
        areaVec = numpy.zeros((3,), numpy.float64)
        p0 = numpy.array(points.GetPoint(ptIds.GetId(0)))
        numPts = ptIds.GetNumberOfIds()
        for i in range(1, numPts - 1):
            p1 = numpy.array(points.GetPoint(ptIds.GetId(i    )))
            p2 = numpy.array(points.GetPoint(ptIds.GetId(i + 1)))
            areaVec += numpy.cross(p1 - p0, p2 - p0)
        return areaVec

    def getPolygonArea(self, ptIds):
        """
        Compute the (scalar) area of a polygon
        @param ptIds list of point indices
        @return area
        """
        areaVec = self.getPolygonVectorArea(ptIds)
        return numpy.linalg.norm(areaVec)

    def getPolygonNormalVector(self, ptIds):
        """
        Compute the normal vector of the polygon
        @param ptIds list of point indices
        @return area
        """
        areaVec = self.getPolygonVectorArea(ptIds)
        area = numpy.linalg.norm(areaVec)
        return areaVec / area

    def getTotalAngle(self, ptId):
        """
        Compute the sum of the angles between this point and the surrounding edges
        @param ptId point Id
        @return total angle in radiants, should be about 2*pi for an internal point
        """
        cellIds = vtk.vtkIdList()

        # get the cell Ids sharing this point
        ptIds = vtk.vtkIdList()
        self.polydata.GetPointCells(ptId, cellIds)

        pt0 = numpy.zeros((3,), numpy.float64)
        edge1 = numpy.zeros((3,), numpy.float64)
        edge2 = numpy.zeros((3,), numpy.float64)
        points = self.polydata.GetPoints()

        # get the point coorindates of ptId
        points.GetPoint(ptId, pt0)

        # iterate over the cells sharing ptId
        angle = 0.
        for iCell in range(cellIds.GetNumberOfIds()):

            # get the points of this cell
            self.polydata.GetCellPoints(cellIds.GetId(iCell), ptIds)

            # find the two edges incident to ptId
            n = ptIds.GetNumberOfIds()
            for j in range(n):
                p1 = ptIds.GetId(j)
                p2 = ptIds.GetId((j + 1) % n)
                if p1 == ptId:
                    self.polydata.GetPoints().GetPoint(p2, edge1)
                elif p2 == ptId:
                    self.polydata.GetPoints().GetPoint(p1, edge2)

            # the two edges
            edge1 -= pt0
            edge2 -= pt0

            # add the angle between pt1 and pt2
            crossProduct = numpy.linalg.norm(numpy.cross(edge1, edge2))
            dotProduct = numpy.dot(edge1, edge2)
            angle += math.atan2(crossProduct, dotProduct)

        return angle

    def collapsePolygon(self, cellId):
        """
        Collapse the vertices of a cell
        @param cellId Id of the polygon
        """
        # points of the cell
        ptIds = vtk.vtkIdList()
        self.polydata.GetCellPoints(cellId, ptIds)
        npts = ptIds.GetNumberOfIds()

        # typically the center of the cell, in some 
        # cases the center of an edge
        center = numpy.zeros( (npts,), numpy.float64 )

        # coordinates of the poinrt
        pt = numpy.zeros( (npts,), numpy.float64 )

        points = self.polydata.GetPoints()
        pointsToMove = []

        # compute the polygon center
        for i in range(npts):
            ptId = ptIds.GetId(i)
            # sum of the angles must be 2*pi for internal points
            totalAngle = self.getTotalAngle(ptId)
            if abs(totalAngle - self.TWOPI) < 0.01:
                # internal, ie non-boundary point
                points.GetPoint(ptId, pt)
                center += pt
                # add point to the list of points to move
                pointsToMove.append(ptId)
            #else: print('*** not an internal node: cellId = {} self.getTotalAngle(ptId) = {}'.format(cellId, totalAngle))
        
        n = len(pointsToMove)
        if n > 0:
            center /= float(len(pointsToMove))

        # move the selected points to the center
        for ptId in pointsToMove:
            points.SetPoint(ptId, center)

        # average the nodal data at the new vertex location
        self.averagePointData(ptIds)

    def coarsen(self, min_cell_area = 1.e-10):
        """
        Coarsen surface by collapsing small polygons
        @param min_cell_area cell area tolerance
        @note operation is in place
        """

        polys = self.polydata.GetPolys()
        numPolys = polys.GetNumberOfCells()
        ptIds = vtk.vtkIdList()
        polys.InitTraversal()
        for cellId in range(numPolys):
            polys.GetNextCell(ptIds)
            area = abs(self.getPolygonArea(ptIds))
            if area < min_cell_area and area > self.EPS:
                self.collapsePolygon(cellId)

        self.deleteZeroPolys()

    def deleteZeroPolys(self):
        """
        Delete all the polygons whose area is zero
        """
        ptIds = vtk.vtkIdList()
        numPolys = self.polydata.GetPolys().GetNumberOfCells()
        for polyId in range(numPolys):
            self.polydata.GetCellPoints(polyId, ptIds)
            if abs(self.getPolygonArea(ptIds)) <= self.EPS:
                self.polydata.DeleteCell(polyId)
            
        self.polydata.RemoveDeletedCells()
        self.polydata.BuildLinks() # not sure if this is required
 
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
            numPts = ptIds.GetNumberOfIds()
            for i in range(numPts):
                ptId = ptIds.GetId(i)
                vals[:] = arr.GetTuple(ptId)
                baryVals += vals
            baryVals /= float(numPts)
            # set the field values to the mid cell values
            for j in range(numPts):
                arr.SetTuple(ptIds.GetId(j), baryVals)

##############################################################################


def printVtkPolyData(pdata):

    points = pdata.GetPoints()
    polys = pdata.GetPolys()
    ptIds = vtk.vtkIdList()
    numPolys = polys.GetNumberOfCells()
    print('Number of polygons: {0}'.format(numPolys))
    polys.InitTraversal()
    for i in range(numPolys):
        numPts = ptIds.GetNumberOfIds()
        print('\tCell {0} has {1} points: '.format(i, numPts))
        for j in range(numPts):
            ptId = ptIds.GetId(j)
            pt = points.GetPoint(ptId)
            print('\t\t{0} -> {1}'.format(ptId, pt))


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
