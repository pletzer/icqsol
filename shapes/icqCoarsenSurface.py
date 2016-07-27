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

        # old point id -> new point id
        self.degeneratePtIdMap = {}

        # required so we can get the connectivity between points and 
        # cells
        self.polydata.BuildLinks()
        self.polydata.BuildCells()

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
        @param return list of points that have been moved
        """
        # points of the cell
        ptIds = vtk.vtkIdList()
        self.polydata.GetCellPoints(cellId, ptIds)
        npts = ptIds.GetNumberOfIds()

        # typically the center of the cell, in some 
        # cases the center of an edge
        center = numpy.zeros( (3,), numpy.float64 )

        # coordinates of the poinrt
        pt = numpy.zeros( (npts,), numpy.float64 )

        points = self.polydata.GetPoints()
        pointsToMove = []
        
        # determine which points are internal/boundary
        internalPointIds = []
        boundaryPointIds = []
        for i in range(npts):
            ptId = ptIds.GetId(i)
            totalAngle = self.getTotalAngle(ptId)
            if abs(totalAngle - self.TWOPI) < 1.e-6: #0.01:
                internalPointIds.append(ptId)
            else:
                boundaryPointIds.append(ptId)
                    
        # compute the (central) point where the cell collapses to
        numBoundaryPoints = len(boundaryPointIds)
        if numBoundaryPoints == 0:
            # all points are internal, average the (internal) points
            for ptId in internalPointIds:
                points.GetPoint(ptId, pt)
                center += pt
            center /= float(len(internalPointIds))
            pointsToMove = internalPointIds
        #elif numBoundaryPoints > 1 and numBoundaryPoints < npts:
        elif numBoundaryPoints > npts and numBoundaryPoints < npts:  # not working well so turning off for the time being
            # average the boundary points and move all points to this position
            for ptId in boundaryPointIds:
                points.GetPoint(ptId, pt)
                center += pt
            center /= float(numBoundaryPoints)
            pointsToMove = internalPointIds + boundaryPointIds
        else:
            # likely all the points are boundary or there are more than two boundary
            # edges, nothing to do
            pass

        # move the selected points to the center
        for ptId in pointsToMove:
            points.SetPoint(ptId, center)

        # store in memory which point ids are degenerate
        for ptId in pointsToMove[1:]:
            self.degeneratePtIdMap[ptId] = pointsToMove[0]

        # average the nodal data at the new vertex location
        self.averagePointData(pointsToMove)

        return pointsToMove

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

                # collapse the cell and get the point Ids that have moved
                movedPointIds = self.collapsePolygon(cellId)

                # clean up the connectivity by removing degenerate point Ids
                self.removeDegeneratePointsInConnectivity(movedPointIds)

        self.deleteZeroPolys()

    def removeDegeneratePointsInConnectivity(self, movedPointIds):
        """
        Remove the degenerate point Ids from the connectivity 
        @param movedPointIds list of point Ids that have moved
        """
        neighCellIds = vtk.vtkIdList()
        neighCellPointIds = vtk.vtkIdList()
        if len(movedPointIds) > 1:
            # all other point Ids will be replaced by newPointId
            newPointId = movedPointIds[0]
            for oldPointId in movedPointIds[1:]:
                # get all the cells that contain oldPointId
                self.polydata.GetPointCells(oldPointId, neighCellIds)
                numNeighCells = neighCellIds.GetNumberOfIds()
                for neighCellId in range(numNeighCells):
                    # extract the ptIds of that cell, requires BuildCells 
                    # to be called
                    self.polydata.GetCellPoints(neighCellId, neighCellPointIds)
                    npts = neighCellPointIds.GetNumberOfIds()
                    for i in range(npts):
                        pi = neighCellPointIds.GetId(i)
                        if pi == oldPointId:
                            # replace point Id
                            self.polydata.ReplaceCellPoint(neighCellId, oldPointId, newPointId)

    def deleteZeroPolys(self):
        """
        Delete all the polygons whose areas are (nearly) zero
        """
        ptIds = vtk.vtkIdList()
        numPolys = self.polydata.GetPolys().GetNumberOfCells()
        for polyId in range(numPolys):
            self.polydata.GetCellPoints(polyId, ptIds)
            if abs(self.getPolygonArea(ptIds)) <= self.EPS:
                self.polydata.DeleteCell(polyId)
            
        self.polydata.RemoveDeletedCells()
        self.polydata.BuildLinks() # not sure if this is required
        self.polydata.BuildCells()
 
    def averagePointData(self, pids):
        """
        Average the field at the point locations and set the nodal field values 
        to the average value
        @param pids list of point Ids
        """
        numPts = len(pids)
        if numPts == 0: return

        for el in range(self.numPointData):
            arr = self.pointData.GetArray(el)
            numComps = arr.GetNumberOfComponents()
            vals = numpy.zeros((numComps,), numpy.float64)
            baryVals = numpy.zeros((numComps,), numpy.float64)
            # mid cell values
            for i in range(numPts):
                ptId = pids[i]
                vals[:] = arr.GetTuple(ptId)
                baryVals += vals
            baryVals /= float(numPts)
            # set the field values to the mid cell values
            for j in range(numPts):
                arr.SetTuple(pids[j], baryVals)

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
