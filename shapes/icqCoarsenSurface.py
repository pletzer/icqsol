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

        # point Ids attached to each polygon. As far as I know we 
        # cannot access the point Ids in random way -- must use
        # the traversal itereation
        self.poly2PointIds = []
        self.polys.InitTraversal()
        ptIds = vtk.vtkIdList()
        for polyId in range(self.numPolys):
            self.polys.GetNextCell(ptIds)
            self.polyAreas[polyId] = self.getPolygonArea(ptIds)
            numPts = ptIds.GetNumberOfIds()
            pids = []
            for j in range(numPts):
                pids.append(ptIds.GetId(j))
            self.poly2PointIds.append(tuple(pids))

        # sell indices sorted by increasing polygon areas
        self.sortedPolyIndices = numpy.argsort(self.polyAreas)

    def getVtkPolyData(self):
        """
        Get the vtkPolyData instance associated with the grid
        @return vtkPolyData instance
        """
        return self.polydata

    def getPolygonArea(self, ptIds):
        """
        Compute the (scalar) area of a polygon
        @param polyPtIds list of point indices
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

    def coarsen(self, min_cell_area = 1.e-10):
        """
        Coarsen surface by removing very small cells
        @param min_cell_area cell area tolerance
        @note operation is in place
        """

        # polygon index in the self.sortedPolyIndices array denoting 
        # the first valid polygon. That is indices smaller than 
        # firstValidPolyIndex should all have zero area
        firstValidPolyIndex = 0

        # flag denoting whether we need to colapse the smallest, non-zero
        # polygon
        colapsePoly = True

        # collection of cell Ids sharing an edge
        cellIds = vtk.vtkIdList()

        while colapsePoly and firstValidPolyIndex < self.numPolys:

            polyId = self.sortedPolyIndices[firstValidPolyIndex]

            if self.polyAreas[polyId] > min_cell_area:

                # we're done
                colapsePoly = False
                # exit the loop

            else:

                # point Ids that span the polygon
                ptIds = self.poly2PointIds[polyId]

                # either the number of points or the number of edges
                # spanning the polygon
                n = len(ptIds)

                # center of the polygon
                barycenter = self.getPolyBarycenter(ptIds)

                # collect the point ids spanning the polygon
                ptIdSet = set()

                # collect the edge to neighbor polygon connectivity
                edge2NeighPoly = {}

                # iterate over edges
                for i in range(n):

                    # the two point Ids spanning the edge
                    ptId1, ptId2 = ptIds[i], ptIds[(i + 1) % n]
                    ptIdSet.add(ptId1)
                    ptIdSet.add(ptId2)

                    # get the other cell that shares this edge
                    self.polydata.GetCellEdgeNeighbors(polyId, ptId1, ptId2, cellIds)

                    if cellIds.GetNumberOfIds() == 2:
                        cellId1 = cellIds.GetId(0)
                        cellId2 = cellIds.GetId(1)
                        if cellId1 == polyId:
                            edge2NeighPoly[(ptId1, ptId2)] = cellId2
                        else:
                            edge2NeighPoly[(ptId1, ptId2)] = cellId1

                # is this an internal polygon whose edges are not on the boundary?
                # (boundary polygons need different treatment, TO DO)
                if len(edge2NeighPoly) == n:

                    # internal polygon

                    for edge, cellId in edge2NeighPoly.items():
                        # the two points spanning the edge
                        p1, p2 = self.points.GetPoint(edge[0]), self.points.GetPoint(edge[1])

                        # area between edge and barycenter
                        a12b = self.getTriangleArea(p1, p2, barycenter)

                        # increase the area of the neighbor polygon
                        self.polyAreas[cellId] += a12b
                        
                    # move the poly boundary points to the barycenter
                    for ptId in ptIds:
                        self.points.SetPoint(ptId, barycenter)

                    # average nodal field to the barycenter position
                    self.averagePointData(ptIds)

                    # set polygon area to zero
                    self.polyAreas[polyId] = 0.

                    # tag this poly for removal
                    self.polydata.Delete(polyId)

                else:

                    # boundary polygon
                    # TO DO 
                    pass

            # look at the next bigger poly
            firstValidPolyIndex += 1

        # now remove the polys tagged for deletion
        self.polydata.RemoveDeletedCells()
        self.polydata.BuildLinks()

    def getPolyBarycenter(self, ptIds):
        """
        Get the center of cloud of points
        @param ptIds list of point Ids
        @return center position
        """
        barycenter = numpy.zeros((3,), numpy.float64)
        point = numpy.zeros((3,), numpy.float64)
        for ptId in ptIds:
            self.polydata.GetPoint(ptId, point)
            barycenter += point
        barycenter /= float(len(ptIds))
        return barycenter


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

    def getTriangleArea(p0, p1, p2):
        """
        Compute the area of the triangle
        @param p0 first vertex
        @param p1 second vertex
        @param p2 third vertex
        @note vertices are in counterclockwise direction
        """
        area = numpy.cross(p1 - p0, p2 - p0)
        return numpy.linalg.norm(area)


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
