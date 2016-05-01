#!/usr/bin/env python

from operator import itemgetter
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

        # required in order to get the cell Ids sharing an edge
        self.polydata.BuildLinks()

        # will need to be able to interpolate the nodal data to the
        # new vertices
        self.pointData = self.polydata.GetPointData()
        self.numPointData = self.pointData.GetNumberOfArrays()

        self.computePolyAreas()

        # required so we can get the connectivity between points and 
        # cells
        self.polydata.BuildLinks()

    def computePolyAreas(self):
        """
        Compute the polygon areas
        """
        # polygons
        polys = self.polydata.GetPolys()

        # number of polygons
        numPolys = polys.GetNumberOfCells()

        # polygon areas -- polygons will be sorted according to 
        # their polygon areas
        self.polyAreas = numpy.zeros((numPolys,), numpy.float64)
        polys.InitTraversal()
        ptIds = vtk.vtkIdList()
        for polyId in range(numPolys):
            polys.GetNextCell(ptIds)
            self.polyAreas[polyId] = self.getPolygonArea(ptIds)

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
        points = self.polydata.GetPoints()
        area = numpy.zeros((3,), numpy.float64)
        p0 = numpy.array(points.GetPoint(ptIds.GetId(0)))
        numPts = ptIds.GetNumberOfIds()
        for i in range(1, numPts - 1):
            p1 = numpy.array(points.GetPoint(ptIds.GetId(i    )))
            p2 = numpy.array(points.GetPoint(ptIds.GetId(i + 1)))
            area += numpy.cross(p1 - p0, p2 - p0)
        return numpy.linalg.norm(area)

    def collapsePolygon(self, polyId, zeroPolyList):
        """
        Collapse the vertices of the cell and adjust the 
        neighboring polygons' area
        @param polyId Id of the polygon
        @param zeroPolyList list of polygons with zero area will be updated
        """
        ptIds = vtk.vtkIdList()
        neighPtIds = vtk.vtkIdList()
        neighPolyIds = vtk.vtkIdList()
        points = self.polydata.GetPoints()

        # points of the cell
        self.polydata.GetCellPoints(polyId, ptIds)

        # iterate over the points to compute the barycenter
        barycenter = numpy.zeros((3,), numpy.float64)
        numPts = ptIds.GetNumberOfIds()
        for i in range(numPts):
            barycenter += points.GetPoint(ptIds.GetId(i))
        barycenter /= float(numPts)

        # point under consideration
        pId = vtk.vtkIdList()
        pId.SetNumberOfIds(1)

        # move each vertex of polyId to the barycenter position
        for i in range(numPts):
            # Id of this point
            pI = ptIds.GetId(i)
            self.polydata.GetPoints().SetPoint(pI, barycenter)

        zeroAreas = set()
        for i in range(numPts):

            # get a list of the polys that have this vertex and 
            # correct their area
            pI = ptIds.GetId(i)
            pId.SetId(0, pI)
            self.polydata.GetCellNeighbors(polyId, pId, neighPolyIds)
            numNeigh = neighPolyIds.GetNumberOfIds()
            for j in range(numNeigh):
                neighPolyId = neighPolyIds.GetId(j)
                self.polydata.GetCellPoints(neighPolyId, neighPtIds)
                # correct the area
                area = self.getPolygonArea(neighPtIds)
                self.polyAreas[neighPolyId] = area
                if area < self.EPS :
                    zeroAreas.add(neighPolyId)

        print '=' * 80
        if len(zeroAreas) != 3:
            # something is not right
            print '*** number of zero polys = ', len(zeroAreas), ' != 3'
            return False

        #for neighPolyId in zeroAreas:
        #    print '*** adding ', neighPolyId, ' to the list of zero polys'
        #    zeroPolyList.append(neighPolyId)

        # this poly has now zero area
        self.polyAreas[polyId] = 0.

        # reset the nodal field values to account for the vertices 
        # having moved
        self.averagePointData(ptIds)

        return True

    def coarsen(self, min_cell_area = 1.e-10):
        """
        Coarsen surface by collapsing small polygons
        @param min_cell_area cell area tolerance
        @note operation is in place
        """

        polyId, polyArea = self.findSmallestPolygon()
        if polyId < 0:
            return

        numPolys = self.polydata.GetPolys().GetNumberOfCells()
        count = -1
        success = True
        while polyArea < min_cell_area and polyId > 0 \
                and count < 10*numPolys:

            count += 1
            print '*** count = ', count

            zeroPolyList = []
            
            # collapse polygon. WILL NEED TO DO SOMETHING 
            # ABOUT VERTICES THAT ARE AT THE BOUNDARY
            success = self.collapsePolygon(polyId, zeroPolyList)
            zeroPolyList.append(polyId)

            # find the polygon with the smallest but non-zero area
            polyId, polyArea = self.findSmallestPolygon()  
            numPolys = self.polydata.GetPolys().GetNumberOfCells()

        # delete the zero polys
        for polyId in range(numPolys):
            if self.polyAreas[polyId] < self.EPS:
                self.polydata.DeleteCell(polyId)
            
        self.polydata.RemoveDeletedCells()
        self.polydata.BuildLinks() # not sure if this is required
        self.computePolyAreas()

    def findSmallestPolygon(self):
        """
        Find the smallest polygon whose area is > 0
        @return polygon Id, area
        """

        # collect the indices where area > 0
        inds = numpy.where(self.polyAreas > self.EPS)

        # make sure there is at least one non-degenerate cell
        if len(inds) == 0:
            return -1, -1

        # find the smallest non-zero cells
        polyId, polyArea = min(enumerate(self.polyAreas[inds]), key=itemgetter(1))

        return polyId, polyArea
 
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
