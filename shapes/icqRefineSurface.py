#!/usr/bin/env python

import math
import numpy
import vtk


class RefineSurface:

    def __init__(self, pdata):
        """
        Constructor
        @param pdata vtkPolyData instance
        """
        self.polydata = vtk.vtkPolyData()
        self.points = vtk.vtkPoints()

        # populate with input data
        points = pdata.GetPoints()
        numPoints = points.GetNumberOfPoints()
        for i in range(numPoints):
            self.points.InsertNextPoint(points.GetPoint(i))
        self.polydata.SetPoints(self.points)

        polys = pdata.GetPolys()
        numPolys = polys.GetNumberOfCells()
        self.polydata.Allocate(numPolys, 1)
        ptIds = vtk.vtkIdList()
        polys.InitTraversal()
        for i in range(numPolys):
            polys.GetNextCell(ptIds)
            self.polydata.InsertNextCell(vtk.VTK_POLYGON, ptIds)

    def getVtkPolyData(self):
        """
        Get the vtkPolyData instance associated with the grid
        @return vtkPolyData instance
        """
        return self.polydata

    def refine(self, max_edge_length):
        """
        Refine each cell by adding points on edges longer than max_edge_length
        @param max_edge_length maximum edge length
        @note operation is in place
        """

        # edge to point Ids map
        edge2PointIds = {}

        cells = []

        # iterate over the polygons
        polys = self.polydata.GetPolys()
        numPolys = polys.GetNumberOfCells()
        ptIds = vtk.vtkIdList()
        polys.InitTraversal()
        for iPoly in range(numPolys):

            polys.GetNextCell(ptIds)

            # number of points spanning the polygon
            numPts = ptIds.GetNumberOfIds()

            # must have at least three points
            if numPts < 3:
                continue

            # compute the normal vector on the polygon's plane. Note:
            # assuming the polygon is planar -- only taking the first
            # 3 points
            p0 = numpy.array(self.points.GetPoint(ptIds.GetId(0)))
            p1 = numpy.array(self.points.GetPoint(ptIds.GetId(1)))
            p2 = numpy.array(self.points.GetPoint(ptIds.GetId(2)))
            p1 -= p0
            p2 -= p0
            normal = numpy.cross(p1, p2)
            normal /= math.sqrt(numpy.dot(normal, normal))

            # compute the 2D basis vectors (uVc and vVec) on the plane
            uVec = p1 / math.sqrt(numpy.dot(p1, p1))
            vVec = numpy.cross(normal, uVec)

            # ordered list of point index to uv coordinates on the plane
            uvs = {}

            # iterate over edges, add middle edge points if the edge's
            # length is > max_edge_length
            for iPoint in range(numPts):

                # start/end point indices of the edge
                ptIdBeg = ptIds.GetId(iPoint)
                ptIdEnd = ptIds.GetId((iPoint + 1) % numPts)
                edge = [ptIdBeg, ptIdEnd]
                edge.sort()
                ptIdBeg, ptIdEnd = edge
                e = tuple(edge)

                if e not in edge2PointIds:

                    # new edge
                    ptBeg = numpy.array(self.points.GetPoint(ptIdBeg))
                    ptEnd = numpy.array(self.points.GetPoint(ptIdEnd))
                    d = ptEnd - ptBeg
                    edgeLength = numpy.sqrt(numpy.dot(d, d))

                    # add new points to the edge
                    numSegs = max(1, int(math.ceil(edgeLength/max_edge_length)))
                    pis = [ptIdBeg]
                    delta = d/float(numSegs)
                    for k in range(1, numSegs):
                        pt = ptBeg + k*delta
                        self.points.InsertNextPoint(pt)
                        ptId = self.points.GetNumberOfPoints() - 1
                        pis.append(ptId)
                    pis.append(ptIdEnd)

                    edge2PointIds[e] = pis

                # compute the location of the edge points in u, v plane coordinates
                for ptId in edge2PointIds[e]:
                    pt = self.points.GetPoint(ptId)
                    u, v = numpy.dot(pt, uVec), numpy.dot(pt, vVec)
                    uvs[ptId] = (u, v)

            # triangulate cell and append to list
            cells += self.triangulate(uvs)

        # build the output vtkPolyData object
        pdata = vtk.vtkPolyData()
        ptIds = vtk.vtkIdList()
        pdata.SetPoints(self.points)
        numPolys = len(cells)
        pdata.Allocate(numPolys, 1)
        for cell in cells:
            numPts = len(cell)
            ptIds.SetNumberOfIds(numPts)
            for j in range(numPts):
                ptIds.SetId(j, cell[j])
            pdata.InsertNextCell(vtk.VTK_POLYGON, ptIds)

        self.polydata = pdata

    def triangulate(self, uvs):
        """
        Triangulate a set of u, v points
        @param uvs dictionary containing point index to u, v plane coordinates map
        @return cells
        """
        uvIds = uvs.keys()
        uvPoints = uvs.values()
        points = vtk.vtkPoints()
        pdata = vtk.vtkPolyData()
        delaunay = vtk.vtkDelaunay2D()
        pt = numpy.zeros((3,), numpy.float64)
        for uv in uvPoints:
            pt[0:2] = uv
            points.InsertNextPoint(pt)
        pdata.SetPoints(points)
        if vtk.VTK_MAJOR_VERSION >= 6:
            delaunay.SetInputData(pdata)
        else:
            delaunay.SetInput(pdata)

        delaunay.Update()

        ugrid = delaunay.GetOutput()
        cells = []
        numCells = ugrid.GetNumberOfCells()
        for i in range(numCells):
            ptIds = ugrid.GetCell(i).GetPointIds()
            cell = []
            for j in range(ptIds.GetNumberOfIds()):
                localId = ptIds.GetId(j)
                ptId = uvIds[localId]
                cell.append(ptId)
            cells.append(cell)
        return cells

##############################################################################


def printVtkPolyData(pdata):

    points = pdata.GetPoints()
    polys = pdata.GetPolys()
    ptIds = vtk.vtkIdList()
    numPolys = polys.GetNumberOfCells()
    print 'Number of polygons: {}'.format(numPolys)
    polys.InitTraversal()
    for i in range(numPolys):
        cell = polys.GetNextCell(ptIds)
        numPts = ptIds.GetNumberOfIds()
        print '\tCell {} has {} points: '.format(i, numPts)
        for j in range(numPts):
            ptId = ptIds.GetId(j)
            pt = points.GetPoint(ptId)
            print '\t\t{} -> {}'.format(ptId, pt)


def testNoRefinement():

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

    rs = RefineSurface(pdata)
    rs.refine(max_edge_length=1.5)
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

    rs = RefineSurface(pdata)
    rs.refine(max_edge_length=1.1)
    #printVtkPolyData(rs.getVtkPolyData())
    assert(rs.getVtkPolyData().GetNumberOfPolys() == 4)
    #rs.refine(max_edge_length=0.1)
    #assert(rs.getVtkPolyData().GetNumberOfPolys() == 104)


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

    rs = RefineSurface(pdata)
    rs.refine(max_edge_length=1.1)
    pdata2 = rs.getVtkPolyData()
    #printVtkPolyData(pdata2)
    assert(rs.getVtkPolyData().GetNumberOfPolys() == 8)
    #rs.refine(max_edge_length=0.1)
    #assert(rs.getVtkPolyData().GetNumberOfPolys() == 208)


if __name__ == '__main__':
    testNoRefinement()
    testAddingThreePointsThenMore()
    testStartingWithTwoCells()
