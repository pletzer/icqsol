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
        for i in range(numPolys):
            polys.GetCell(i, ptIds)
            self.polydata.InsertNextCell(vtk.VTK_POLYGON, ptIds)

    def refine(self, max_edge_length):
        """
        Refine each cell by adding points on edges longer than max_edge_length
        @param max_edge_length maximum edge length
        @return refined vtkPolyData instance
        """

        # edge to point Ids map
        edge2PointIds = {}

        cells = []

        # iterate over the polygons
        polys = self.polydata.GetPolys()
        numPolys = self.polydata.GetNumberOfPolys()
        ptIds = vtk.vtkIdList()
        polys.InitTraversal()
        for i in range(numPolys):

            polys.GetNextCell(ptIds)

            # number of points spanning the polygon
            numPts = ptIds.GetNumberOfIds()

            # must have at least three points
            if numPts < 3:
                continue

            # compute the normal vector of the polygon plane. Note:
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
            for i in range(numPts):

                # start/end point indices of the edge
                ptIdBeg = ptIds.GetId(i)
                ptIdEnd = ptIds.GetId((i + 1) % numPts)
                edge = [ptIdBeg, ptIdEnd]
                edge.sort()
                ptIdBeg = edge[0]
                ptIdEnd = edge[1]

                e = tuple(edge)
                if e not in edge2PointIds:
                    # new edge
                    ptBeg = numpy.array(self.points.GetPoint(ptIdBeg))
                    ptEnd = numpy.array(self.points.GetPoint(ptIdEnd))
                    d = ptEnd - ptBeg
                    edgeLength = numpy.sqrt(numpy.dot(d, d))
                    # add new points to the edge
                    numSegs = int(math.ceil(edgeLength / max_edge_length))
                    pis = [ptIdBeg]
                    u, v = numpy.dot(ptBeg, uVec), numpy.dot(ptBeg, vVec)
                    uvs[ptIdBeg] = (u, v)
                    for k in range(1, numSegs):
                        pt = ptBeg + k * d / numSegs
                        self.points.InsertNextPoint(pt)
                        ptId = self.points.GetNumberOfPoints() - 1
                        pis.append(ptId)
                        u, v = numpy.dot(pt, uVec), numpy.dot(pt, vVec)
                        uvs[ptId] = (u, v)
                    edge2PointIds[e] = pis

            # triangulate cell and append to list
            cells += self.triangulate(uvs)

        # build the output vtkPolyData object
        ptIds = vtk.vtkIdList()
        numCells = len(cells)
        for i in range(numCells):
            cell = cells[i]
            numPts = len(cell)
            ptIds.SetNumberOfIds(numPts)
            for j in range(numPts):
                ptIds.SetId(j, cell[j])
            self.polydata.InsertNextCell(vtk.VTK_POLYGON, ptIds)

        return self.polydata

    def triangulate(self, uvs):
        """
        Triangulate a set of u, v points
        @param uvs point index to u, v plane coordinates dictionary
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
    pdata2 = rs.refine(max_edge_length=1.5)
    assert(pdata2.GetNumberOfPolys() == 1)


def testAddingThreePoints():

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
    pdata2 = rs.refine(max_edge_length=1.1)
    assert(pdata2.GetNumberOfPolys() == 4)


if __name__ == '__main__':
    testNoRefinement()
    testAddingThreePoints()
