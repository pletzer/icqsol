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
        self.polydataOutput = vtk.vtkPolyData()
        self.pointsOutput = vtk.vtkPoints()
        
        # populate with input data
        points = pdata.GetPoints()
        numPoints = points.GetNumberOfPoints()
        for i in range(numPoints):
            self.pointsOutput.InsertNextPoint(points.GetPoint(i))
        self.polydataOutput.SetPoints(self.pointsOutput)
        
        polys = pdata.GetPolys()
        numPolys = polys.GetNumberOfCells()
        self.polydataOutput.Allocate(numPolys, 1)
        ptIds = vtk.vtkIdList()
        for i in range(numPolys):
            polys.GetCell(i, ptIds)
            self.polydataOutput.InsertNextCell(vtk.VTK_POLYGON, ptIds)

    def refine(self, maxEdgeLength):
        """
        Refine each cell by adding points in the middle of each edge
        @param maxEdgeLength maximum edge length
        @return refined vtkPolyData instance
        """

        # edge to point Ids map
        edge2PointIds = {}

        cells = []
        
        # iterate over the polygons
        polys = self.polydataOutput.GetPolys()
        numPolys = self.polydataOutput.GetNumberOfPolys()
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
            # assuming the polygon is planar
            p0 = numpy.array(self.pointsOutput.GetPoint(ptIds.GetId(0)))
            p1 = numpy.array(self.pointsOutput.GetPoint(ptIds.GetId(1)))
            p2 = numpy.array(self.pointsOutput.GetPoint(ptIds.GetId(2)))
            p1 -= p0
            p2 -= p0
            normal = numpy.cross(p1, p2)
            
            # compute the 2D basis vectors (uVc and vVec) on the plane
            uVec = p1
            vVec = numpy.cross(uVec, normal)
            
            # ordered list of point index to uv coordinates on the plane
            uvs = []
            
            # iterate over edges, add middle edge points if the edge's
            # length is > maxEdgeLength
            for i in range(numPts):
                
                # start/end point indices of the edge
                ptIdBeg = ptIds.GetId(i)
                ptIdEnd = ptIds.GetId((i + 1) % numPts)
                edge = [ptIdBeg, ptIdEnd]
                edge.sort()
                ptIdBeg = edge[0]
                ptIdEnd = edge[1]
                
                e = tuple(edge)
                if not e in edge2PointIds:
                    # new edge
                    ptBeg = numpy.array(self.pointsOutput.GetPoint(ptIdBeg))
                    ptEnd = numpy.array(self.pointsOutput.GetPoint(ptIdEnd))
                    d = ptEnd - ptBeg
                    edgeLength = numpy.sqrt(numpy.dot(d, d))
                    # add new points to the edge
                    numSegs = int(math.ceil(edgeLength / maxEdgeLength))
                    pis = [ptIdBeg]
                    uvs.append([numpy.dot(ptBeg, uVec), numpy.dot(ptBeg, vVec)])
                    for k in range(1, numSegs - 1):
                        pt = ptBeg + k * d / numSegs
                        self.pointsOutput.InsertNextPoint(pt)
                        ptId = self.pointsOutput.GetNumberOfPoints() - 1
                        pis.append(ptId)
                        u, v = numpy.dot(pt, uVec), numpy.dot(pt, vVec)
                        uvs.append((ptId, [u, v]))
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
            self.polydataOutput.InsertNextCell(i, ptIds)
                
        return self.polydataOutput
    
    def triangulate(self, uvs):
        """
        Triangulate set of u, v points
        @param point idex to uvs plane coordinates dictionary
        @return cells 
        """
        points = vtk.vtkPoints()
        pdata = vtk.vtkPolyData()
        delaunay = vtk.vtkDelaunay2D()
        pt = numpy.zeros((3,), numpy.float64)
        for iuv in uvs:
            pt[0:2] = iuv[1]
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
        ptIds = vtk.vtkIdList()
        for i in range(numCells):
            ugrid.getCell(i, ptIds)
            cell = []
            for j in range(ptIds.getNumberOfIds()):
                localId = ptIds.GetId(j)
                ptId = uvs[localId][0]
                cell.append(ptId)
            cells.append(cell)
        return cells

##############################################################################


def test1():
    
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
    
    #print pdata

    rs = RefineSurface(pdata)
    pdata2 = rs.refine(1.1)
    assert(pdata2.GetNumberOfPolys() == 1)


if __name__ == '__main__':
    test1()
