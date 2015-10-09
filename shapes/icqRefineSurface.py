#!/usr/bin/env python

import math
import numpy
import operator
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
    
    def computeNormal(self, verts):
        """
        Compute the normal to the polygon by taking the first three vertices
        that are not degenerate
        @param vertices as numpy arrays
        @return normal vector or (0., 0., 0.), tangential lengths, and area
        """
        pa = verts[0]
        pb = verts[1] - pa
        for i in range(2, len(verts)):
            pc = verts[i] - pa
            area = numpy.cross(pb, pc)
            areaVal = numpy.sqrt(numpy.dot(area, area))
            if areaVal > 0.0:
                return area/areaVal, pb, pc, areaVal
        return numpy.array([0.,0.,0.]), pb, pc, 0.0

    def refine(self, max_edge_length):
        """
        Refine each cell by adding points on edges longer than max_edge_length
        @param max_edge_length maximum edge length
        @note operation is in place
        """

        # iterate over the polys
        polys = self.polydata.GetPolys()
        ptIds = vtk.vtkIdList()
        polys.InitTraversal()
        cells = []
        edges = set()
        for iPoly in range(polys.GetNumberOfCells()):
            
            polys.GetNextCell(ptIds)
            
            # compute the two tangential unit vectors
            uVec, vVec, normal = self.computeUVNormal(self.points, ptIds)
            if normal.dot(normal) == 0:
                # zero area polygon, nothing to do
                continue
            
            # collect the point ids of the polygon
            polyPtIds = []

            # iterate over edges
            numNodes = ptIds.GetNumberOfIds()
            for iNode in range(numNodes):
                i0 = ptIds.GetId(iNode)
                i1 = ptIds.GetId((iNode + 1) % numNodes)
                edge = [i0, i1]
                edge.sort()
                edge = tuple(edge)
                
                polyPtIds.append(i0)

                if edge in edges:
                    # edge has already been split, go to next edge
                    continue
                
                edges.add(edge)
                
                i0, i1 = edge
                p0 = numpy.array(self.points.GetPoint(edge[0]))
                p1 = numpy.array(self.points.GetPoint(edge[1]))
                dEdge = p1 - p0
                edgeLength = numpy.sqrt(dEdge.dot(dEdge))
                numSegs = max(1, math.ceil(edgeLength/max_edge_length))

                # add points along edge
                pt = p0
                d10 = (p1 - p0) / float(numSegs)
                for iSeg in range(numSegs - 1):
                    pt += d10
                    ptId = self.points.GetNumberOfPoints()
                    # insert point
                    self.points.InsertNextPoint(pt)
                    polyPtIds.append(ptId)

            # triangulate the cell
            polyCells = self.triangulatePolygon(uVec, vVec, self.points, polyPtIds)
            cells += polyCells

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

        # Reset the polydata struct
        self.polydata = pdata
            
    def computeUVNormal(self, points, ptIds):
        """
        Compute the two tangential unit vectors and the normal vector
        @param points vtkPoints instance
        @param ptIds point indices
        @return u vector, v vector, normal
        """
        uVec = numpy.zeros((3,), numpy.float64)
        vVec = numpy.zeros((3,), numpy.float64)
        numPts = ptIds.GetNumberOfIds()
        if numPts < 3:
            return uVec, vVec, uVec
        p0 = numpy.array(points.GetPoint(ptIds.GetId(0)))
        for i in range(1, numPts - 1):
            dp1 = numpy.array(points.GetPoint(ptIds.GetId(i))) - p0
            dp2 = numpy.array(points.GetPoint(ptIds.GetId(i + 1))) - p0
            perp = numpy.cross(dp1, dp2)
            pDotp = perp.dot(perp)
            if pDotp > 0:
                normal = perp / numpy.sqrt(pDotp)
                uVec = dp1 / numpy.sqrt(dp1.dot(dp1))
                vVec = numpy.cross(normal, uVec)
                return uVec, vVec, normal
        return uVec, vVec, numpy.zeros((3,), numpy.float64)

    def triangulatePolygon(self, uVec, vVec, points, polyPtIds):
        """
        Triangulate polygon using the uVec x vVec projection
        @param uVec unit vector tangential to the polygon
        @param vVec second unit vector tangential to the polygon
        @param points vtkPoints object
        @param polyPtIds list of polygon's point indices
        @return list of cells (list of point indices)
        """
        
        numPts = len(polyPtIds)
        if numPts < 3:
            return []
        elif numPts == 3:
            # no need to do any triangulation, just return the cell
            return [polyPtIds]
        
        pts = vtk.vtkPoints()
        pts.SetNumberOfPoints(numPts)
        
        # project each point onto the plane
        pt = numpy.zeros((3,), numpy.float64)
        for i in range(numPts):
            p = numpy.array(points.GetPoint(polyPtIds[i]))
            pt[0:2] = p.dot(uVec), p.dot(vVec)
            pts.SetPoint(i, pt)
        
        pdata = vtk.vtkPolyData()
        pdata.SetPoints(pts)
        
        delaunay = vtk.vtkDelaunay2D()
        delaunay.SetTolerance(1.e-5)
        delaunay.SetAlpha(0.0)
        if vtk.VTK_MAJOR_VERSION >= 6:
            delaunay.SetInputData(pdata)
        else:
            delaunay.SetInput(pdata)
        
        delaunay.Update()
        
        ugrid = delaunay.GetOutput()
        cells = []
        numCells = ugrid.GetNumberOfCells()
        if numCells == 0:
            # Delaunay triangulation failed
            # not sure why this is happening...
            return [polyPtIds]
        else:
            for i in range(numCells):
                ptIds = ugrid.GetCell(i).GetPointIds()
                cell = []
                for j in range(ptIds.GetNumberOfIds()):
                    localId = ptIds.GetId(j)
                    ptId = polyPtIds[localId]
                    cell.append(ptId)
                cells.append(cell)
        return cells
    
    
    def refine2(self, max_edge_length):
        """
        Refine each cell by adding points on edges longer than max_edge_length
        @param max_edge_length maximum edge length
        @note operation is in place
        """
        
        print '***** max_edge_length = ', max_edge_length
        print '***** input: number of polys = ', self.polydata.GetNumberOfPolys()

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
                print '+++++ numPts = ', numPts
                continue
            
            verts = [ numpy.array(self.points.GetPoint(ptIds.GetId(i))) \
                     for i in range(numPts)]
            normal, p1, p2, area = self.computeNormal(verts)
            if area == 0:
                print '????? area = ', area
                # degenerate cell
                continue

            # compute the 2D basis vectors (uVec and vVec) on the plane
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
                    
                    # point indices
                    pis = [ptIdBeg]
                    
                    # increment along edge
                    delta = d/float(numSegs)
                    
                    for k in range(1, numSegs):
                        pt = ptBeg + k*delta
                        
                        # insert new point
                        print '**** inserting point ', pt, ' along edge ', e, ' in cell ', ptIds
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
                
                # add middle point for fat triangles
                #if area >  * numpy.sqrt(p1.dot(p1)) * numpy.sqrt(p2.dot(p2)):
                #    print '*** adding middle point'
                #    pMid = reduce(operator.add, verts) / float(numPts)
                #    self.points.InsertNextPoint(pMid)
                #    ptId = self.points.GetNumberOfPoints() - 1
                #    uvs[ptId] = (u, v)

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
        
        print '.....output: number of polys: ', pdata.GetNumberOfPolys()
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
