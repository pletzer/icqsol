#!/usr/bin/env python

from __future__ import print_function
import math
import numpy
import vtk
import triangle


class RefineSurface:

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

    def refine(self, max_edge_length):
        """
        Refine each cell by adding points on edges longer than max_edge_length
        @param max_edge_length maximum edge length (> 0)
        @note operation is in place
        """

        # no operation if max_edge_length is not strictly positive
        if max_edge_length <= 0:
            return

        # iterate over the polys
        polys = self.polydata.GetPolys()
        ptIds = vtk.vtkIdList()
        polys.InitTraversal()
        cells = []
        edge2PtIds = {}
        newCellBegEndInds = []
        iNewBeg = 0
        for iPoly in range(polys.GetNumberOfCells()):

            polys.GetNextCell(ptIds)
            numPolyPts = ptIds.GetNumberOfIds()
            polyPtIds = []

            if numPolyPts < 3:
                # need at least three points, next polygon
                continue

            # compute the two tangential unit vectors
            uVec, vVec, normal = self.computeUVNormal(ptIds)
            if normal.dot(normal) == 0:
                # zero area polygon, nothing to do
                continue

            # iterate over edges
            for iNode in range(numPolyPts):
                i0 = ptIds.GetId(iNode)
                i1 = ptIds.GetId((iNode + 1) % numPolyPts)
                edge = (i0, i1)
                edgeCompl = (i1, i0)

                edgePtIds = edge2PtIds.get(edgeCompl, [])
                if len(edgePtIds) > 0:
                    # edge has already been split, take edge and reverse order
                    edgePtIds.reverse()
                    edgePtIds = [i0] + edgePtIds[:-1]
                    polyPtIds += edgePtIds
                    continue

                p0 = numpy.array(self.points.GetPoint(i0))
                p1 = numpy.array(self.points.GetPoint(i1))
                dEdge = p1 - p0
                edgeLength = numpy.sqrt(dEdge.dot(dEdge))
                numSegs = int(max(1, math.ceil(edgeLength/max_edge_length)))

                # add points along edge
                d10 = (p1 - p0) / float(numSegs)
                edgePtIds.append(i0)
                for iSeg in range(1, numSegs):
                    pt = p0 + iSeg*d10
                    # index of point to add
                    ptId = self.points.GetNumberOfPoints()
                    # insert point
                    self.points.InsertNextPoint(pt)
                    edgePtIds.append(ptId)
                    # interpolate the field values along the edge
                    w0 = (numSegs - iSeg)/float(numSegs)
                    w1 = iSeg/float(numSegs)
                    for name in self.pointData:
                        v0 = numpy.array(self.pointData[name].GetTuple(i0))
                        v1 = numpy.array(self.pointData[name].GetTuple(i1))
                        interpTuple = w0*v0 + w1*v1
                        self.pointData[name].InsertNextTuple(interpTuple)

                edge2PtIds[edge] = edgePtIds
                polyPtIds += edgePtIds

            # triangulate the cell
            polyCells = self.triangulatePolygon(uVec, vVec, polyPtIds,
                                                max_edge_length)
            cells += polyCells
            iNewEnd = iNewBeg + len(polyCells)
            # all the cell data have the same value inside this polygon
            # store the srta/end cell indices 
            newCellBegEndInds.append( (iNewBeg, iNewEnd) )
            iNewBeg = iNewEnd

        # build the output vtkPolyData object
        newPolyData = vtk.vtkPolyData()

        # set points
        newPolyData.SetPoints(self.points)

        # set all the point data
        for name in self.pointData:
            newPolyData.GetPointData().AddArray(self.pointData[name])

        # build the connectivity
        numPolys = len(cells)
        newPolyData.Allocate(numPolys, 1)
        ptIds = vtk.vtkIdList()
        for cell in cells:
            numPolyPts = len(cell)
            ptIds.SetNumberOfIds(numPolyPts)
            for j in range(numPolyPts):
                ptIds.SetId(j, cell[j])
            newPolyData.InsertNextCell(vtk.VTK_POLYGON, ptIds)

        numTuples = self.polydata.GetPolys().GetNumberOfCells()
        newNumTuples = newPolyData.GetPolys().GetNumberOfCells()

        # add the cell arrays
        cd = self.polydata.GetCellData()
        numCellArrays = cd.GetNumberOfArrays()
        for el in range(numCellArrays):
            arr = cd.GetArray(el)
            name = arr.GetName()
            numComps = arr.GetNumberOfComponents()

            newArr = self.cellData[name]
            newArr.SetNumberOfComponents(numComps)
            newArr.SetNumberOfTuples(newNumTuples)
            for i in range(numTuples):
                vals = arr.GetTuple(i)
                iNewBeg, iNewEnd = newCellBegEndInds[i]
                for j in range(iNewBeg, iNewEnd):
                    newArr.SetTuple(j, vals)
            newPolyData.GetCellData().AddArray(newArr)

        # Reset the polydata struct
        self.polydata = newPolyData

    def computeUVNormal(self, ptIds):
        """
        Compute the two tangential unit vectors and the normal vector
        @param ptIds point indices
        @return u vector, v vector, normal
        """
        uVec = numpy.zeros((3,), numpy.float64)
        vVec = numpy.zeros((3,), numpy.float64)
        normal = numpy.zeros((3,), numpy.float64)
        numPts = ptIds.GetNumberOfIds()
        if numPts < 3:
            return uVec, vVec, normal
        p0 = numpy.array(self.points.GetPoint(ptIds.GetId(0)))
        for i in range(1, numPts - 1):
            dp1 = numpy.array(self.points.GetPoint(ptIds.GetId(i))) - p0
            dp2 = numpy.array(self.points.GetPoint(ptIds.GetId(i + 1))) - p0
            perp = numpy.cross(dp1, dp2)
            pDotp = perp.dot(perp)
            if pDotp > 0:
                normal = perp / numpy.sqrt(pDotp)
                uVec = dp1 / numpy.sqrt(dp1.dot(dp1))
                vVec = numpy.cross(normal, uVec)
                return uVec, vVec, normal
        return uVec, vVec, normal

    def triangulatePolygon(self, uVec, vVec, polyPtIds, max_edge_length):
        """
        Triangulate polygon using the uVec x vVec projection
        @param uVec unit vector tangential to the polygon (input)
        @param vVec second unit vector tangential to the polygon (input)
        @param polyPtIds list of polygon's point indices (input and output)
        @param max_edge_length maximum edge length (input)
        @return list of cells (list of point indices)
        """

        numPolyPts = len(polyPtIds)
        if numPolyPts < 3:
            # need at least three points
            return []

        area = self.getPolygonArea(polyPtIds)
        if area < 1.e-15:
            return []

        self.removeDegeneratePoints(polyPtIds)

        # reset
        numPolyPts = len(polyPtIds)
        if numPolyPts < 3:
            return []

        # reference position
        p0 = numpy.array(self.points.GetPoint(polyPtIds[0]))

        pts = []
        pointData = {}
        for name in self.pointData:
            pointData[name] = []
        # project each point onto the plane
        for i in range(numPolyPts):
            ptId = polyPtIds[i]
            pos = numpy.array(self.points.GetPoint(ptId))
            pos -= p0
            pts.append((pos.dot(uVec), pos.dot(vVec)))
            for name, pd in pointData.items():
                pd.append(tuple(self.pointData[name].GetTuple(ptId)))

        # list of segments
        segs = [(i, (i + 1)%numPolyPts) for i in range(numPolyPts)]

        tri = triangle.Triangle()
        tri.set_points(pts)
        tri.set_segments(segs)
        attrs, names, components = self.pointDataToAttributes(pointData)
        tri.set_attributes(attrs)

        # internal points will be added if triangle area exceeds threshold
        maxArea = None
        if max_edge_length < float('inf') and max_edge_length > 0.:
            maxArea = 0.5 * max_edge_length**2

        # triangulate
        # p: triangulate a straight planar graph
        # z: zero based indexing
        # Q: quiet mode
        tri.triangulate(area=maxArea, mode='pzQ')

        nodes = tri.get_nodes()
        polyCells = tri.get_triangles()
        interpolatedAttrs = tri.get_attributes()

        # add internal vertices
        pIndex2PtId = {}
        for pIndex in range(len(pts), len(nodes)):
            ptId = self.points.GetNumberOfPoints()
            u, v = nodes[pIndex][0]
            p = p0 + u*uVec + v*vVec
            # insert new point
            self.points.InsertNextPoint(p)
            polyPtIds.append(ptId)
            pIndex2PtId[pIndex] = ptId

        # add point data inside the cell
        for name, pd in self.pointData.items():
            numTuples = pd.GetNumberOfTuples()
            newSize = numTuples + len(nodes) - len(pts)
            success = pd.Resize(newSize)
            if success != 1:
                msg = 'Failed to resize vtkDoubleArray associated with field {0}'.format(name)
                raise MemoryError(msg)
            pd.SetNumberOfTuples(newSize)
            pd.SetName(name)

            # add the new internal point data
            iBeg = names.index(name)
            iEnd = len(names) - names[::-1].index(name)
            for pIndex in range(len(pts), len(nodes)):
                ptId = pIndex2PtId[pIndex]
                vals = interpolatedAttrs[pIndex][iBeg:iEnd]
                pd.SetTuple(ptId, vals)

        cells = []
        for c in polyCells:
            ia, ib, ic = c[0]
            cell = [polyPtIds[ia], polyPtIds[ib], polyPtIds[ic]]
            cells.append(cell)

        return cells

    def getPolygonArea(self, polyPtIds):
        """
        Compute the (scalar) area of a polygon
        @param polyPtIds list of point indices
        @return area
        """
        area = numpy.zeros((3,), numpy.float64)
        p0 = numpy.array(self.points.GetPoint(polyPtIds[0]))
        numPolyPts = len(polyPtIds)
        for i in range(1, numPolyPts - 1):
            p1 = numpy.array(self.points.GetPoint(polyPtIds[i]))
            p2 = numpy.array(self.points.GetPoint(polyPtIds[i + 1]))
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
        
    def pointDataToAttributes(self, pointData):
        """
        Convert point data dictionary 
        into an attribute array which can be passed to the triangle program
        @param pointData {name: [(comp0, comp1, ...), ...], ...}
        @return attributes, names, and components as separate lists
        """
        attrs = []
        names = []
        components = []

        if len(pointData) == 0:
            return attrs, names, components

        scalarNames = list(pointData.keys())
        scalarData = list(pointData.values())
        
        # set the name and component lists, which is 
        # the same for each node
        for j in range(len(scalarNames)):
            name = scalarNames[j]
            data = scalarData[j]
            for k in range(len(data[0])):
                names.append(name)
                components.append(k)
                
        numPts = len(scalarData[0])
        for i in range(numPts):
            a = []

            # iterate ove the scalar field names
            for j in range(len(scalarNames)):
                data = scalarData[j]
          
                # iterate over the components
                for k in range(len(data[0])):
                    a.append(data[i][k])

            attrs.append(tuple(a))

        return attrs, names, components

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
    assert(rs.getVtkPolyData().GetNumberOfPolys() == 4)


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
    assert(rs.getVtkPolyData().GetNumberOfPolys() == 8)


if __name__ == '__main__':
    testNoRefinement()
    testAddingThreePointsThenMore()
    testStartingWithTwoCells()
