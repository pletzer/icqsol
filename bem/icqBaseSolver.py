#!/usr/bin/env python

from __future__ import print_function
import vtk
import numpy
from icqsol.shapes.icqRefineSurface import RefineSurface
from icqsol.util.icqDataFetcher import getArrayIndexFromNameAndProjectOntoCells

class BaseSolver:

    def __init__(self, pdata, max_edge_length, order=5):
        """
        Constructor
        @param pdata instance of vtkPolyData
        @param max_edge_length maximum edge length, used to turn
                               polygons into triangles
        """

        # triangulate
        rs = RefineSurface(pdata)
        rs.refine(max_edge_length=max_edge_length)
        self.pdata = rs.getVtkPolyData()

        # store the point indices for each cell
        self.ptIdList = []
        ptIds = vtk.vtkIdList()
        polys = self.pdata.GetPolys()
        polys.InitTraversal()
        for i in range(polys.GetNumberOfCells()):
            polys.GetNextCell(ptIds)
            assert(ptIds.GetNumberOfIds() == 3)
            self.ptIdList.append([ptIds.GetId(0),
                                  ptIds.GetId(1),
                                  ptIds.GetId(2)])

        self.points = self.pdata.GetPoints()
        self.polys = self.pdata.GetPolys()
        self.numTriangles = self.polys.GetNumberOfCells()

        # order of the integration, method dependent
        self.order = order

        # set in the derived classes
        self.responseName = 'NO-SET'
        self.sourceName = 'NOT-SET'

    def getVtkPolyData(self):
        """
        Get the (modified) vtkPolyData object
        @return object
        """
        return self.pdata

    def getPoints(self):
        """
        Get grid points
        @return points
        """
        points = self.pdata.GetPoints()
        numPoints = points.GetNumberOfPoints()
        res = numpy.zeros((numPoints, 3), numpy.float64)
        for i in range(numPoints):
            res[i, :] = points.GetPoint(i)
        return res

    def getCells(self):
        """
        Get cell connectivity
        @return cells
        """
        polys = self.pdata.GetPolys()
        numCells = polys.GetNumberOfCells()
        res = numpy.zeros((numCells, 3), numpy.int)
        polys.InitTraversal()
        ptIds = vtk.vtkIdList()
        for i in range(numCells):
            polys.GetNextCell(ptIds)
            res[i, :] = ptIds.GetId(0), ptIds.GetId(1), ptIds.GetId(2)
        return res
    
    def setResponseFieldName(self, name):
        """
        Set the name of the response field
        @param name name
        """
        self.responseName = name

    def setSourceFieldName(self, name):
        """
        Set the name of the source field
        @param name name
        """
        self.sourceName = name

    def getSourceArrayIndex(self):
        """
        Get the source field index, projecting onto triangles is need be
        """
        srcIndex = getArrayIndexFromNameAndProjectOntoCells(self.pdata, self.sourceName)
        if srcIndex < 0:
            msg = 'ERROR: could not find any cell field named {0}!'.format(self.sourceName)
            raise RuntimeError(msg)
        return srcIndex

    def getSourceArray(self, srcIndex):
        """
        Set the source array 
        @param srcIndex source array index
        """
        srcArray = self.pdata.GetCellData().GetArray(srcIndex)
        n = self.pdata.GetNumberOfPolys()
        src = numpy.zeros((n,), numpy.float64)
        for i in range(n):
            src[i] = srcArray.GetComponent(i, 0)
        return src

    def addResponseField(self, rsp):
        """
        Add the response field to the polydata
        @param rsp response numpy array
        """
        rspData = vtk.vtkDoubleArray()
        rspData.SetNumberOfComponents(1)
        n = len(rsp)
        rspData.SetNumberOfTuples(n)
        rspData.SetName(self.responseName)
        for i in range(n):
            rspData.SetTuple(i, [rsp[i]])
        self.pdata.GetCellData().AddArray(rspData)

    def setSourceFromExpression(self, expression):
        """
        Set the source from expression
        @param expression expression of x, y, and z
        """
        from math import sqrt, pi, sin, cos, tan, log, exp
        
        n = self.pdata.GetNumberOfPolys()
        sourceData = vtk.vtkDoubleArray()
        sourceData.SetNumberOfComponents(1)
        sourceData.SetNumberOfTuples(n)
        sourceData.SetName(self.sourceName)
        midPoint = numpy.zeros((3,), numpy.float64)
        ptIds = vtk.vtkIdList()
        cells = self.pdata.GetPolys()
        cells.InitTraversal()
        for i in range(n):
            cell = cells.GetNextCell(ptIds)
            npts = ptIds.GetNumberOfIds()
            midPoint *= 0 # reset
            for j in range(npts):
                midPoint += self.points.GetPoint(ptIds.GetId(j))
            midPoint /= float(npts)
            x, y, z = midPoint
            v = eval(expression)
            sourceData.SetTuple(i, [v])
        self.pdata.GetCellData().AddArray(sourceData)

