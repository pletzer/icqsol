#!/usr/bin/env python

from __future__ import print_function
import vtk
import numpy

def getArrayIndexFromName(data, name):
    """
    Get the array index from its name
    @param data either a vtkCellData or a vtkPointData object
    @param name name
    @return index >= 0 if name exists or -1 if it does not
    """
    numArrays = data.GetNumberOfArrays()
    index = -1
    for i in range(numArrays):
        arr = data.GetArray(i)
        if arr.GetName() == name:
            index = i
            break
    return index

def getArrayIndexFromNameAndProjectOntoCells(pdata, name):
    """
    Get array index from its name, if the array 
    is a point data type then project onto cells
    @param pdata polydata instance
    @param name name
    @return index >= 0 if name exists or -1 if it does not
    """
    cellData = pdata.GetCellData()
    pointData = pdata.GetPointData()
    index = getArrayIndexFromName(cellData, name)
    if index < 0:
        # Maybe a point array?
        index2 = getArrayIndexFromName(pointData, name)

        if index2 >= 0:
            # Project from points to cells
            pointArr = pointData.GetArray(index2)
            numCells = pdata.GetNumberOfPolys()
            cellArr = vtk.vtkDoubleArray()
            cellArr.SetName(name) # same name as the point array
            cellArr.SetNumberOfComponents(1)
            cellArr.SetNumberOfTuples(numCells)
            cells = pdata.GetPolys()
            ptIds = vtk.vtkIdList()
            cells.InitTraversal()
            for i in range(numCells):
                cells.GetNextCell(ptIds)
                # We know for sure that all cells are triangles
                ia, ib, ic = ptIds.GetId(0), \
                    ptIds.GetId(1), ptIds.GetId(2)
                va = pointArr.GetComponent(ia, 0)
                vb = pointArr.GetComponent(ib, 0)
                vc = pointArr.GetComponent(ic, 0)
                cellArr.SetComponent(i, 0, (va + vb + vc)/3.)
            # Add the cell array
            cellData.AddArray(cellArr)
            return getPotentialArrayIndexFromName(name)

    return index


