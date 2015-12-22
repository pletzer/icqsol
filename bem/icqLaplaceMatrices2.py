#!/usr/bin/env python

import vtk
import numpy
from icqsol.shapes.icqRefineSurface import RefineSurface
from icqsol.bem.icqPotentialIntegrals import PotentialIntegrals
from icqsol.bem.icqQuadrature import triangleQuadrature

FOUR_PI = 4. * numpy.pi

def getIntegralOneOverROff(xObs, paSrc, pbSrc, pcSrc, order):
    def green(x):
        r = xObs - x
        return 1.0/numpy.sqrt(r.dot(r))
    return triangleQuadrature(order, paSrc, pbSrc, pcSrc, green)

def getIntegralMinusOneOverRCubeOff(xObs, paSrc, pbSrc, pcSrc, normalSrc, order):
    def kreen(x):
        r = xObs - x
        return normalSrc.dot(r)/numpy.sqrt(r.dot(r))**3
    return triangleQuadrature(order, paSrc, pbSrc, pcSrc, kreen)


class LaplaceMatrices2:

    def __init__(self, pdata, max_edge_length, order):
        """
        Constructor
        @param pdata instance of vtkPolyData
        @param max_edge_length maximum edge length, used to turn
                               polygons into triangles
        @param order order of the Gauss quadrature scheme
        """

        assert(order > 0 and order <= 5)
        
        self.normalDerivativeName = 'normal_derivative'
        self.potentialName = 'potential'

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

        shp = (self.numTriangles, self.numTriangles)
        self.gMat = numpy.zeros(shp, numpy.float64)
        self.kMat = numpy.zeros(shp, numpy.float64)
        
        self.order = order
        self.__computeMatrices()

    def setPotentialName(self, name):
        self.potentialName = name

    def setNormalDerivativeName(self, name):
        self.normalDerivativeName = name

    def __computeCoupling(self, iObs, jSrc):
                
        # observer
        cellObs = self.ptIdList[iObs]
        paObs = numpy.array(self.points.GetPoint(cellObs[0]))
        pbObs = numpy.array(self.points.GetPoint(cellObs[1]))
        pcObs = numpy.array(self.points.GetPoint(cellObs[2]))
        xObs = (paObs + pbObs + pcObs) / 3.0

        # source
        cellSrc = self.ptIdList[jSrc]
        paSrc = numpy.array(self.points.GetPoint(cellSrc[0]))
        pbSrc = numpy.array(self.points.GetPoint(cellSrc[1]))
        pcSrc = numpy.array(self.points.GetPoint(cellSrc[2]))
        
        normal = numpy.cross(pbSrc - paSrc, pcSrc - paSrc)
        normal /= numpy.linalg.norm(normal)
        elev = (xObs - paSrc).dot(normal)
        
        if iObs == jSrc:
            
            # singular term
            
            pot0ab = PotentialIntegrals(xObs, paSrc, pbSrc, self.order)
            pot0bc = PotentialIntegrals(xObs, pbSrc, pcSrc, self.order)
            pot0ca = PotentialIntegrals(xObs, pcSrc, paSrc, self.order)
        
            self.gMat[iObs, jSrc] = pot0ab.getIntegralOneOverR(elev) + \
                                    pot0bc.getIntegralOneOverR(elev) + \
                                    pot0ca.getIntegralOneOverR(elev)
            self.gMat[iObs, jSrc] /= FOUR_PI
        
            self.kMat[iObs, jSrc] = 0.0
        
        else:
        
            # off diagonal term 

            self.gMat[iObs, jSrc] = getIntegralOneOverROff(xObs,
                paSrc, pbSrc, pcSrc, self.order)
            self.gMat[iObs, jSrc] /= FOUR_PI
            
            self.kMat[iObs, jSrc] = getIntegralMinusOneOverRCubeOff(xObs,
                paSrc, pbSrc, pcSrc, normal, self.order)
            self.kMat[iObs, jSrc] /= FOUR_PI

    def __computeMatrices(self):

        # iterate over the observer triangles
        for iObs in range(self.numTriangles):

            # iterate over source triangles
            for jSrc in range(0, self.numTriangles):
            #for jSrc in range(0, self.numTriangles):
                
                self.__computeCoupling(iObs, jSrc)
                #self.gMat[jSrc, iObs] = self.gMat[iObs, jSrc]
                #self.kMat[jSrc, iObs] = self.kMat[iObs, jSrc]

    def getGreenMatrix(self):
        """
        Return the Green function matrix
        @return matrix
        """
        return self.gMat

    def getNormalDerivativeGreenMatrix(self):
        """
        Return the normal derivative of the Green function matrix
        @return matrix
        """
        return self.kMat

    def computeNeumannFromDirichlet(self, dirichletExpr):
        """
        Get the Neumann boundary values from the Dirichlet boundary conditions
        @param dirichletExpr expression for the potential values
        @return response
        """
        from math import pi, sin, cos, log, exp, sqrt
        
        kMat = self.getNormalDerivativeGreenMatrix()
        n = kMat.shape[0]
        
        v = numpy.zeros((n,), numpy.float64)
        
        # subtract residue
        for i in range(n):
            kMat[i, i] -= 0.5
        
        # evaluate potential on each triangle center
        ptIds = vtk.vtkIdList()
        polys = self.pdata.GetPolys()
        points = self.pdata.GetPoints()
        polys.InitTraversal()
        for i in range(n):
            
            polys.GetNextCell(ptIds)
            ia, ib, ic = ptIds.GetId(0), ptIds.GetId(1), ptIds.GetId(2)
            pa = numpy.array(points.GetPoint(ia))
            pb = numpy.array(points.GetPoint(ib))
            pc = numpy.array(points.GetPoint(ic))
            pMid = (pa + pb + pc)/3.
            x, y, z = pMid
            
            # set the Dirichlet value
            v[i] = eval(dirichletExpr)
        
        gMat = self.getGreenMatrix()
        
        # solve
        gM1 = numpy.linalg.inv(gMat)
        
        normalDerivative = gM1.dot(kMat).dot(v)
        
        # add field
        cellData = self.pdata.GetCellData()
        
        potentialData = vtk.vtkDoubleArray()
        potentialData.SetNumberOfComponents(1)
        potentialData.SetNumberOfTuples(n)
        potentialData.SetName(self.potentialName)

        normalDerivData = vtk.vtkDoubleArray()
        normalDerivData.SetNumberOfComponents(1)
        normalDerivData.SetNumberOfTuples(n)
        normalDerivData.SetName(self.normalDerivativeName)

        for i in range(n):
            potentialData.SetTuple(i, [v[i],])
            normalDerivData.SetTuple(i, [normalDerivative[i],])

        cellData.AddArray(potentialData)
        cellData.AddArray(normalDerivData)
        
        return normalDerivative 


###############################################################################


def testSingTriangle():

    "Single triangle"

    h = 0.1
    # create set of points
    points = vtk.vtkPoints()
    points.SetNumberOfPoints(3)
    points.SetPoint(0, [1., -1.*h/3., -1.*h/3.])
    points.SetPoint(1, [1., 2.*h/3., -1.*h/3.])
    points.SetPoint(2, [1., -1.*h/3., 2.*h/3.])

    # create vtkPolyData object
    pdata = vtk.vtkPolyData()
    pdata.SetPoints(points)
    ptIds = vtk.vtkIdList()
    ptIds.SetNumberOfIds(3)
    ptIds.SetId(0, 0)
    ptIds.SetId(1, 1)
    ptIds.SetId(2, 2)
    pdata.Allocate(1, 1)
    pdata.InsertNextCell(vtk.VTK_POLYGON, ptIds)

    for order in range(1, 6):
        lslm = LaplaceMatrices2(pdata,
                                max_edge_length=1000.,
                                order=order)
        print 'order = ', order
        print 'g matrix: ', lslm.getGreenMatrix()
        print 'k matrix: ', lslm.getNormalDerivativeGreenMatrix()

def testTwoTrianglesCoplanar():

    "Two triangles"

    # create set of points
    points = vtk.vtkPoints()
    points.SetNumberOfPoints(4)
    points.SetPoint(0, [0., 0., 0.])
    points.SetPoint(1, [1., 0., 0.])
    points.SetPoint(2, [0., 1., 0.])
    points.SetPoint(3, [1., 1., 0.])

    # create vtkPolyData object
    pdata = vtk.vtkPolyData()
    pdata.SetPoints(points)
    
    pdata.Allocate(2, 1)
    ptIds = vtk.vtkIdList()
    ptIds.SetNumberOfIds(3)
    
    ptIds.SetId(0, 0)
    ptIds.SetId(1, 1)
    ptIds.SetId(2, 2)
    pdata.InsertNextCell(vtk.VTK_POLYGON, ptIds)
    ptIds.SetId(0, 1)
    ptIds.SetId(1, 3)
    ptIds.SetId(2, 2)
    pdata.InsertNextCell(vtk.VTK_POLYGON, ptIds)

    for order in range(1, 6):
        lslm = LaplaceMatrices2(pdata,
                                max_edge_length=1000.,
                                order=order)
        print 'order = ', order
        print 'g matrix: ', lslm.getGreenMatrix()
        print 'k matrix: ', lslm.getNormalDerivativeGreenMatrix()

def testTwoTriangles():

    "Two triangles"

    # create set of points
    points = vtk.vtkPoints()
    points.SetNumberOfPoints(4)
    points.SetPoint(0, [0., 0., 0.])
    points.SetPoint(1, [1., 0., 0.])
    points.SetPoint(2, [0., 1., 0.])
    points.SetPoint(3, [0., 0., 1.])

    # create vtkPolyData object
    pdata = vtk.vtkPolyData()
    pdata.SetPoints(points)
    
    pdata.Allocate(2, 1)
    ptIds = vtk.vtkIdList()
    ptIds.SetNumberOfIds(3)
    
    ptIds.SetId(0, 0)
    ptIds.SetId(1, 1)
    ptIds.SetId(2, 3)
    pdata.InsertNextCell(vtk.VTK_POLYGON, ptIds)
    ptIds.SetId(0, 0)
    ptIds.SetId(1, 3)
    ptIds.SetId(2, 2)
    pdata.InsertNextCell(vtk.VTK_POLYGON, ptIds)

    for order in range(1, 6):
        lslm = LaplaceMatrices2(pdata,
                                max_edge_length=1000.,
                                order=order)
        print 'order = ', order
        print 'g matrix: ', lslm.getGreenMatrix()
        print 'k matrix: ', lslm.getNormalDerivativeGreenMatrix()

if __name__ == '__main__':
    #testSingleTriangle()
    #testTwoTrianglesCoplanar()
    testTwoTriangles()
