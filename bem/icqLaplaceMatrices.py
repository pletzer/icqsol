#!/usr/bin/env python

import vtk
import numpy
from icqsol.shapes.icqRefineSurface import RefineSurface
from icqsol.bem.icqPotentialIntegrals import PotentialIntegrals
from icqsol.bem.icqQuadrature import gaussPtsAndWeights

FOUR_PI = 4. * numpy.pi

class LaplaceMatrices:

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
            
    def getVtkPolyData(self):
        """
        Get the (modified) vtkPolyData object
        @return object
        """
        return self.pdata

    def __computeMatrices(self):

        # iterate over the source triangles
        for jSrc in range(self.numTriangles):

            cellSrc = self.ptIdList[jSrc]
            paSrc = numpy.array(self.points.GetPoint(cellSrc[0]))
            pbSrc = numpy.array(self.points.GetPoint(cellSrc[1]))
            pcSrc = numpy.array(self.points.GetPoint(cellSrc[2]))

            pb2Src = pbSrc - paSrc
            pc2Src = pcSrc - paSrc
            areaSrcVec = numpy.cross(pb2Src, pc2Src)
            areaSrc = numpy.linalg.norm(areaSrcVec)
            normalSrc = areaSrcVec / areaSrc

            # iterate the observer triangles
            for iObs in range(self.numTriangles):

                cellObs = self.ptIdList[iObs]
                paObs = numpy.array(self.points.GetPoint(cellObs[0]))
                pbObs = numpy.array(self.points.GetPoint(cellObs[1]))
                pcObs = numpy.array(self.points.GetPoint(cellObs[2]))

                # observer is at mid point
                xObs = (paObs + pbObs + pcObs) / 3.0

                elev = (xObs - paSrc).dot(normalSrc)
                
                self.gMat[iObs, jSrc] = 0.0
                self.kMat[iObs, jSrc] = 0.0
               
                if iObs == jSrc:
            
                    # singular term
                    pot0ab = PotentialIntegrals(xObs, paSrc, pbSrc, self.order)
                    pot0bc = PotentialIntegrals(xObs, pbSrc, pcSrc, self.order)
                    pot0ca = PotentialIntegrals(xObs, pcSrc, paSrc, self.order)
        
                    self.gMat[iObs, jSrc] += pot0ab.getIntegralOneOverR(elev) + \
                                    pot0bc.getIntegralOneOverR(elev) + \
                                    pot0ca.getIntegralOneOverR(elev)
                    self.gMat[iObs, jSrc] /= FOUR_PI

                    # no contribution for kMat
        
                else:
        
                    # off diagonal term 
            
                    # Gauss wuadrature order estimate
                    normDistance = numpy.linalg.norm((paSrc + pbSrc + pcSrc)/3. - xObs) \
                        / numpy.sqrt(areaSrc)
                    offDiagonalOrder = min(8, max(1, int(8 * 2 / normDistance)))

                    # Gauss quadrature weights
                    gpws = gaussPtsAndWeights[offDiagonalOrder]

                    # number of Gauss points
                    npts = gpws.shape[1]

                    # triangle positions and weights
                    xsis, etas, weights = gpws[0, :], gpws[1, :], gpws[2, :]

                    for k in range(npts):
                        dr = xObs - paSrc - pb2Src*xsis[k] - pc2Src*etas[k]
                        drNorm = numpy.sqrt(dr.dot(dr))
                        self.gMat[iObs, jSrc] += weights[k] / drNorm
                        self.kMat[iObs, jSrc] += weights[k] * normalSrc.dot(dr) / drNorm**3

                    self.gMat[iObs, jSrc] *= 0.5 * areaSrc / FOUR_PI
                    self.kMat[iObs, jSrc] *= 0.5 * areaSrc / FOUR_PI


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

    def computeNeumannFromDirichlet(self, dirichletExpr):
        """
        Get the Neumann boundary values from the Dirichlet boundary conditions
        @param dirichletExpr expression for the potential values
        @return response
        """
        from math import pi, sin, cos, log, exp, sqrt
        
        # copy the matrix
        kMat = self.getNormalDerivativeGreenMatrix().copy()
        n = kMat.shape[0]
        
        v = numpy.zeros((n,), numpy.float64)
        
        # add residue. Note: Green function is 1/(4*pi*R), del^2 G = - delta
        for i in range(n):
            kMat[i, i] += 0.5

        pointArray = self.getPoints()
        cellArray = self.getCells()

        for i in range(n):
            ia, ib, ic = cellArray[i, :]
            x, y, z = (pointArray[ia, :] + pointArray[ib, :] + pointArray[ic, :])/3.            
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
        lslm = LaplaceMatrices(pdata,
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
        lslm = LaplaceMatrices(pdata,
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
        lslm = LaplaceMatrices(pdata,
                               max_edge_length=1000.,
                               order=order)
        print 'order = ', order
        print 'g matrix: ', lslm.getGreenMatrix()
        print 'k matrix: ', lslm.getNormalDerivativeGreenMatrix()

if __name__ == '__main__':
    #testSingleTriangle()
    #testTwoTrianglesCoplanar()
    testTwoTriangles()
