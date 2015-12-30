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
        
        self.normalDerivativeJumpName = 'normal_derivative_jump'
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
        
        self.order = order
        self.__computeMatrices()

    def setPotentialName(self, name):
        self.potentialName = name

    def setNormalDerivativeJumpName(self, name):
        self.normalDerivativeJumpName             
            
    def getVtkPolyData(self):
        """
        Get the (modified) vtkPolyData object
        @return object
        """
        return self.pdata

    def __computeMatrices(self):

        # iterate over the source triangles
        for jSrc in range(self.numTriangles):

            ia, ib, ic = self.ptIdList[jSrc]

            # The triangle vertex positions
            paSrc = numpy.array(self.points.GetPoint(ia))
            pbSrc = numpy.array(self.points.GetPoint(ib))
            pcSrc = numpy.array(self.points.GetPoint(ic))

            # The triangle's normal vector and area at the center of the triangle
            pb2Src = pbSrc - paSrc
            pc2Src = pcSrc - paSrc
            areaSrcVec = numpy.cross(pb2Src, pc2Src)
            areaSrc = numpy.linalg.norm(areaSrcVec)
            normalSrc = areaSrcVec / areaSrc

            # iterate over the observer triangles
            for iObs in range(self.numTriangles):

                cellObs = self.ptIdList[iObs]
                paObs = numpy.array(self.points.GetPoint(cellObs[0]))
                pbObs = numpy.array(self.points.GetPoint(cellObs[1]))
                pcObs = numpy.array(self.points.GetPoint(cellObs[2]))

                # Observer is at mid point
                xObs = (paObs + pbObs + pcObs) / 3.0
                
                if iObs == jSrc:
            
                    # Singular term
                    pot0ab = PotentialIntegrals(xObs, paSrc, pbSrc, self.order)
                    pot0bc = PotentialIntegrals(xObs, pbSrc, pcSrc, self.order)
                    pot0ca = PotentialIntegrals(xObs, pcSrc, paSrc, self.order)
        
                    self.gMat[iObs, jSrc] = pot0ab.getIntegralOneOverR() + \
                                            pot0bc.getIntegralOneOverR() + \
                                            pot0ca.getIntegralOneOverR()
                    self.gMat[iObs, jSrc] /= (-FOUR_PI)
        
                else:
        
                    #
                    # Off diagonal term
                    # 
            
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
                    self.gMat[iObs, jSrc] *= 0.5 * areaSrc / (-FOUR_PI)

    def getGreenMatrix(self):
        """
        Return the Green function matrix
        @return matrix
        """
        return self.gMat

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

    def computeNeumannJumpFromDirichlet(self, dirichletExpr):
        """
        Get the jump of the normal potential derivative from the Dirichlet boundary conditions
        @param dirichletExpr expression for the potential values
        @return response
        """
        from math import pi, sin, cos, log, exp, sqrt
                
        
        pointArray = self.getPoints()
        cellArray = self.getCells()
        n = cellArray.shape[0]

        # Set the potential.
        v = numpy.zeros((n,), numpy.float64)
        for i in range(n):
            ia, ib, ic = cellArray[i, :]
            x, y, z = (pointArray[ia, :] + pointArray[ib, :] + pointArray[ic, :])/3.            
            # set the value
            v[i] = eval(dirichletExpr)
        
        gMat = self.getGreenMatrix()
        
        normalDerivativeJump = numpy.linalg.inv(gMat).dot(v)

        # add field
        cellData = self.pdata.GetCellData()
        
        potentialData = vtk.vtkDoubleArray()
        potentialData.SetNumberOfComponents(1)
        potentialData.SetNumberOfTuples(n)
        potentialData.SetName(self.potentialName)

        normalDerivData = vtk.vtkDoubleArray()
        normalDerivData.SetNumberOfComponents(1)
        normalDerivData.SetNumberOfTuples(n)
        normalDerivData.SetName(self.normalDerivativeJumpName)

        for i in range(n):
            potentialData.SetTuple(i, [v[i],])
            normalDerivData.SetTuple(i, [normalDerivativeJump[i],])

        cellData.AddArray(potentialData)
        cellData.AddArray(normalDerivData)
        
        return normalDerivativeJump




###############################################################################


def testSingleTriangle():

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

if __name__ == '__main__':
    testSingleTriangle()
    testTwoTrianglesCoplanar()
    testTwoTriangles()
