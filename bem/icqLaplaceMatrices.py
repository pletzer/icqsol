#!/usr/bin/env python

import vtk
import numpy
from icqsol.shapes.icqRefineSurface import RefineSurface
from icqsol.bem.icqQuadrature import triangleQuadrature
from scipy.special import sph_harm


class SrcFunc:

    """
    Source functor (singular)
    """

    def __init__(self, m, n, center):
        self.m = m
        self.n = n
        self.center = center

    def __call__(self, xp):
        x = xp - self.center
        r = numpy.sqrt(x.dot(x))
        rho = numpy.sqrt(x[0]**2 + x[1]**2)
        theta = numpy.arcsin(rho/r)
        phi = numpy.arccos(x[0]/rho)
        # assumes r is further away than observer from center
        return sph_harm(self.m, self.n, phi, theta).conjugate() / \
            r**(self.n + 1)


class SrcNormalDerivFun:

    """
    Source functor for normal derivative of Green function
    """
    def __init__(self, normal, xobs):
        self.normal = normal
        self.xobs = xobs
    
    def __call__(self, x):
        d = self.xobs - x
        r = numpy.sqrt(d.dot(d)) + 1.e-10
        return self.normal.dot(d)/(4.*numpy.pi * r**1.5)


class ObsFunc:
    """
    Observer functor
    """

    def __init__(self, m, n, center):
        self.m = m
        self.n = n
        self.center = center

    def __call__(self, xp):
        x = xp - self.center
        r = numpy.sqrt(x.dot(x))
        rho = numpy.sqrt(x[0]**2 + x[1]**2)
        theta = numpy.arcsin(rho/r)
        phi = numpy.arccos(x[0]/rho)
        # assumes r is closer than any source location from center
        return sph_harm(self.m, self.n, phi, theta) * \
            r**(self.n)/(2*self.n + 1)


class LaplaceMatrices:

    def __init__(self, pdata, max_edge_length, order, maxN):
        """
        Constructor
        @param pdata instance of vtkPolyData
        @param max_edge_length maximum edge length, used to turn
                               polygons into triangles
        @param order order of the expansion (in the range 1 to 8)
        @param maxN maximum n in spherical harmonic expansion
                            (> 0, ~ 10 is a good number)
        """

        # triangulate
        rs = RefineSurface(pdata)
        rs.refine(max_edge_length=max_edge_length)
        self.pdata = rs.getVtkPolyData()

        polys = self.pdata.GetPolys()
        numTriangles = polys.GetNumberOfCells()

        shp = (numTriangles, numTriangles)
        self.singleLayerMatrix = numpy.zeros(shp, numpy.complex)
        self.doubleLayerMatrix = numpy.zeros(shp, numpy.complex)
        self.__computeMatricesByExpansion(order, maxN)

    def __computeMatricesByExpansion(self, order, maxN):

        polys = self.pdata.GetPolys()
        points = self.pdata.GetPoints()
        numTriangles = polys.GetNumberOfCells()

        ptIdsObs = vtk.vtkIdList()
        ptIdsSrc = vtk.vtkIdList()

        # iterate over observer points
        polys.InitTraversal()
        for iTriangleObs in range(numTriangles):

            polys.GetNextCell(ptIdsObs)

            ia = ptIdsObs.GetId(0)
            ib = ptIdsObs.GetId(1)
            ic = ptIdsObs.GetId(2)
            paObs = numpy.array(points.GetPoint(ia))
            pbObs = numpy.array(points.GetPoint(ib))
            pcObs = numpy.array(points.GetPoint(ic))

            # triangle mid point is the observer position
            pObs = (paObs + pbObs + pcObs)/3.

            # edges
            dp1Obs = pbObs - paObs
            dp2Obs = pcObs - paObs
            dp3Obs = pcObs - pbObs

            # normal to the triangle
            normal = numpy.cross(dp1Obs, dp2Obs)
            normal /= numpy.sqrt(normal.dot(normal))

            # max edge length
            edgeLength = numpy.sqrt(max(dp1Obs.dot(dp1Obs),
                                        dp2Obs.dot(dp2Obs),
                                        dp3Obs.dot(dp3Obs)))

            # set the coordinate reference position distance >~ cell size.
            # Make it too small and the quadrature will not converge.
            # Make it too large and many more expansion terms are needed

            distance = 1.5
            center = pObs + distance*edgeLength*normal

            # iterate over source triangles
            for iTriangleSrc in range(numTriangles):
            
                polys.GetCell(iTriangleSrc, ptIdsSrc)
            
                ia = ptIdsSrc.GetId(0)
                ib = ptIdsSrc.GetId(1)
                ic = ptIdsSrc.GetId(2)
                paSrc = numpy.array(points.GetPoint(ia))
                pbSrc = numpy.array(points.GetPoint(ib))
                pcSrc = numpy.array(points.GetPoint(ic))
                #print '*** paSrc, pbSrc, pcSrc = ', paSrc, pbSrc, pcSrc
            
                normalSrc = numpy.cross(pbSrc - paSrc, pcSrc - paSrc)
                normalSrc /= numpy.sqrt(normalSrc.dot(normalSrc))
            
                # normal derivative of Green function
                kSrc = SrcNormalDerivFun(normalSrc, pObs)
                
                beta = triangleQuadrature(order,
                                          paSrc, pbSrc, pcSrc,
                                          kSrc)
                self.doubleLayerMatrix[iTriangleObs, iTriangleSrc] = beta

                # expand in spherical harmonics
                for n in range(maxN):
                    for m in range(-n, n+1):

                        fObs = ObsFunc(m, n, center)
                        fObsVal = fObs(pObs)

                        # Green functor
                        gSrc = SrcFunc(m, n, center)

                        # evaluate the integrals
                        alpha = triangleQuadrature(order,
                                                   paSrc, pbSrc, pcSrc,
                                                   gSrc)
                        
                        self.singleLayerMatrix[iTriangleObs, iTriangleSrc] += \
                            fObsVal*alpha
    

    def getSingleLayerMatrix(self):
        """
        Return the single layer coupling matrix
        @return matrix
        """
        return self.singleLayerMatrix

    def getDoubleLayerMatrix(self):
        """
        Return the double layer coupling matrix
        @return matrix
        """
        return self.doubleLayerMatrix

    def computeNeumannFromDirichlet(self, dirichletExpr):
        """
        Get the Neumann boundary values from the Dirichlet boundary conditions
        @param dirichletExpr expression for the potential values
        @return response
        """
        from math import pi, sin, cos, log, exp, sqrt
        
        kMat = self.getDoubleLayerMatrix()
        n = kMat.shape[0]
        
        v = numpy.zeros((n,), numpy.complex)
        # add residue
        for i in range(n):
            kMat[i, i] += 0.5
        
        # evaluate potential on each triangle center
        ptIds = vtk.vtkIdList()
        polys = self.pdata.GetPolys()
        points = self.pdata.GetPoints()
        for i in range(n):
            
            polys.GetCell(i, ptIds)
            ia, ib, ic = ptIds.GetId(0), ptIds.GetId(1), ptIds.GetId(2)
            pa = numpy.array(points.GetPoint(ia))
            pb = numpy.array(points.GetPoint(ib))
            pc = numpy.array(points.GetPoint(ic))
            pMid = (pa + pb + pc)/3.
            x, y, z = pMid
            
            v[i] = eval(dirichletExpr)
        
        gMat = self.getSingleLayerMatrix()
        print 'g = ', gMat
        print 'k = ', kMat
        
        # solve
        gM1 = numpy.linalg.inv(gMat)
        return gM1.dot(kMat).dot(v)


###############################################################################


def test():

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

    for order in range(1, 9):
        lslm = LaplaceMatrices(pdata,
                               max_edge_length=1000.,
                               order=order, maxN=20)
        print 'order = ', order
        print 'single layer matrix: ', lslm.getSingleLayerMatrix()
        print 'double layer matrix: ', lslm.getDoubleLayerMatrix()

if __name__ == '__main__':
    test()
