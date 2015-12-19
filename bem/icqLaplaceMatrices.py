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
                            (>= 0, ~10 is a good number)
        """

        assert(order > 0 and order <= 8)
        assert(maxN >= 0)

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
        self.maxN = maxN
        self.__computeMatrices()
    
    def __computeNormalDerivativeGreenCoupling(self, iObs, jSrc):
        
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
        
        normalSrc = numpy.cross(pbSrc - paSrc, pcSrc - paSrc)
        lengthSqr = normalSrc.dot(normalSrc)
        assert(lengthSqr > 0.)
        normalSrc /= numpy.sqrt(lengthSqr)

        if iObs != jSrc:
            
            # use standard Gauss quadrature
            
            def kreen(x):
                r = xObs - x
                return normalSrc.dot(r)/(4. * numpy.pi * numpy.sqrt(r.dot(r))**3)
            
            self.kMat[iObs, jSrc] = triangleQuadrature(self.order,
                                                       paSrc, pbSrc, pcSrc,
                                                                    kreen)
        
        else:

            # zero contribution [normal and (xObs - xSrc) are orthogonal]
            self.kMat[iObs, jSrc] = 0.0
            

    def __computeGreenCoupling(self, iObs, jSrc):
                
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

        if iObs != jSrc:
        
            # use standard Gauss quadrature
            
            def green(x):
                r = xObs - x
                return 1.0/(4. * numpy.pi * numpy.sqrt(r.dot(r)))
        
            alpha = triangleQuadrature(self.order, paSrc, pbSrc, pcSrc, green)
            self.gMat[iObs, jSrc] = alpha
        
        else:
            
            # quadrature by expansion in spherical harmonics
            
            # normal to the triangle
            dp1Obs = pbObs - paObs
            dp2Obs = pcObs - paObs
            dp3Obs = pcObs - pbObs
            normal = numpy.cross(dp1Obs, dp2Obs)
            lengthSqr = normal.dot(normal)
            assert(lengthSqr > 0.)
            normal /= numpy.sqrt(lengthSqr)
            
            # max edge length of the source triangle
            edgeLength = numpy.sqrt(max(dp1Obs.dot(dp1Obs),
                                        dp2Obs.dot(dp2Obs),
                                        dp3Obs.dot(dp3Obs)))
                                        
            # set the coordinate reference position distance >~ cell size.
            # Make it too small and the quadrature will not converge.
            # Make it too large and many more expansion terms are needed

            distance = 1.5
            center = xObs + distance*edgeLength*normal

            alpha = 0.0
            for n in range(self.maxN):
                for m in range(-n, n+1):
                    
                    fObs = ObsFunc(m, n, center)
                    fObsVal = fObs(xObs)
                    
                    # Green functor
                    gSrc = SrcFunc(m, n, center)
                    
                    # evaluate the integrals
                    alpha += fObsVal * \
                                triangleQuadrature(self.order,
                                                   paSrc, pbSrc, pcSrc,
                                                   gSrc)

                    # add contribution
            self.gMat[iObs, jSrc] = alpha.real

    def __computeMatrices(self):

        # iterate over the observer triangles
        for iObs in range(self.numTriangles):

            # iterate over source triangles
            for jSrc in range(self.numTriangles):
                
                self.__computeGreenCoupling(iObs, jSrc)
                self.__computeNormalDerivativeGreenCoupling(iObs, jSrc)

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
        # add residue
        for i in range(n):
            kMat[i, i] -= 0.5
        
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
        
        gMat = self.getGreenMatrix()
        print 'g = ', gMat
        print 'k = ', kMat
        
        # solve
        gM1 = numpy.linalg.inv(gMat)
        return gM1.dot(kMat).dot(v)


###############################################################################


def test():

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

    for order in range(1, 9):
        lslm = LaplaceMatrices(pdata,
                               max_edge_length=1000.,
                               order=order, maxN=20)
        print 'order = ', order
        print 'g matrix: ', lslm.getGreenMatrix()
        print 'k matrix: ', lslm.getNormalDerivativeGreenMatrix()


def test2():

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

    for order in range(1, 9):
        lslm = LaplaceMatrices(pdata,
                               max_edge_length=1000.,
                               order=order, maxN=20)
        print 'order = ', order
        print 'g matrix: ', lslm.getGreenMatrix()
        print 'k matrix: ', lslm.getNormalDerivativeGreenMatrix()

if __name__ == '__main__':
    #test()
    test2()
