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
        return sph_harm(self.m, self.n, phi, theta).conjugate()/r**(self.n + 1)

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
        return sph_harm(self.m, self.n, phi, theta) * r**(self.n)/(2*self.n + 1)


class LaplaceSingleLayerMatrix:

    def __init__(self, pdata, max_edge_length, order, maxN):
        """
        Constructor
        @param pdata instance of vtkPolyData
        @param max_edge_length maximum edge length, used to turn polygons into triangles
        @param order order of the expansion (in the range 1 to 8)
        @param maxN maximum n in spherical harmonic expansion (> 0, ~ 10 is a good number)
        """
        
        # triangulate
        rs = RefineSurface(pdata)
        rs.refine(max_edge_length=max_edge_length)
        self.pdata = rs.getVtkPolyData()
        
        numPoints = self.pdata.GetPoints().GetNumberOfPoints()
        polys = self.pdata.GetPolys()
        numTriangles = polys.GetNumberOfCells()
        
        self.matrix = numpy.zeros((numTriangles, numTriangles), numpy.complex)
        self.__computeMatrixByExpansion(order, maxN)
    
   
        
    def __computeMatrixByExpansion(self, order, maxN):
        
        polys = self.pdata.GetPolys()
        points = self.pdata.GetPoints()
        numTriangles = polys.GetNumberOfCells()
        
        ptIdsObs = vtk.vtkIdList()
        ptIdsSrc = vtk.vtkIdList()

        # iterate over observer points
        polys.InitTraversal()
        for iTriangleObs in range(numTriangles):
            
            polys.GetNextCell(ptIdsObs)
            
            ia, ib, ic = ptIdsObs.GetId(0), ptIdsObs.GetId(1), ptIdsObs.GetId(2)
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
                                  
            # set the coordinate reference position
            # distance >~ cell size, make it too small and the quadrature will not
            # converge. Make it too large and many more expansion terms are needed
            
            distance = 1.5
            center = pObs + distance*edgeLength*normal
        
            # expand in spherical harmonics
            for n in range(maxN):
                for m in range(-n, n+1):
    
                    fObs = ObsFunc(m, n, center)
                    fObsVal = fObs(pObs)
                    
                    fSrc = SrcFunc(m, n, center)
    
                    # iterate over source triangles
                    for iTriangleSrc in range(numTriangles):
                        
                        polys.GetCell(iTriangleSrc, ptIdsSrc)
                        
                        ia, ib, ic = ptIdsSrc.GetId(0), ptIdsSrc.GetId(1), ptIdsSrc.GetId(2)
                        paSrc = numpy.array(points.GetPoint(ia))
                        pbSrc = numpy.array(points.GetPoint(ib))
                        pcSrc = numpy.array(points.GetPoint(ic))
    
                        # evaluate the integral
                        alphaMN = triangleQuadrature(order, paSrc, pbSrc, pcSrc, fSrc)
                        self.matrix[iTriangleObs, iTriangleSrc] += fObsVal * alphaMN
        


    def getMatrix(self):
        """
        Return the coupling matrix
        @return matrix
        """
        return self.matrix

#############################################################################################

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
        lslm = LaplaceSingleLayerMatrix(pdata, max_edge_length=1000., order=order, maxN=20)
        print 'order = ', order, ' matrix = ', lslm.getMatrix()

if __name__ == '__main__':
    test()