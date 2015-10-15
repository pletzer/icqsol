#!/usr/bin/env python

import vtk
import numpy
from icqsol.shapes.icqRefineSurface import RefineSurface
from icqsol.bem.icqQuadrature import triangleQuadrature
from scipy.special import sph_harm

class SrcFunc:
    
    def __init__(self, m, n):
        self.m = m
        self.n = n
    
    def __call__(self, x):
        r = numpy.sqrt(x.dot(x))
        rho = numpy.sqrt(x[0]**2 + x[1]**2)
        theta = numpy.arcsin(rho/r)
        phi = numpy.arccos(x[0]/rho)
        return sph_harm(self.m, self.n, phi, theta).conjugate()/r**(self.n + 1)

class ObsFunc:
    
    def __init__(self, m, n):
        self.m = m
        self.n = n
    
    def __call__(self, x):
        r = numpy.sqrt(x.dot(x))
        rho = numpy.sqrt(x[0]**2 + x[1]**2)
        theta = numpy.arcsin(rho/r)
        phi = numpy.arccos(x[0]/rho)
        return sph_harm(self.m, self.n, phi, theta) * r**(self.n)/(2*self.n + 1)


class LaplaceSingleLayerMatrix:

    def __init__(self, pdata, max_edge_length, order, maxN):
        """
        Constructor
        @param pdata instance of vtkPolyData
        """
        # triangulate
        rs = RefineSurface(pdata)
        rs.refine(max_edge_length=max_edge_length)
        self.pdata = rs.getVtkPolyData()
        
        numPoints = self.pdata.GetPoints().GetNumberOfPoints()
        polys = self.pdata.GetPolys()
        numTriangles = polys.GetNumberOfCells()
        
        self.matrix = numpy.zeros((numTriangles, numTriangles), numpy.complex)

        ptIds = vtk.vtkIdList()
        ptIds2 = vtk.vtkIdList()
        points = self.pdata.GetPoints()
        polys.InitTraversal()
        # iterate over the source triangles
        for iTriangle in range(numTriangles):
            polys.GetNextCell(ptIds)
            ia, ib, ic = ptIds.GetId(0), ptIds.GetId(1), ptIds.GetId(2)
            pa = numpy.array(points.GetPoint(ia))
            pb = numpy.array(points.GetPoint(ib))
            pc = numpy.array(points.GetPoint(ic))
            # expand in spherical harmonics
            for n in range(maxN):
                for m in range(-n, n+1):
                    f = SrcFunc(m, n)
                    f2 = ObsFunc(m, n)
                    alphaMN = triangleQuadrature(order, pa, pb, pc, f)
                    # iterate over the observer positions
                    for iTriangle2 in range(numTriangles):
                        polys.GetCell(iTriangle2, ptIds2)
                        ia2, ib2, ic2 = ptIds2.GetId(0), ptIds2.GetId(1), ptIds2.GetId(2)
                        pa2 = numpy.array(points.GetPoint(ia2))
                        pb2 = numpy.array(points.GetPoint(ib2))
                        pc2 = numpy.array(points.GetPoint(ic2))
                        
                        normal = numpy.cross(pb2 - pa2, pc2 - pa2)
                        normal /= numpy.sqrt(normal.dot(normal))
                        # observer is at mid point
                        edgeLength = numpy.sqrt(max((pb2-pa2).dot(pb2-pa2), 
                                                    (pc2-pb2).dot(pc2-pb2),
                                                    (pa2-pc2).dot(pa2-pc2)))
                        pMid = (pa2 + pb2 + pc2)/3. - edgeLength*0.0*normal

                        self.matrix[iTriangle2, iTriangle] += f2(pMid)*alphaMN

    def getMatrix(self):
        return self.matrix

#############################################################################################

def test():
    
    # create set of points
    points = vtk.vtkPoints()
    points.SetNumberOfPoints(3)
    points.SetPoint(0, [1., -1./3., -1./3.])
    points.SetPoint(1, [1., 2./3., -1./3.])
    points.SetPoint(2, [1., -1./3., 2./3.])

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

    for order in range(2, 9):
        lslm = LaplaceSingleLayerMatrix(pdata, max_edge_length=1000., order=order, maxN=10)
        print 'order = ', order, ' matrix = ', lslm.getMatrix()

if __name__ == '__main__':
    test()