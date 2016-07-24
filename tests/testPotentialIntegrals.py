#!/usr/bin/env python

from icqsol.bem.icqPotentialIntegrals import PotentialIntegrals
from icqsol.bem.icqLaplaceMatrices import getFullyQualifiedSharedLibraryName
import numpy

def testObserverOnA(order):
    paSrc = numpy.array([0., 0., 0.])
    pbSrc = numpy.array([1., 0., 0.])
    pcSrc = numpy.array([0., 1., 0.])
    xObs = paSrc
    integral = PotentialIntegrals(xObs, pbSrc, pcSrc, order).getIntegralOneOverR()
    exact = numpy.sqrt(2.) * numpy.arcsinh(1.)
    print('testObserverOnA: order = {0} integral = {1} exact = {2} error = {3}'.format(\
        order, integral, exact, integral - exact))

def testObserverOnB(order):
    paSrc = numpy.array([0., 0., 0.])
    pbSrc = numpy.array([1., 0., 0.])
    pcSrc = numpy.array([0., 1., 0.])
    xObs = pbSrc
    integral = PotentialIntegrals(xObs, pcSrc, paSrc, order).getIntegralOneOverR()
    exact = numpy.arcsinh(1.)
    print('testObserverOnB: order = {0} integral = {1} exact = {2} error = {3}'.format(\
        order, integral, exact, integral - exact))

def testObserverOnC(order):
    paSrc = numpy.array([0., 0., 0.])
    pbSrc = numpy.array([1., 0., 0.])
    pcSrc = numpy.array([0., 1., 0.])
    xObs = pcSrc
    integral = PotentialIntegrals(xObs, paSrc, pbSrc, order).getIntegralOneOverR()
    exact = numpy.arcsinh(1.)
    print('testObserverOnC: order = {0} integral = {1} exact = {2} error = {3}'.format(\
        order, integral, exact, integral - exact))

def testOffDiagonal2Triangles():
    import vtk
    import sys
    import pkg_resources
    PY_MAJOR_VERSION = sys.version_info[0]
    from ctypes import cdll, c_long, POINTER, c_double

    pdata = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    points.SetNumberOfPoints(4)
    points.SetPoint(0, (0., 0., 0.))
    points.SetPoint(1, (1., 0., 0.))
    points.SetPoint(2, (1., 1., 0.))
    points.SetPoint(3, (0., 1., 0.))
    pdata.SetPoints(points)
    pdata.Allocate(2, 1)
    ptIds = vtk.vtkIdList()
    ptIds.SetNumberOfIds(3)
    ptIds.SetId(0, 0); ptIds.SetId(1, 1); ptIds.SetId(2, 3)
    pdata.InsertNextCell(vtk.VTK_POLYGON, ptIds)
    ptIds.SetId(0, 1); ptIds.SetId(1, 2); ptIds.SetId(2, 3)
    pdata.InsertNextCell(vtk.VTK_POLYGON, ptIds)
    addr = int(pdata.GetAddressAsString('vtkPolyData')[5:], 0)
    gMat = numpy.zeros((2,2), numpy.float64)
    if PY_MAJOR_VERSION < 3:
        fullyQualifiedLibName = pkg_resources.resource_filename('icqsol', 'icqLaplaceMatricesCpp.so')
    else:
        libName = pkg_resources.resource_filename('icqsol', 'icqLaplaceMatricesCpp')
        fullyQualifiedLibName = getFullyQualifiedSharedLibraryName(libName)
    lib = cdll.LoadLibrary(fullyQualifiedLibName)
    lib.computeOffDiagonalTerms(c_long(addr),
                                gMat.ctypes.data_as(POINTER(c_double)))
    exact = numpy.array([[0, -0.07635909342383773],[-0.07635909342383773, 0]])
    print(gMat)
    print(exact)
    print('error: {0}'.format(gMat - exact))

if __name__ == '__main__':
    for order in range(1, 6):
        testObserverOnA(order)
        testObserverOnB(order)
        testObserverOnC(order)
        print('-'*80)
    testOffDiagonal2Triangles()

