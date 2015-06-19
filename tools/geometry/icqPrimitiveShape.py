# standard python modules
import re

# extensions
import vtk
import numpy

from icqBaseShape import BaseShape

class PrimitiveShape(BaseShape):

  def __init__(self,):
    """
    Constructor
    """
    self.surfaceFuncs = []
    self.evalFunc = None

  def setSurfaceFunctions(self, funcs):
    """
    Set the surface functions
    @param funcs a list of 3-tuples, [(x(u,v), y(u,v), z(u,v)), ...] 
           with 0 <= u, v <= 1
    """
    self.surfaceFuncs = funcs

  def setEvaluateFunction(self, func):
    """
    Set the function of (x, y, z) that returns > 0 inside and 0 outside
    """
    self.evalFunc = func

  def evaluate(self, pts):
    """
    Evaluate the inside/outside function 
    @param pts array of points
    """
    return self.evalFunc(pts)

################################################################################
def test():

  # create a cylinder terminated by two half-spheres
  radius = 1.0
  radius2 = radius**2
  length = 0.5

  #
  # define the surface of the cylinder
  #
  def xCylSurf(u, v):
    return radius*numpy.cos(2*numpy.pi*u)

  def yCylSurf(u, v):
    return length*v

  def zCylSurf(u, v):
    return radius*numpy.sin(2*numpy.pi*u)

  #
  # define the low end surface
  #
  def xLoSurf(u, v):
    return -radius*numpy.sin(numpy.pi*u)*numpy.cos(numpy.pi*v)

  def yLoSurf(u, v):
    return -radius*numpy.sin(numpy.pi*u)*numpy.sin(numpy.pi*v)

  def zLoSurf(u, v):
    return radius*numpy.cos(numpy.pi*u)

  #
  # define the high end surface
  #
  def xHiSurf(u, v):
    return -radius*numpy.sin(numpy.pi*u)*numpy.cos(numpy.pi*v)

  def yHiSurf(u, v):
    return length + radius*numpy.sin(numpy.pi*u)*numpy.sin(numpy.pi*v)

  def zHiSurf(u, v):
    return radius*numpy.cos(numpy.pi*u)

  def objectVolume(x, y, z):
    rho2 = x**2 + z**2
    r12 = rho2 + y**2
    r22 = rho2 + (y-length)**2
    res = y < length
    res &= y > 0
    res |= (r12 < radius2)
    res |= (r22 < radius2)
    return res

  s = PrimitiveShape()
  s.setSurfaceFunctions([(xCylSurf, yCylSurf, zCylSurf), 
                         (xLoSurf, yLoSurf, zLoSurf),
                         (xHiSurf, yHiSurf, zHiSurf)])
  s.setEvaluateFunction(objectVolume)

if __name__ == '__main__': test()
