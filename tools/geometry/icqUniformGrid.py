#!/usr/bin/env python

import numpy

def uniformIndexGrid2D(loBound, hiBound):
    """
    Construct a linspace grid in 2D
    @param loBound low end of the box
    @param hiBound high end of the box
    @return ii, jj arrays
    """
    ns = numpy.array([hiBound[i] - loBound[i] for i in range(2)])
    ns1 = ns + 1
    xs = numpy.array([loBound[0] + i \
                      for i in range(ns1[0])], numpy.int)
    ys = numpy.array([loBound[1] + j \
                      for j in range(ns1[1])], numpy.int)
    
    xx = numpy.outer(xs, numpy.ones( (ns1[1],), numpy.int )).reshape(ns1)
    yy = numpy.outer(numpy.ones( (ns1[0],), numpy.int ), ys).reshape(ns1)
    
    return xx, yy


def uniformGrid2D(loBound, hiBound, ns):
    """
    Construct a linspace grid in 2D
    @param loBound low end of the box
    @param hiBound high end of the box
    @param ns number of cells in x, and y
    @return xx, yy arrays
    """
    ns1 = [n + 1 for n in ns]
    xs = numpy.linspace(loBound[0], hiBound[0], ns1[0])
    ys = numpy.linspace(loBound[1], hiBound[1], ns1[1])
    
    xx = numpy.outer(xs, numpy.ones( (ns1[1],) )).reshape(ns1)
    yy = numpy.outer(numpy.ones( (ns1[0],) ), ys).reshape(ns1)
                      
    return xx, yy

def uniformGrid3D(loBound, hiBound, ns):
  """
  Construct a linspace grid in 3D
  @param loBound low end of the box
  @param hiBound high end of the box
  @param ns number of cells in x, y, and z
  @return xxx, yyy, zzz arrays
  """
  ns1 = [n + 1 for n in ns]
  xs = numpy.linspace(loBound[0], hiBound[0], ns1[0])
  ys = numpy.linspace(loBound[1], hiBound[1], ns1[1])
  zs = numpy.linspace(loBound[2], hiBound[2], ns1[2])

  xxx = numpy.outer(xs, numpy.ones( (ns1[1], ns1[2]) )).reshape(ns1)
  yyy = numpy.outer( numpy.outer(numpy.ones((ns1[0],)), ys),
                     numpy.ones((ns1[2],)) ).reshape(ns1)
  zzz = numpy.outer(numpy.ones( (ns1[0], ns1[1]) ), zs).reshape(ns1)
  
  return xxx, yyy, zzz

def uniformGrid(loBound, hiBound, ns):

  """
  Construct a linspace grid
  @param loBound low end of the box
  @param hiBound high end of the box
  @param ns number of cells in x, y, and optionally in z
  @return xx, yy in 2D or xxx, yyy, and zzz arrays in 3D
  """
  ndims = len(loBound)
  if ndims == 3:
    return uniformGrid3D(loBound, hiBound, ns)
  elif ndims == 2:
    return uniformGrid2D(loBound, hiBound, ns)
  else:
    raise NotImplementedError

################################################################################
def test():
  xxx, yyy, zzz = uniformGrid(loBound=(0., 0., 0.), 
                              hiBound=(1., 2., 3.), ns=(1,2,3))
  print 'xxx = ', xxx.flat[:]
  print 'yyy = ', yyy.flat[:]
  print 'zzz = ', zzz.flat[:]

if __name__ == '__main__': test()




