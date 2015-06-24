#!/usr/bin/env python

import numpy

def uniformGrid(loBound, hiBound, ns):
  """
  Construct a linspace grid
  @param loBound low end of the box
  @param hiBound high end of the box
  @param ns number of cells in x, y, and z
  """
  ns1 = [n + 1 for n in ns]
  xs = numpy.linspace(loBound[0], hiBound[0], ns1[0])
  ys = numpy.linspace(loBound[1], hiBound[1], ns1[1])
  zs = numpy.linspace(loBound[2], hiBound[2], ns1[2])

  xxx = numpy.outer(numpy.ones( (ns1[2], ns1[1]) ), xs).reshape(ns1)
  yyy = numpy.outer( numpy.outer(numpy.ones((ns1[2],)), ys), 
                     numpy.ones((ns1[0],)) ).reshape(ns1)
  zzz = numpy.outer(zs, numpy.ones( (ns1[1], ns1[0]) )).reshape(ns1)
  
  return xxx, yyy, zzz

################################################################################
def test():
  xxx, yyy, zzz = uniformGrid(loBound=(0., 0., 0.), 
                              hiBound=(1., 2., 3.), ns=(1,2,3))
  print 'xxx = ', xxx.flat[:]
  print 'yyy = ', yyy.flat[:]
  print 'zzz = ', zzz.flat[:]

if __name__ == '__main__': test()




