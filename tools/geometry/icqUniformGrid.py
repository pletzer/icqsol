#!/usr/bin/env python

import numpy

def uniformGrid(loBound, hiBound, ns):
  """
  Construct a linspace grid
  @param loBound low end of the box
  @param hiBound high end of the box
  @param ns number of cells in x, y, and z
  """
  xs = numpy.linspace(loBound[0],hiBound[0], ns[0])
  ys = numpy.linspace(loBound[1],hiBound[1], ns[1])
  zs = numpy.linspace(loBound[2],hiBound[2], ns[2])

  xxx = numpy.outer(numpy.ones( (ns[2], ns[1]) ), xs)
  yyy = numpy.outer( numpy.outer(numpy.ones((ns[2],)), ys), numpy.ones((ns[0],)) )
  zzz = numpy.outer(zs, numpy.ones( (ns[1], ns[0]) ))
  
  return xxx, yyy, zzz
 



