#!/usr/bin/env python

"""
@brief A class for testing if a point is inside a shape
@author pletzer@psu.edu
"""

import numpy

class Inside:
    
    def __init__(self, points, polys, domainMins, domainMaxs):
        """
        Constructor
        @param points list of points
        @param polys connectivty
        """
        
        self.eps = 1.23456789e-14
        
        self.polys = polys
        self.points = points
        self.domainMins = domainMins
        self.domainMaxs = domainMaxs
        
        # must have at least one point
        self.ndims = len(points[0])
        
        # matrix and right and side vector
        self.mat = numpy.zeros( (self.ndims, self.ndims), numpy.float64 )
        self.b = numpy.zeros( (self.ndims,), numpy.float64 )
        
        # the optimal direction for shooting the ray
        self.direction = float('inf') * numpy.ones((self.ndims,), numpy.float64)
        
        # find the box corners
        self.xmins = +float('inf') * numpy.ones(self.ndims, numpy.float64)
        self.xmaxs = -float('inf') * numpy.ones(self.ndims, numpy.float64)
        for i in range(self.ndims):
            xmin = min([p[i] for p in points])
            xmax = max([p[i] for p in points])
            self.xmins[i] = min(self.xmins[i], xmin)
            self.xmaxs[i] = max(self.xmaxs[i], xmax)
    
    def isInside(self, point):
        """
        Determine if a point is inside the shape
        @param point point
        @return True if inside, False otherwise
        """
        
        # quick check, point must be inside box
        outsideBox = reduce(lambda x,y: x or y,
                      [(point[i] < self.xmins[i] - self.eps) or \
                       (point[i] > self.xmaxs[i] + self.eps)
                             for i in range(self.ndims)])
                             
        if outsideBox:
            return False
                             
        # any direction will do but things will run faster if the direction
        # points to the nearest domain box face (fewer intersections to
        #  compute)
        self.computeOptimalDirection(point)
                             
        # set the first column in our matrix (independent of the poly)
        self.mat[:, 0] = self.direction
                             
        # compute the number of intersections between the ray and the polys
        numberOfIntersections = 0
        for poly in self.polys:
                                 
            # compute the overlap betwen the ray box and the face box
            if not self.areBoxesOverlapping(point, poly):
                # intersection is impossible, don't bother...
                continue
                                 
            lmbda, xis = self.computeIntersection(point, poly)
            if lmbda > 0.0 + self.eps:
                # the direction is ok
                sums = 0.0
                rayIntersects = True
                for i in range(len(xis)):
                    rayIntersects &= (xis[i] > self.eps) and \
                                (xis[i] < 1 - sums - self.eps)
                    sums += xis[i]
                    if rayIntersects:
                        numberOfIntersections += 1
                             
        # even number is outside (False), odd number means inside (True)
        return (numberOfIntersections % 2)
    
    def areBoxesOverlapping(self, point, poly):
        
        for i in range(self.ndims):
            
            xminFace = min( [self.points[j][i] for j in poly] )
            xmaxFace = max( [self.points[j][i] for j in poly] )
            
            p, d = point[i], self.direction[i]
            xminRay = min(p, p + d)
            xmaxRay = max(p, p + d)
            
            if (xmaxRay < xminFace + self.eps) or \
                (xmaxFace < xminRay + self.eps):
                    # no overlap
                    return False
        return True
    
    def computeOptimalDirection(self, point):
        """
        Compute the direction of the ray to to the nearest
        domain box face. This will update self.direction
        @param point starting point of the ray
        """
        # iterate over the faces of the box then find the minimum distance between
        # the point and the face.
        minDistance = float('inf')
        
        for pm in (-1, 1):
            pls = (1 + pm)/2.
            mns = (1 - pm)/2.
            for axis in range(self.ndims):
                # the normal vector contains very small values in place of
                # zero elements
                
                #  in order to avoid issues with ray hitting exactly a node
                normal = 100 * self.eps * numpy.random.rand( self.ndims )
                normal[axis] = pm
                distance = pls*(self.domainMaxs[axis] - point[axis]) + mns*(point[axis] - self.domainMins[axis])
                if distance < minDistance:
                    # expand a little beyond the domain (1.1)
                    self.direction = normal * (1.1* distance)
                    minDistance = distance
    
    def computeIntersection(self, point, poly):
        self.b = self.points[poly[0]] - point
        for i in range(self.ndims - 1):
            self.mat[:, 1 + i] = self.points[poly[0]] - self.points[poly[1 + i]]
        solution = numpy.linalg.solve(self.mat, self.b)
        return solution[0], solution[1:]


#################################################################################

def test2d():
    points = [numpy.array([0., 0.]),
              numpy.array([1., 0.]),
              numpy.array([0., 1.])]
    polys = [(0, 1), (1, 2), (2, 0)]
    inside = Inside(points, polys,
                    domainMins=numpy.array([0., 0.]),
                    domainMaxs=numpy.array([1., 1.]))

    # really outside
    assert(not inside.isInside(numpy.array([-0.1, -0.2])))

    # inside
    assert(inside.isInside(numpy.array([+0.1, +0.2])))

    # outside
    assert(not inside.isInside(numpy.array([-0.1, +0.2])))
    assert(not inside.isInside(numpy.array([+0.1, -0.2])))

    # on the face
    assert(not inside.isInside(numpy.array([-0.000001, +0.2])))

if __name__ == '__main__':
    test2d()
