#!/usr/bin/env python

"""
@brief A class for testing if a point is inside a shape
@author alexander@gokliya.net
"""

from __future__ import print_function
import numpy
import math
from functools import reduce


class Inside:

    def __init__(self, shape):
        """
        Constructor
        @param shape instance of Shape
        @param polys connectivty
        """
        points, polys, count = shape.toVerticesAndPolygons()
        # must have at least one point
        self.ndims = len(points[0])

        # find the box corners
        self.xmins = numpy.array([min([p[i] for p in points]) for
                                 i in range(self.ndims)])
        self.xmaxs = numpy.array([max([p[i] for p in points]) for
                                 i in range(self.ndims)])

        self.eps = 1.23456789e-14

        self.polys = polys
        self.points = [numpy.array(p) for p in points]

        # matrix and right and side vector
        self.mat = numpy.zeros((self.ndims, self.ndims), numpy.float64)
        self.b = numpy.zeros((self.ndims,), numpy.float64)

        # the optimal direction for shooting the ray
        self.direction = float('inf') * \
            numpy.ones((self.ndims,), numpy.float64)

        # min/max box corners of the ray
        self.xminRay = float('inf') * numpy.ones((self.ndims,), numpy.float64)
        self.xmaxRay = float('inf') * numpy.ones((self.ndims,), numpy.float64)

    def isInside(self, point, minDistance):
        """
        Determine if a point is inside the shape
        @param point point
        @param minDistance a point is declared inside if is at least minDistance
               away from the face
        @return +1 if inside, -1 if outside, and 0 if indefinite
        """

        # quick check, point must be inside box
        if reduce(lambda x, y: x or y,
                  [(point[i] < self.xmins[i] - minDistance) or
                   (point[i] > self.xmaxs[i] + minDistance)
                   for i in range(self.ndims)]):
            # outside box
            return -1

        # any direction will do but things will run faster if the direction
        # points to the nearest domain box face (fewer intersections to
        #  compute)
        self.computeOptimalDirection(point)

        # set the first column in our matrix (independent of the poly)
        self.mat[:, 0] = self.direction

        # compute the number of intersections between the ray and the polys
        numberOfIntersections = 0
        for poly in self.polys:

            isTriangle = 1
            numPoints = len(poly)
            if numPoints == 4:
                isTriangle = 0

            # compute the overlap betwen the ray and face boxes
            if not self.areBoxesOverlapping(point, poly, minDistance):
                # intersection is impossible, don't bother...
                continue

            lmbda, xis = self.computeIntersection(point, poly)

            if lmbda > -self.eps:

                # the direction is ok

                sums = 0.0
                rayIntersects = True
                for i in range(len(xis)):
                    rayIntersects &= (xis[i] >= 0.0)
                    rayIntersects &= (xis[i] < 1.0 - isTriangle*sums)
                    sums += xis[i]

                if rayIntersects:

                    target = lmbda * self.direction
                    distancePointToFace = math.sqrt(numpy.dot(target, target))
                    if distancePointToFace <= minDistance:
                        # marginal, cannot say
                        return 0
                    else:
                        numberOfIntersections += 1

        # even number is outside, odd number means inside
        return 2*(numberOfIntersections % 2) - 1

    def areBoxesOverlapping(self, point, poly, minDistance):
        """
        Determine if there is overlap between the ray and the poly min/max
        coordinates
        @param point starting point of the ray
        @param poly the face
        @param minDistance the distance below which we cannot say for sure
               whether there is an overlap or not
        """

        for i in range(self.ndims):

            xminFace = min([self.points[j][i] for j in poly])
            xmaxFace = max([self.points[j][i] for j in poly])

            if (self.xmaxRay[i] < xminFace - minDistance) or \
               (self.xminRay[i] > xmaxFace + minDistance):
                # no overlap possible
                return False
        # there is a strong possibility of overlap
        return True

    def computeOptimalDirection(self, point):
        """
        Compute the direction of the ray to the nearest
        domain box face and perturb that direction slightly 
        to minimize the probability that the ray hits a  node
        @param point starting point of the ray
        @note this will update self.direction
        """
        # iterate over the faces of the box then find the minimum
        # distance between the point and the face.
        minDistance = float('inf')

        # high/low side of the box
        for pm in (-1, 1):
            pls = (1 + pm)/2.  # 0 on the low side, 1 on the high side
            mns = (1 - pm)/2.  # 1 on the low side, 0 on the high side

            # iterate over the axes
            for axis in range(self.ndims):
                # the normal vector contains very small values in place of
                # zeros in order to avoid issues with rays hitting the exact
                # location of a node
                normal = numpy.array([(i + 1)*1.23456789e-8 for i
                                      in range(self.ndims)])
                normal[axis] = pm

                distance = pls*(self.xmaxs[axis] - point[axis]) + \
                    mns*(point[axis] - self.xmins[axis])
                if distance < minDistance:
                    self.direction = normal * max(distance, minDistance)
                    minDistance = distance

        for i in range(self.ndims):
            p, d = point[i], self.direction[i]
            self.xminRay[i] = min(p, p + d)
            self.xmaxRay[i] = max(p, p + d)

    def computeIntersection(self, point, poly):
        """
        Compute the intersection with either a triangle or a quadrilateral
        @param point starting point of the ray
        @param poly either a triangle or quad
        @return lmbda (linear parametric coordinate), xis (triangle parametric coordinates)
        """
        ip0 = poly[0]
        self.b = self.points[ip0] - point
        self.mat[:, 1] = self.points[ip0] - self.points[poly[1]]
        self.mat[:, 2] = self.points[ip0] - self.points[poly[-1]]
        solution = numpy.linalg.solve(self.mat, self.b)
        return solution[0], solution[1:]


##############################################################################
def test():

    from icqsol.shapes.ShapeManager import ShapeManager

    shape_mgr = ShapeManager('vtk', 'POLYDATA')
    shp = shape_mgr.createShape('sphere', origin=(0., 0., 0.),
                                radius=1.0, n_theta=6, n_phi=3)

    inside = Inside(shp)

    pt = numpy.array([0., 0., 0.])
    assert(inside.isInside(pt, 0.) == 1)

    pt = numpy.array([1.01, 0., 0.])
    assert(inside.isInside(pt, 0.) == -1)

    pt = numpy.array([0., 0.9999999999, 0.])
    assert(inside.isInside(pt, 0.0) == 1)
    assert(inside.isInside(pt, 0.01) == 0)

    pt = numpy.array([0., 1.0000000001, 0.])
    assert(inside.isInside(pt, 0.) == -1)
    assert(inside.isInside(pt, 0.01) == 0)

if __name__ == '__main__':
    test()
