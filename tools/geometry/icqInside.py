#!/usr/bin/env python

"""
@brief A class for testing if a point is inside a shape
@author pletzer@psu.edu
"""

import numpy


class Inside:

    def __init__(self, shape):
        """
        Constructor
        @param shape instance of Shape
        @param polys connectivty
        """
        points, polys, count = shape.csg.toVerticesAndPolygons()
        # must have at least one point
        self.ndims = len(points[0])

        # find the box corners
        self.xmins = numpy.array([min([p[i] for p in points]) for i in range(self.ndims)])
        self.xmaxs = numpy.array([max([p[i] for p in points]) for i in range(self.ndims)])

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

    def isInside(self, point):
        """
        Determine if a point is inside the shape
        @param point point
        @return +1 if inside, -1 if outside, and 0 if indefinite
        """

        # quick check, point must be inside box
        outsideBox = reduce(lambda x, y: x or y,
                            [(point[i] < self.xmins[i] - self.eps) or
                             (point[i] > self.xmaxs[i] + self.eps)
                             for i in range(self.ndims)])

        if outsideBox:
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

            triangle = 1
            numPoints = len(poly)
            if numPoints == 4:
                triangle = 0

            # compute the overlap betwen the ray and face boxes
            if not self.areBoxesOverlapping(point, poly):
                # intersection is impossible, don't bother...
                continue

            lmbda, xis = self.computeIntersection(point, poly)

            if lmbda > -self.eps:

                # the direction is ok

                sums = 0.0
                rayIntersects = True
                for i in range(len(xis)):
                    rayIntersects &= (xis[i] >= 0.0) #-self.eps)
                    rayIntersects &= (xis[i] < 1 - triangle*sums) # - self.eps)
                    sums += xis[i]

                if rayIntersects:

                    if abs(lmbda) <= self.eps:
                        # marginal, cannot say
                        return 0
                    else:
                        numberOfIntersections += 1

        # even number is outside, odd number means inside
        return 2*(numberOfIntersections % 2) - 1

    def areBoxesOverlapping(self, point, poly):

        for i in range(self.ndims):

            xminFace = min([self.points[j][i] for j in poly])
            xmaxFace = max([self.points[j][i] for j in poly])

            if (self.xmaxRay[i] < xminFace - self.eps) or \
                (self.xminRay[i] > xmaxFace + self.eps):
                    # no overlap
                    return False
        return True

    def computeOptimalDirection(self, point):
        """
        Compute the direction of the ray to the nearest
        domain box face. This will update self.direction
        @param point starting point of the ray
        """
        # iterate over the faces of the box then find the minimum
        # distance between the point and the face.
        minDistance = float('inf')
        for pm in (-1, 1):
            pls = (1 + pm)/2.
            mns = (1 - pm)/2.
            for axis in range(self.ndims):
                # the normal vector contains very small values in place of
                # zeros in order to avoid issues with rays hitting the exact
                # location of a node
                normal = numpy.array([self.eps*(100. + i) for i in range(self.ndims)])
                normal[axis] = pm
                distance = pls*(self.xmaxs[axis] -
                                point[axis]) + mns*(point[axis] -
                                                    self.xmins[axis])
                if distance < minDistance:
                    # expand a little beyond the domain (1.1)
                    self.direction = normal * (1.1 * distance)
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
        """
        ip0 = poly[0]
        self.b = self.points[ip0] - point
        self.mat[:, 1] = self.points[ip0] - self.points[poly[1]]
        self.mat[:, 2] = self.points[ip0] - self.points[poly[-1]]
        solution = numpy.linalg.solve(self.mat, self.b)
        return solution[0], solution[1:]

##############################################################################
def test():

    from icqsol.tools.geometry.icqSphere import Sphere
    shp = Sphere(origin=(0., 0., 0.), radius=1.0, n_theta=6, n_phi=3)
    #shp.debug()

    inside = Inside(shp)

    pt = numpy.array([0., 0., 0.])
    assert(inside.isInside(pt) == 1)

    pt = numpy.array([1.01, 0., 0.])
    assert(inside.isInside(pt) == -1)

    pt = numpy.array([0., 0.9999999999, 0.])
    assert(inside.isInside(pt) == 1)

    pt = numpy.array([0., 1.0000000001, 0.])
    assert(inside.isInside(pt) == -1)



if __name__ == '__main__':
    test()
