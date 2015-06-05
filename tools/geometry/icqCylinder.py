#!/usr/bin/env python

import vtk

class Cylinder(vtk.vtkImplicitBoolean):

  def __init__(self, radius, origin, length):
    """
    Constructor
    @param radius radius
    @param origin center of the cylinder in the x, y directions
    @param length length of the cylinder in the z direction
    """

    # infinite cylinder in the y direction
    self.cyl = vtk.vtkCylinder()
    self.cyl.SetRadius(radius)
    self.cyl.SetCenter(origin)

    # rotate the cylinder so the axis is in z
    self.cylTransform = vtk.vtkGeneralTransform()
    self.cylTransform.RotateX(90.0)
    self.cyl.SetTransform(self.cylTransform)

    # cut the cylinder
    self.planeLo = vtk.vtkPlane()
    self.planeLo.SetOrigin(origin[0], origin[1], origin[2] - 0.5*length)
    self.planeLo.SetNormal(0., 0., -1.)
    self.planeHi = vtk.vtkPlane()
    self.planeHi.SetOrigin(origin[0], origin[1], origin[2] + 0.5*length)
    self.planeHi.SetNormal(0., 0., 1.)

    # combine
    self.SetOperationTypeToIntersection()
    self.AddFunction(self.cyl)
    self.AddFunction(self.planeLo)
    self.AddFunction(self.planeHi)

    self.loBounds = numpy.array([origin[0] - radius, origin[1] - radius, origin[2] - 0.5*length])
    self.hiBounds = numpy.array([origin[0] + radius, origin[1] + radius, origin[2] + 0.5*length])

  def getBounds(self): 
    """
    Get min/max bounds
    @return low bound, hi bound
    """
    return self.loBounds, self.hiBounds
