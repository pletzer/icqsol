#!/usr/bin/env python

"""
Test intersection operation
@author pletzer@psu.edu
"""

from icqsol.tools.geometry.icqSphere import Sphere
from icqsol.tools.geometry.icqBox import Box

# box-box containment
a = Box(origin=(0., 0., 0.), lengths=(1., 1., 0.2))
b = Box(origin=(0.5, 0.5, 0.), lengths=(1., 1., 0.2))
s = a.getBoundarySurfaceInside(b)

a.save('testInsideOtherShape_BoxBox_a.vtk', 
	   file_format='vtk', file_type='ascii')
b.save('testInsideOtherShape_BoxBox_b.vtk', 
	   file_format='vtk', file_type='ascii')
s.save('testInsideOtherShape_BoxBox_result.vtk', 
	   file_format='vtk', file_type='ascii')

# box - sphere
a = Box(origin=(0., 0., 0.), lengths=(1., 1., 0.2))
b = Sphere(origin=(0.5, 0.5, 0.), radius=0.5)
s = a.getBoundarySurfaceInside(b)

a.save('testInsideOtherShape_BoxSphere_a.vtk', 
	   file_format='vtk', file_type='ascii')
b.save('testInsideOtherShape_BoxSphere_b.vtk', 
	   file_format='vtk', file_type='ascii')
s.save('testInsideOtherShape_BoxSphere_result.vtk', 
	   file_format='vtk', file_type='ascii')


