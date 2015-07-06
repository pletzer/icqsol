from icqsol.tools.geometry.icqSphere import Sphere

"""
Test conversion from a shape to a list of polygons
@author pletzer@psu.edu
"""

sphere = Sphere(origin = [0., 0., 0.], radius = 1.0)
print sphere.toPolygons()

