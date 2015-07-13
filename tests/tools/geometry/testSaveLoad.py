"""
Test save/load methods
@author pletzer@psu.edu
"""

from icqsol.tools.geometry.icqShape import Shape
from icqsol.tools.geometry.icqBox import Box
from icqsol.tools.geometry.icqSphere import Sphere

shp = Box(origin = [0., 0., 0.], lengths = [1., 1., 1.],)
shp.save(file_name='testSaveLoad.vtk', file_format='vtk', file_type='binary')

shp2 = Shape()
shp2.load('testSaveLoad.vtk')

