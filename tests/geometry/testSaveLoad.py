"""
Test save/load methods
@author pletzer@psu.edu
"""

from icqsol.tools.geometry.icqShape import Shape
from icqsol.tools.geometry.icqBox import Box

shp = Box(origin=[0., 0., 0.], lengths=[1., 1., 1.],)
shp.save(file_name='testSaveLoad.vtk', file_format='vtk', file_type='ascii')

shp2 = Shape.load('testSaveLoad.vtk')
shp2.save(file_name='testSaveLoad2.vtk', file_format='vtk', file_type='ascii')
