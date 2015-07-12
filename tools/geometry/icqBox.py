#!/usr/bin/env  python

import  vtk
import  numpy

from  icqShape  import  Shape
from csg.core import CSG

class  Box(Shape):

    def  __init__(self,  origin,  lengths):
        """
        Constructor
        @param  radius  radius
        @param  origin/low  end  of  the  box
        @param  lengths  lengths  in  x,  y,  and  z
        """
        center = [(origin[i] + lengths[i])/2.0 for i in range(len(origin))]
        radius = [le/2. for le in lengths]
        shp = CSG.cube(center = center, 
                      radius = radius)

        self.polygons = shp.polygons

################################################################################
def  test():

    box = Box(origin=(0.,  0.,  0.),  lengths=(0.5,  1.,  2.),)
    box.save('box.vtk',  file_format='vtk',  file_type='ascii')
    box.show()

if  __name__ == '__main__':
    test()

