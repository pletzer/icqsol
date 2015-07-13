#!/usr/bin/env  python

import  vtk
import  numpy

from  icqShape  import  Shape
from csg.core import CSG
from csg.geom import Vector

class  Cone(Shape):

    def  __init__(self,  radius,  origin,  lengths,
                        n_theta=32):
        """
        Constructor
        @param  radius  radius
        @param  origin  location  of  the  focal  point
        @param  lengths  lengths  of  the  cone  
        @param  n_theta  number  of  theta  cells
        """

        ori = Vector(origin[0], origin[1], origin[2])
        end = Vector(origin[0] + lengths[0],
                     origin[1] + lengths[1],
                     origin[2] + lengths[2])
        shp = CSG.cone(start = ori, 
                       end = end,
                       radius = radius,
                       slices = n_theta)

        Shape.__init__(self, csg=shp)

################################################################################
def  test():

    con  =  Cone(radius=1.0,  origin=(0.,  0.,  0.),  lengths=[1., 0., 0.],
                 n_theta=8,)
    con.save('con.vtk',  file_format='vtk',  file_type='ascii')
    con.show()

if  __name__  ==  '__main__':
    test()

