#!/usr/bin/env python

from icqShape import Shape


def CompositeShape(shape_tuples=[], expression=''):
    """
    @param shape_tuples list of (variable_name, shape) pairs
    @param expression expression involving +, -, and * operations. The name 
           of the variables must match that used i
    """
    for i in range(len(shape_tuples)):
        exec(shape_tuples[i][0] + ' = shape_tuples[i][1]')
    return eval(expression)

###############################################################################

def test():
    from icqSphere import Sphere
    s1 = Sphere(radius=1, origin=(0., 0., 0.))
    s2 = Sphere(radius=1.2, origin=(0.8, 0., 0.))
    s3 = CompositeShape([('s1', s1), ('s2', s2)], 's1 + s2')
    s3.save('s3.vtk', file_format='vtk', file_type='ascii')

if __name__ == '__main__': test()