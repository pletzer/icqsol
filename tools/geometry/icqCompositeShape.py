#!/usr/bin/env python

import re
from icqShape import Shape


class ShapeComposer( Shape ):

    def __init__( self, shape_tuples=None ):
        self.shape_tuples = shape_tuples
        self.csg = None

    def __add__( self, other ):
        """
        Union
        @param other Shape instance
        @return CompositeShape
        """
        self.csg = Shape( self.csg + other.csg )
        return self.csg

    def __sub__( self, other ):
        """
        Removal
        @param other Shape instance
        @return CompositeShape
        """
        self.csg = Shape( self.csg - other.csg )
        return self.csg

    def __mul__( self, other ):
        """
        Intersection
        @param other Shape instance
        @return CompositeShape
        """
        self.csg = Shape( self.csg * other.csg )
        return self.csg

    def compose( self, expression ):
        expr = expression
        for shape_tuple in self.shape_tuples:
            expr_var, shape = shape_tuple
            expr = re.sub( r'%s' % expr_var, shape, expr )
        return eval( expr )
