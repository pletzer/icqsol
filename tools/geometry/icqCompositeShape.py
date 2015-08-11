#!/usr/bin/env python

import re
from icqShape import Shape


class CompositeShapeBuilder( Shape ):

    def __init__( self, shape_tuples=None ):
        """
        Constructor
        @param shape_tuples list of tuples (expr_var, shape)
        """
        self.shape_tuples = shape_tuples

    def compose( self, expression ):
        if expression is None or self.shape_tupels is None:
            return None
        expr = expression
        for expr_var, shape in self.shape_tuples:
            # Return the string obtained by replacing the leftmost
            # non-overlapping occurrences of expr_var in expr by
            # the replacement shape.
            expr = re.sub( expr_var, shape, expr )
        return eval( expr )
