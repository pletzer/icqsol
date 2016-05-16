#!/usr/bin/env python
"""
@brief A base class for constructing shapes
"""
from __future__ import print_function
import os
import re
import vtk
import numpy
# We need the following to handle expressions received from callers.
from numpy import linspace
from math import sqrt, sin, cos, tan, log, exp, pi, asin, acos, atan, atan2, e

from csg.geom import Vector, Vertex, Polygon, BSPNode
from csg.core import CSG
from icqsol.shapes.icqShape import Box, Cone, Cylinder, Sphere
from icqsol.shapes.icqShape import DEFAULTS, CompositeShape
from icqsol.color.icqColorMap import ColorMap
from icqsol.shapes.icqRefineSurface import RefineSurface
from icqsol.shapes.icqCoarsenSurface import CoarsenSurface

LOCATIONS = ['POINT', 'CELL']
VTK_DATASET_TYPES = ['STRUCTURED_GRID', 'POLYDATA', 'UNSTRUCTURED_GRID']
FILE_FORMATS = ['ply', 'vtk']


class ShapeManager(object):

    def __init__(self, file_format=None, vtk_dataset_type=None):
        """
        This class incorporates features from 2 primary classes: CSG and
        VTK.  It uses CSG to assemble shapes, which has no knowledge of
        fields.  VTK data objects are grids with fields typically attached
        to them, but not always.

        If file_format is 'vtk', this class's methods work with the POLYDATA
        vtk_dataset_type, so other types are converted as a first step.
        """
        self.file_format = file_format
        self.vtk_dataset_type = vtk_dataset_type
        self.vtk_geometry_filter = None
        if self.file_format is None:
            self.reader = None
            self.writer = None
        elif self.file_format == 'vtk':
            assert self.vtk_dataset_type in VTK_DATASET_TYPES, \
                "Invalid vtk_dataset_type %s" % str(self.vtk_dataset_type)
            self.setVtkGeometryFilter()
            self.setReader(self.file_format, self.vtk_dataset_type)
            self.setWriter(self.file_format, self.vtk_dataset_type)
        else:
            # We have a PLY file
            self.setReader(self.file_format)
            self.setWriter(self.file_format)

    def createShape(self, type, origin=None, lengths=None, radius=None,
                    angle=None, n_theta=None, n_phi=None):
        """
        Create a primitive shape which can be one of box, cone, cylinder or
        sphere.
        @param type the type of shape: box, cone, cylinder or sphere
        @param origin an (x,y,z) tuple consisting of float origin coordinates
        @param lengths (optional) float lengths in the (x,y,z) directions
        @param radius (optional) float radius
        @param angle (optional) float angle
        @param n_theta (optional) number of longitudes (if applicable)
        @param n_phi (optional) number of latitudes (if applicable)
        """
        # Set defaults if necessary.
        if origin is None:
            origin = DEFAULTS['origin']
        if lengths is None:
            lengths = DEFAULTS['lengths']
        if radius is None:
            radius = DEFAULTS['radius']
        if angle is None:
            angle = DEFAULTS['angle']
        if n_theta is None:
            n_theta = DEFAULTS['n_theta']
        if n_phi is None:
            n_phi = DEFAULTS['n_phi']
        # Create the specified shape.
        if type == 'box':
            return Box(origin, lengths)
        if type == 'cone':
            return Cone(radius, origin, lengths, n_theta)
        if type == 'cylinder':
            return Cylinder(radius, origin, lengths, n_theta)
        if type == 'sphere':
            return Sphere(radius, origin, n_theta, n_phi)
        return None

    def addTextureToVtkPolyData(self, vtk_poly_data,
                                texture_file,
                                max_edge_length=float('inf'),
                                texture_file_format=''):
        """
        Add texture to a vtkPolyData object
        @param vtk_poly_data, VTKPolyData object
        @param texture_file texture file (jpg or png)
        @param max_edge_length maximum edge length, refine if need be
        @param texture_file_format texture file format (jpg or png)
        @return new vtkPolyDataObject
        """

        # Select the reader. The user can either explicitly set the format
        # or the reader can be selected implicitly from the file suffix.
        reader = None
        if texture_file_format.lower().find('jp') >= 0 \
             or texture_file.lower().find('.jp') > 0:
                reader = vtk.vtkJPEGReader()
        elif texture_file_format.lower().find('png') >= 0 \
             or texture_file.lower().find('.png') > 0:
                reader = vtk.vtkPNGReader()
        if not reader:
            msg = 'Could not choose reader from file format {0} or file name {1}'
            raise NotImplementedError(
                        msg.format(texture_file_format, texture_file)
                                      )

        reader.SetFileName(texture_file)
        reader.Update()

        imageData = reader.GetOutput()

        # number of pixels
        n0, n1, one = imageData.GetDimensions()

        # box bounds
        xmin, xmax, ymin, ymax, zmin, zmax = vtk_poly_data.GetBounds()

        # mid position of the box
        xmid = (xmin + xmax)/2.
        ymid = (ymin + ymax)/2.
        zmid = (zmin + zmax)/2.
        xlen = xmax - xmin
        ylen = ymax - ymin
        zlen = zmax - zmin
        midPos = numpy.array([xmid, ymid, zmid])
        
        # extrude a little to ensure that the projection box
        # lie beyond the object
        xlen *= 1.01
        ylen *= 1.01
        zlen *= 1.01        
        xmin = xmid - 0.5*xlen
        xmax = xmid + 0.5*xlen
        ymin = ymid - 0.5*ylen
        ymax = ymid + 0.5*ylen
        zmin = zmid - 0.5*zlen
        zmax = zmid + 0.5*zlen

        # face normals of the box
        faceUnitVecs = [numpy.array([+1., 0., 0.]),
                        numpy.array([0., +1., 0.]),
                        numpy.array([0., 0., +1.]),
                        numpy.array([-1., 0., 0.]),
                        numpy.array([0., -1., 0.]),
                        numpy.array([0., 0., -1.])]

        # the index that does not vary (index where the above
        # normals are non-zero)
        constIndex = [0, 1, 2, 0, 1, 2]

        # the parametric u, v unit vector tangential to each face. u cross v
        # gives the normal
        uvVecs = [(numpy.array([0., +1., 0.]), numpy.array([0., 0., +1.])),
                  (numpy.array([-1., 0., 0.]), numpy.array([0., 0., +1.])),
                  (numpy.array([-1., 0., 0.]), numpy.array([0., -1., 0.])),
                  (numpy.array([0., 0., -1.]), numpy.array([0., -1., 0.])),
                  (numpy.array([0., 0., -1.]), numpy.array([+1., 0., 0.])),
                  (numpy.array([0., +1., 0.]), numpy.array([+1., 0., 0.]))]
                  
        uvLengths = [(ylen, zlen),
                     (xlen, zlen),
                     (xlen, ylen),
                     (zlen, ylen),
                     (zlen, xlen),
                     (ylen, xlen)]

        # the low end start points for each face
        startPoints = [numpy.array([xmax, ymin, zmin]),
                       numpy.array([xmax, ymax, zmin]),
                       numpy.array([xmax, ymax, zmax]),
                       numpy.array([xmin, ymax, zmax]),
                       numpy.array([xmin, ymin, zmax]),
                       numpy.array([xmin, ymin, zmin])]

        # d0 is the number of local indices on each tile (along u and v)
        d = min(n0 // 4, n1 // 3)

        def getImageIndices(xyz):
            """
            Map a position on the object to a set of two indices
            which can be used to retrieve the color from the texture
            file
            @param xyz a point
            @return i0, i1 indices, 0 <= i0 < n0 and 0 <= i1 < n1
            """

            # direction of the ray
            direction = xyz - midPos
            # make sure it is not zero
            direction += 1.e-10*numpy.array([1.234567, 2.3456789, 3.45678])

            # find  the face that is most aligned to the ray
            faceIndex = numpy.argmax([fUnit.dot(direction) for
                                      fUnit in faceUnitVecs])

            uVec = uvVecs[faceIndex][0]
            vVec = uvVecs[faceIndex][1]
            startPos = startPoints[faceIndex]

            # find the intersection point of the ray with the face
            ci = constIndex[faceIndex]
            lambd = (startPos[ci] - midPos[ci])/direction[ci]
            projectedPoint = midPos + lambd*direction

            # global tile indices, the tiles are arranged in staircase fashion
            tileI = (faceIndex + 1) // 2
            tileJ = faceIndex // 2

            # local indices on the tile
            j0 = int(d * uVec.dot(projectedPoint - startPos) / uvLengths[faceIndex][0])
            j1 = int(d * vVec.dot(projectedPoint - startPos) / uvLengths[faceIndex][1])

            # the image indices
            i0, i1 = tileI*d + j0, tileJ*d + j1

            return i0, i1

        # Refine if need be.
        pdata = self.refineVtkPolyData(vtk_poly_data,
                                       max_edge_length=max_edge_length)

        # set the colors from the texture file
        rgbArray = vtk.vtkUnsignedCharArray()
        rgbArray.SetNumberOfComponents(3)
        numPoints = pdata.GetNumberOfPoints()
        rgbArray.SetNumberOfTuples(numPoints)
        rgbArray.SetName('Colors')
        point_data_array = imageData.GetPointData().GetArray(0)
        for i in range(numPoints):
            xyz = pdata.GetPoint(i)
            i0, i1 = getImageIndices(xyz)
            rgba = point_data_array.GetTuple(i0 + n0*i1)
            rgbArray.SetTuple(i, rgba[:3])

        pdata.GetPointData().SetScalars(rgbArray)

        return pdata

    def getFieldCentering(self, vtk_poly_data, field_name):
        """
        Get array centering
        @param vtk_poly_data, vtkPolyData instance
        @param field_name name of the field to integrate
        return True if array is nodal, False if cell centered        
        """
        isPoint = False
        array = vtk_poly_data.GetPointData().GetScalars(field_name)
        if array is None:
            array = vtk_poly_data.GetCellData().GetScalars(field_name)
        else:
            isPoint = True
        # Bail out if field was not found
        if array is None:
            raise NotImplementedError(
              'Could not find field "{0}"!'.format(field_name))
        return isPoint

    def getFieldRange(self, vtk_poly_data, field_name,
                      field_component=0):
        """
        Get the min/max field values
        @param vtk_poly_data, vtkPolyData instance
        @param field_name name of the field to integrate
        @param field_component field component
        @return min, max values        
        """
        # Get array and centering.
        isPoint = self.getFieldCentering(vtk_poly_data, field_name)
        array = None
        if isPoint:
            array = vtk_poly_data.GetPointData().GetScalars(field_name)
        else:
            array = vtk_poly_data.GetCellData().GetScalars(field_name)
        numComps = array.GetNumberOfComponents()
        assert field_component < numComps, \
            "Field component {0} must be < {1}".format(field_component, numComps)
        return array.GetRange(field_component)
        
    def integrateSurfaceField(self, vtk_poly_data, field_name,
                              field_component=0):
        """
        Surface integral of a field (point or cell)
        @param vtk_poly_data, vtkPolyData instance
        @param field_name name of the field to integrate
        @param field_component field component
        @return surface integral
        """
        res = 0
        # Determine if the field is point or cell centered.
        isPoint = self.getFieldCentering(vtk_poly_data, field_name)
        if isPoint:
            array = vtk_poly_data.GetPointData().GetScalars(field_name)
        else:
            array = vtk_poly_data.GetCellData().GetScalars(field_name)
        numComps = array.GetNumberOfComponents()
        assert field_component < numComps, \
            "Field component {0} must be < {1}".format(field_component, numComps)
        
        # Iterate over all the polys.
        points = vtk_poly_data.GetPoints()
        cells = vtk_poly_data.GetPolys()
        numCells = cells.GetNumberOfCells()
        ptIds = vtk.vtkIdList()
        cells.InitTraversal()
        if isPoint:
            # Point array. Average the nodal field to get the 
            # cell centered value
            for i in range(numCells):
                cell = cells.GetNextCell(ptIds)
                npts = ptIds.GetNumberOfIds()
                if npts < 3: 
                    continue
                ia = ptIds.GetId(0)
                pa = numpy.array(points.GetPoint(ia))
                fa = array.GetComponent(ia, field_component)
                for j in range(1, npts - 1):
                    ib, ic = ptIds.GetId(j), ptIds.GetId(j + 1)
                    pb = numpy.array(points.GetPoint(ib))
                    pc = numpy.array(points.GetPoint(ic))
                    pb -= pa
                    pc -= pa
                    areaTimesTwo = numpy.linalg.norm(numpy.cross(pb, pc))
                    fb = array.GetComponent(ib, field_component)
                    fc = array.GetComponent(ic, field_component)
                    res += areaTimesTwo * (fa + fb + fc) / 6.0
        else:
            # Cell centered field
            for i in range(numCells):
                cell = cells.GetNextCell(ptIds)
                npts = ptIds.GetNumberOfIds()
                if npts < 3:
                    continue
                ia = ptIds.GetId(0)
                pa = numpy.array(points.GetPoint(ia))
                f = array.GetComponent(i, field_component)
                for j in range(1, npts - 1):
                    ib, ic = ptIds.GetId(j), ptIds.GetId(j + 1)
                    pb = numpy.array(points.GetPoint(ib))
                    pc = numpy.array(points.GetPoint(ic))
                    pb -= pa
                    pc -= pa
                    areaTimesTwo = numpy.linalg.norm(numpy.cross(pb, pc))
                    res += areaTimesTwo * f / 2.0
                    
        return res

    def addSurfaceFieldFromExpressionToVtkPolyData(self, vtk_poly_data,
                                                   field_name,
                                                   expression,
                                                   time_points,
                                                   max_edge_length=float('inf'),
                                                   location='POINT'):
        """
        Add a surface field to a shape using an expression consisting of
        legal variables x,y,z (shape point coordinates) and t (time).
        @param vtk_poly_data, VTKPolyData converted from shape
        @param field_name, the name of the surface field
        @param expression, expression of variables x, y, z, and t
        @param time_points, list of floating point values defining
               snapshots in a time sequence
        @param max_edge_length maximum edge length, refine if need be
        @param location location of field within cell, either 'POINT' or 'CELL'
        @return vtkPolyData instance
        """
        # Sometimes the received expression has been mangled.
        for op in [('__gt__', '>'), ('__lt__', '<'),
                   ('__ge__', '>='), ('__le__', '<='),
                   ('__eq__', '=='), ('__ne__', '!=')]:
            expression = re.sub(op[0], op[1], expression)
        
        # Refine if need be.
        pdata = self.refineVtkPolyData(vtk_poly_data,
                                       max_edge_length=max_edge_length)
        # The field name cannot have spaces.
        valid_field_name = field_name.replace(' ', '_')
        # Get the points from the shape.
        points = pdata.GetPoints()
        num_points = points.GetNumberOfPoints()
        # Define the data.
        data = vtk.vtkDoubleArray()
        data.SetName(valid_field_name)
        # Handle time points.
        num_time_points = len(time_points)
        # Update the data.
        data.SetNumberOfComponents(num_time_points)
        if location.upper() == 'POINT':
            data.SetNumberOfTuples(num_points)
            # Set the surface field values.
            for i in range(num_points):
                x, y, z = points.GetPoint(i)
                for j in range(num_time_points):
                    t = time_points[j]
                    field_value = eval(expression)
                    data.SetComponent(i, j, field_value)
            # Add the field.
            pdata.GetPointData().AddArray(data)
        elif location.upper() == 'CELL':
            num_cells = pdata.GetNumberOfCells()
            data.SetNumberOfTuples(num_cells)
            # Set the surface field values.
            ptIds = vtk.vtkIdList()
            for i in range(num_cells):
                pdata.GetCellPoints(i, ptIds)
                avrgPosition = numpy.zeros((3,), numpy.float64)
                num_pts = ptIds.GetNumberOfIds()
                for k in range(num_pts):
                    avrgPosition += numpy.array(
                        points.GetPoint(ptIds.GetId(k)))
                avrgPosition /= num_pts
                x, y, z = avrgPosition
                for j in range(num_time_points):
                    t = time_points[j]
                    field_value = eval(expression)
                    data.SetComponent(i, j, field_value)
            # Add the field.
            pdata.GetCellData().AddArray(data)

        return pdata

    def addSurfaceFieldFromExpressionToShape(self, shape, field_name,
                                             expression,
                                             time_points,
                                             max_edge_length=float('inf'),
                                             location='POINT'):
        """
        Add a surface field to a shape using an expression consisting of
        legal variables x,y,z (shape point coordinates) and t (time).
        @param shape
        @param field_name the name of the surface field
        @param expression expression consisting containing variables
                          x, y, z, and t
        @param time_points list of floating point values defining
               snapshots in a time sequence
        @param max_edge_length max edge length for refinement
        @param location location of field within cell, either 'POINT' or 'CELL'
                        (capitalization does not matter)
        @return pdata VTKPolyData converted from shape with added surface field
        """
        # Get the points from the shape.
        pdataInput = self.shapeToVTKPolyData(shape)

        # Refine if need be.
        pdata = self.refineVtkPolyData(pdataInput,
                                       max_edge_length=max_edge_length)
        return self.addSurfaceFieldFromExpressionToVtkPolyData(pdata,
                                                               field_name,
                                                               expression,
                                                               time_points,
                                                               max_edge_length,
                                                               location)

    def colorSurfaceField(self, vtk_poly_data, color_map,
                          field_name='', field_component=0):
        """
        Color a selected surface field of a shape.
        @param vtk_poly_data data defining the shape with
                             surface field information
        @param color_map, name of the color map to use
        @param field_name, the name of the surface field to color
        @param field_component field component
        @return colored_vtk_poly_data, color applied to vtk_poly_data
        """
        # Copy the received vtk_poly_data, creating another
        # vtkPolyData object with the same points and cells.
        vtk_poly_data_copy = vtk.vtkPolyData()
        vtk_poly_data_copy.SetPoints(vtk_poly_data.GetPoints())
        vtk_poly_data_copy.SetPolys(vtk_poly_data.GetPolys())
        # Get information from the data.
        isPoint = self.getFieldCentering(vtk_poly_data, field_name)
        if isPoint:
            array = vtk_poly_data.GetPointData().GetScalars(field_name)
        else:
            array = vtk_poly_data.GetCellData().GetScalars(field_name)
        numComps = array.GetNumberOfComponents()
        assert field_component < numComps, \
            "Field component should be < {0}".format(numComps)
        # Get the min/max field values.
        fmin, fmax = array.GetRange()
        # Prepare for coloring.
        rgbArray = vtk.vtkUnsignedCharArray()
        rgbArray.SetNumberOfComponents(3)
        numElements = array.GetNumberOfTuples()
        rgbArray.SetNumberOfTuples(numElements)
        rgbArray.SetName('Colors')
        # Build the color map.
        colorMap = ColorMap(fmin, fmax)
        colorMethod = eval('colorMap.' + color_map)
        # Set the colors.
        for i in range(numElements):
            f = array.GetComponent(i, field_component)
            rgbArray.SetTuple(i, colorMethod(f))
        # Attach the array.
        if isPoint:
            vtk_poly_data_copy.GetPointData().SetScalars(rgbArray)
        else:
            vtk_poly_data_copy.GetCellData().SetScalars(rgbArray)
        return vtk_poly_data_copy

    def cleanPolygon(self, verts, min_triangle_area):
        """
        Remove vertices that are degenerate
        @param verts list of Vertex instances (input and output)
        @param min_triangle_area minimum triangle area
        """
        numVerts = len(verts)
        if numVerts < 3:
            # Nothing to do.
            return
        pa = verts[0].pos
        area = verts[1].pos.minus(pa).cross(verts[2].pos.minus(pa)).length()
        if area > min_triangle_area:
            # We can compute the normal.
            return
        else:
            # Remove vertex and recursively call this method.
            del verts[1]
            self.cleanPolygon(verts, min_triangle_area)

    def refineShape(self, shape, refine=1):
        """
        Refine a shape by inserting one vertex to each edge and one
        vertex to the center of the polygon
        @param shape shape
        @param refine number of refinements
        @return new shape
        """
        s = shape.clone()
        for i in range(refine):
            s = s.refine()
        return s

    def refineVtkPolyData(self, polydata, max_edge_length):
        """
        Refine a vtkPolyData object by adding points along cell edges
        @param polydata vtkPolyData instance
        @param max_edge_length maximum edge length, edges smaller than
                               this value will not be segmented
        @return vtkPolyData instance
        """
        rs = RefineSurface(polydata)
        rs.refine(max_edge_length=max_edge_length)
        return rs.getVtkPolyData()

    def coarsenVtkPolyData(self, polydata, min_cell_area):
        """
        Coarsen a vtkPolyData object by coalescing cells
        @param polydata vtkPolyData instance
        @param min_cell_area minimum cell area beyond which cells
               will be merged
        @return vtkPolyData instance
        """
        cs = CoarsenSurface(polydata)
        cs.coarsen(min_cell_area=min_cell_area)
        return cs.getVtkPolyData()

    def cloneShape(self, shape):
        """
        Clone shape
        @param shape shape
        @return new shape
        """
        return shape.clone()

    def composeShapes(self, shape_tuples=[], expression=''):
        """
        Compose shapes into a more complex shape.
        @param shape_tuples list of (variable_name, shape) pairs
        @param expression expression involving +, -, and * operations.
        @return new shape
        """
        return CompositeShape(shape_tuples, expression)

    def getBoundarySurfaceInsideShape(self, shape, other):
        """
        Return the portion of the surface that is inside another shape
        @param shape
        @param other other shape
        @return shape
        """
        a = BSPNode(shape.clone().polygons)
        b = BSPNode(other.clone().polygons)
        b.invert()
        a.clipTo(b)
        return CSG.fromPolygons(a.allPolygons())

    def loadAsVtkData(self, file_name):
        """
        Load a subclass of a vtkData object from a file
        where a subclass is one of vtkStructuredPoints(),
        vtkStructuredGrid(), vtkPolyData(), vtkRectilinearGrid()
        or vtkUnstructuredGrid()
        @param file_name suffix should be .ply or .vtk
        @return a subclass of a VtkData object
        """
        if not os.path.exists(file_name):
            raise IOError('File {0} not found'.format(file_name))
        self.reader.SetFileName(file_name)
        self.reader.Update()
        return self.reader.GetOutput()

    def loadAsVtkPolyData(self, file_name):
        """
        Load a vtkPolyData object from a file.
        @param file_name, suffix should be .ply or .vtk
        @return a vtkPolyData object
        """
        vtk_data = self.loadAsVtkData(file_name)
        if self.file_format == 'vtk':
            if self.vtk_dataset_type == 'POLYDATA':
                vtk_poly_data = vtk_data
            else:
                vtk_poly_data = self.convertToPolyData(vtk_data)
        else:
            # We have a PLY file.
            vtk_poly_data = vtk_data
        return vtk_poly_data

    def loadAsShape(self, file_name):
        """
        Load an in-memory shape from a file
        @param file_name suffix should be .ply or .vtk
        @return a Shape object
        """
        vtk_poly_data = self.loadAsVtkPolyData(file_name)
        return self.shapeFromVTKPolyData(vtk_poly_data)

    def rotateVtkPolyData(self, pdata, axis=(1., 0., 0.), angleDeg=0.0):
        """
        Rotate vtkPolyData object along given axis
        @param pdata vtkPolyData instance (modified on output)
        @param axis rotation axis
        @param angleDeg angle in degrees
        """
        transform = vtk.vtkTransform()
        transform.RotateWXYZ(angleDeg, axis[0], axis[1], axis[2])
        transformFilter = vtk.vtkTransformPolyDataFilter()
        transformFilter.SetTransform(transform)
        if vtk.VTK_MAJOR_VERSION >= 6:
            transformFilter.SetInputData(pdata)
        else:
            transformFilter.SetInput(pdata)
        transformFilter.Update()
        pdata.DeepCopy(transformFilter.GetOutput())

    def translateVtkPolyData(self, pdata, displ=(0., 0., 0.)):
        """
        Translate vtkPolyData
        @param pdata vtkPolyData instance (modified on output)
        @param displ displacement vector
        """
        transform = vtk.vtkTransform()
        transform.Translate(displ[0], displ[1], displ[2])
        transformFilter = vtk.vtkTransformPolyDataFilter()
        transformFilter.SetTransform(transform)
        if vtk.VTK_MAJOR_VERSION >= 6:
            transformFilter.SetInputData(pdata)
        else:
            transformFilter.SetInput(pdata)
        transformFilter.Update()
        pdata.DeepCopy(transformFilter.GetOutput())

    def scaleVtkPolyData(self, pdata, factors=(1., 1., 1.)):
        """
        Translate vtkPolyData
        @param pdata vtkPolyData instance (modified on output)
        @param amplification vector
        """
        transform = vtk.vtkTransform()
        transform.Scale(factors[0], factors[1], factors[2])
        transformFilter = vtk.vtkTransformPolyDataFilter()
        transformFilter.SetTransform(transform)
        if vtk.VTK_MAJOR_VERSION >= 6:
            transformFilter.SetInputData(pdata)
        else:
            transformFilter.SetInput(pdata)
        transformFilter.Update()
        pdata.DeepCopy(transformFilter.GetOutput())

    def rotateShape(self, shape, axis=(1., 0., 0.), angleDeg=0.0):
        """
        Rotate along axis
        @param shape
        @param axis rotation axis
        @param angleDeg angle in degrees
        """
        return shape.rotate(axis, angleDeg)

    def saveShape(self, shape, file_name, file_type, normals=True):
        """
        Save the shape to a file
        @param shape for saving
        @param file_name file name
        @param file_type either 'ascii' or 'binary'
        @param normals resolve features (corners) and save
                       normal vectors if True
        """
        vtk_poly_data = self.shapeToVTKPolyData(shape)
        self.saveVtkPolyData(vtk_poly_data, file_name, file_type)

    def saveVtkPolyData(self, vtk_poly_data, file_name,
                        file_type, normals=True):
        """
        Save the vtk_poly_data to a file
        @param vtk_poly_data for saving
        @param file_name file name
        @param file_type either 'ascii' or 'binary'
        @param normals resolve features (corners) and save
                       normal vectors if True
        """
        # compute the vertex normals
        if normals:
            vtk_pdata_save = self.computeVertexNormals(vtk_poly_data)
        else:
            vtk_pdata_save = vtk_poly_data

        self.writer.SetFileName(file_name)
        if file_type.lower() == 'ascii':
            self.writer.SetFileTypeToASCII()
        else:
            self.writer.SetFileTypeToBinary()
        if vtk.VTK_MAJOR_VERSION >= 6:
            self.writer.SetInputData(vtk_pdata_save)
        else:
            self.writer.SetInput(vtk_pdata_save)
        self.writer.Write()
        self.writer.Update()

    def shapeFromPolygons(self, polys):
        return CSG.fromPolygons(polys)

    def shapeToPolygons(self, shape):
        return shape.toPolygons()

    def shapeFromVTKPolyData(self, pdata, min_cell_area=1.e-8):
        """
        Create a shape from a VTK PolyData object
        @param pdata vtkPolyData instance
        @param min_cell_area tolerance for cell areas
        @return shape
        @note field data will get lost
        """
        # Store the cell connectivity as CSG polygons.
        numCells = pdata.GetNumberOfPolys()
        cells = pdata.GetPolys()
        cells.InitTraversal()
        ptIds = vtk.vtkIdList()
        polygons = []
        for i in range(numCells):
            cells.GetNextCell(ptIds)
            npts = ptIds.GetNumberOfIds()
            verts = []
            for j in range(npts):
                pointIndex = ptIds.GetId(j)
                pt = pdata.GetPoint(pointIndex)
                v = Vertex(Vector(pt[0], pt[1], pt[2]))
                verts.append(v)
            self.cleanPolygon(verts, min_cell_area)
            if len(verts) >= 3:
                polygons.append(Polygon(verts))
        # Instantiate the shape.
        return CSG.fromPolygons(polygons)

    def shapeToVTKPolyData(self, shape):
        """
        Convert shape to a VTK polydata object
        """
        verts, polys, count = shape.toVerticesAndPolygons()

        points = vtk.vtkPoints()
        numPoints = len(verts)
        points.SetNumberOfPoints(numPoints)
        for i in range(numPoints):
            points.SetPoint(i, verts[i])

        pdata = vtk.vtkPolyData()
        pdata.SetPoints(points)
        
        # build the connectivity
        numPolys = len(polys)
        pdata.Allocate(numPolys, 1)
        ptIds = vtk.vtkIdList()
        for poly in polys:
            numPolyPts = len(poly)
            ptIds.SetNumberOfIds(numPolyPts)
            for j in range(numPolyPts):
                ptIds.SetId(j, poly[j])
            pdata.InsertNextCell(vtk.VTK_POLYGON, ptIds)

        return pdata

    def computeVertexNormals(self, pdata, min_feature_angle=60.0):
        """
        Compute the vertex normals
        @param pdata vtkPolyData instance
        @param min_feature_angle angle above which a sharp feature
                                 is detected. Additional points
                                 will be added to resolve the feature
        @return new vtkPolyData object
        """
        normalsFilter = vtk.vtkPolyDataNormals()
        normalsFilter.SplittingOn()
        normalsFilter.SetFeatureAngle(min_feature_angle)
        normalsFilter.ComputePointNormalsOn()
        if vtk.VTK_MAJOR_VERSION >= 6:
            normalsFilter.SetInputData(pdata)
        else:
            normalsFilter.SetInput(pdata)
        pdata_new = normalsFilter.GetOutput()
        normalsFilter.Update()
        return pdata_new

    def showShape(self, shape, windowSizeX=600, windowSizeY=400, filename=''):
        """
        Show the boundary surface or write image to file
        @param shape CSG instance
        @param windowSizeX number of pixels in x
        @param windowSizeY number of pixels in y
        @param filename write to a file if this keyword
               is present and a non-empty string
        """
        pdata = self.shapeToVTKPolyData(shape)
        self.showVtkPolyData(pdata, windowSizeX, windowSizeY, filename)

    def showVtkPolyData(self, pdata, windowSizeX=600,
                        windowSizeY=400, filename=''):
        """
        Show the boundary surface or write image to file
        @param pdata vtkPolyData instance
        @param windowSizeX number of pixels in x
        @param windowSizeY number of pixels in y
        @param filename write to a file if this keyword
        is present and a non-empty string
        """
        # Create a rendering window and renderer.
        try:
            ren = vtk.vtkRenderer()
            renWin = vtk.vtkRenderWindow()
            iren = vtk.vtkRenderWindowInteractor()
            camera = vtk.vtkCamera()
            mapper = vtk.vtkPolyDataMapper()
            actor = vtk.vtkActor()
            axes = [vtk.vtkArrowSource(),
                    vtk.vtkArrowSource(),
                    vtk.vtkArrowSource()]
            axesTransf = [vtk.vtkTransform(),
                          vtk.vtkTransform(),
                          vtk.vtkTransform()]
            axesTPD = [vtk.vtkTransformPolyDataFilter(),
                       vtk.vtkTransformPolyDataFilter(),
                       vtk.vtkTransformPolyDataFilter()]
            axesMappers = [vtk.vtkPolyDataMapper(),
                           vtk.vtkPolyDataMapper(),
                           vtk.vtkPolyDataMapper()]
            axesActors = [vtk.vtkActor(), vtk.vtkActor(), vtk.vtkActor()]

            renderLarge = vtk.vtkRenderLargeImage()
        except:
            print('WARNING: Cannot call show method -- likely missing VTK')
            return
        renWin.AddRenderer(ren)
        renWin.SetSize(windowSizeX, windowSizeY)
        # Create a renderwindowinteractor.
        iren.SetRenderWindow(renWin)
        # Camera
        xmin, xmax, ymin, ymax, zmin, zmax = pdata.GetBounds()
        lo = numpy.array([xmin, ymin, zmin])
        hi = numpy.array([xmax, ymax, zmax])
        camera.SetFocalPoint(hi)
        center = 0.5*(lo + hi)
        camera.SetPosition(center + hi - lo)
        camera.Zoom(1.0)
        ren.SetActiveCamera(camera)
        # Mapper.
        if vtk.VTK_MAJOR_VERSION >= 6:
            mapper.SetInputData(pdata)
        else:
            mapper.SetInput(pdata)
        # Actor.
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(1, 1, 1)
        # Add axes.
        axesColrs = [(1., 0., 0.,), (0., 1., 0.,), (0., 0., 1.,)]
        for a in axes:
            a.SetShaftRadius(0.01)
            a.SetTipLength(0.2)
            a.SetTipRadius(0.03)
        for at in axesTransf:
            at.PostMultiply()
        # Rotate the y and z arrows (initially along x).
        axesTransf[1].RotateZ(90.0)
        axesTransf[2].RotateY(-90.0)
        # Scale.
        for i in range(3):
            factor = hi[i] - lo[i]
            scale = [1., 1., 1.]
            scale[i] = factor
            axesTransf[i].Scale(scale)
        # Translate to loBounds.
        for at in axesTransf:
            at.Translate(lo)
        for i in range(3):
            axesTPD[i].SetInputConnection(axes[i].GetOutputPort())
            axesTPD[i].SetTransform(axesTransf[i])
            axesMappers[i].SetInputConnection(axesTPD[i].GetOutputPort())
            axesActors[i].SetMapper(axesMappers[i])
            axesActors[i].GetProperty().SetColor(axesColrs[i])
            ren.AddActor(axesActors[i])
        # Assign actor to the renderer.
        ren.AddActor(actor)
        # Write to file.
        writer = None
        if filename:
            if filename.lower().find('.png') > 0:
                writer = vtk.vtkPNGWriter()
            elif filename.lower().find('.jp') > 0:
                writer = vtk.vtkJPEGWriter()
            elif filename.lower().find('.tiff') > 0:
                writer = vtk.vtkTIFFWriter()
            if writer:
                renderLarge.SetInput(ren)
                renderLarge.SetMagnification(1)
                renderLarge.Update()
                writer.SetFileName(filename)
                writer.SetInputConnection(renderLarge.GetOutputPort())
                writer.Write()
        else:
            # Fire up interactor.
            iren.Initialize()
            renWin.Render()
            iren.Start()

    def translateShape(self, shape, disp=(0., 0., 0.)):
        """
        Translate
        @param shape
        @param disp displacement
        """
        shape.translate(disp)

    def convertToPolyData(self, vtk_data):
        """
        Convert vtk_data to vtk.vtkPolyData().
        @param vtk_data data whose vtk_dataset_type is one of
                         the VTK_DATASET_TYPES besides POLYDATA
        """
        if vtk.VTK_MAJOR_VERSION >= 6:
            self.vtk_geometry_filter.SetInputData(vtk_data)
        else:
            self.vtk_geometry_filter.SetInput(vtk_data)
        self.vtk_geometry_filter.Update()
        vtk_poly_data = self.vtk_geometry_filter.GetOutput()
        return vtk_poly_data

    def getVtkGeometryFilter(self):
        return self.vtk_geometry_filter

    def setVtkGeometryFilter(self):
        """
        Set the vtk_geometry_filter.  We coerce STRUCTURED_GRID and
        UNSTRUCTURED_GRID vtk_dataset_types to be POLYDATA, so we only
        need a vtk.vtkGeometryFilter().
        """
        self.vtk_geometry_filter = vtk.vtkGeometryFilter()

    def chooseReader(self, file_format, vtk_dataset_type=None):
        """
        Return a reader based on file_format and possibly vtk_dataset_type.
        @param vtk_dataset_type, None or one of the VTK_DATASET_TYPES
        """
        # Handle .ply files.
        if file_format == 'ply':
            return vtk.vtkPLYReader()
        # Handle .vtk files.
        if vtk_dataset_type == 'STRUCTURED_GRID':
            return vtk.vtkStructuredGridReader()
        elif vtk_dataset_type == 'POLYDATA':
            return vtk.vtkPolyDataReader()
        elif vtk_dataset_type == 'UNSTRUCTURED_GRID':
            return vtk.vtkUnstructuredGridReader()

    def getReader(self):
        return self.reader

    def setReader(self, file_format=None, vtk_dataset_type=None):
        """
        Set the reader.
        @param file_format, one of vtk or ply - default to self.file_format
        @param vtk_dataset_type, one of the VTK_VTK_DATASET_TYPES
        """
        if file_format is None:
            file_format = self.file_format
        assert file_format in FILE_FORMATS, \
            "Invalid file_format %s" % str(file_format)
        if file_format == 'vtk':
            if vtk_dataset_type is None:
                vtk_dataset_type = self.vtk_dataset_type
            assert vtk_dataset_type in VTK_DATASET_TYPES, \
                "Invalid vtk_dataset_type %s" % str(vtk_dataset_type)
        self.reader = self.chooseReader(file_format, vtk_dataset_type)

    def chooseWriter(self, file_format, vtk_dataset_type):
        """
        Return a writer based on file_format and possibly vtk_dataset_type.
        @param vtk_dataset_type, None or one of the VTK_DATASET_TYPES
        """
        if file_format == 'ply':
            return vtk.vtkPLYWriter()
        # For now we'll just return the POLYDATA writer since methods work
        # only with that vtk_dataset_type.
        return vtk.vtkPolyDataWriter()
        if vtk_dataset_type == 'STRUCTURED_GRID':
            return vtk.vtkStructuredGridWriter()
        elif vtk_dataset_type == 'POLYDATA':
            return vtk.vtkPolyDataWriter()
        elif vtk_dataset_type == 'UNSTRUCTURED_GRID':
            return vtk.vtkUnstructuredGridWriter()

    def getWriter(self):
        return self.writer

    def setWriter(self, file_format=None, vtk_dataset_type=None):
        """
        Set the writer.
        @param file_format, one of vtk or ply - defaults to self.file_format
        @param vtk_dataset_type, None or one of the VTK_VTK_DATASET_TYPES
        """
        if file_format is None:
            file_format = self.file_format
        assert file_format in FILE_FORMATS, \
            "Invalid file_format %s" % str(file_format)
        if file_format == 'vtk':
            if vtk_dataset_type is None:
                vtk_dataset_type = self.vtk_dataset_type
            assert vtk_dataset_type in VTK_DATASET_TYPES, \
                "Invalid vtk_dataset_type %s" % str(vtk_dataset_type)
        self.writer = self.chooseWriter(file_format, vtk_dataset_type)


###############################################################################
def testPrimitiveShapes():
    shape_mgr = ShapeManager(file_format='vtk', vtk_dataset_type='POLYDATA')

    # Box
    box = shape_mgr.createShape('box', origin=(0.,  0.,  0.),
                                lengths=(0.5,  1.,  2.))
    shape_mgr.saveShape(shape=box, file_name='box.vtk', file_type='ascii')
    shape_mgr.showShape(box)

    # Cone
    con = shape_mgr.createShape('cone', radius=1.0, origin=(0.,  0.,  0.),
                                lengths=[1., 0., 0.], n_theta=8)
    shape_mgr.saveShape(shape=con, file_name='con.vtk', file_type='ascii')
    shape_mgr.showShape(con)

    # Cylinder
    cyl = shape_mgr.createShape('cylinder', radius=1.0, origin=(0., 0., 0.),
                                lengths=(1., 0., 0.), n_theta=8)
    shape_mgr.saveShape(shape=cyl, file_name='cyl.vtk', file_type='ascii')
    shape_mgr.showShape(cyl)

    # Sphere
    sph = shape_mgr.createShape('sphere', radius=1.0, origin=(0., 0., 0.),
                                n_theta=8, n_phi=4)
    shape_mgr.saveShape(shape=sph, file_name='sph.vtk', file_type='ascii')
    shape_mgr.showShape(sph)


def testSaveLoad():
    shape_mgr = ShapeManager(file_format='vtk', vtk_dataset_type='POLYDATA')
    s = shape_mgr.createShape('sphere', radius=0.7, origin=(0., 0., 0.),
                              n_theta=8, n_phi=4)
    shape_mgr.saveShape(shape=s, file_name='t.vtk', file_type='ascii')
    s2 = shape_mgr.loadAsShape('t.vtk')
    s3 = shape_mgr.createShape('box', origin=(0.1, 0.2, 0.3),
                               lengths=(1.1, 1.2, 1.3))
    s4 = s2 + s3


def testConstructiveGeometry():
    shape_mgr = ShapeManager(file_format='ply')
    s1 = shape_mgr.createShape('sphere', radius=0.7, origin=(0., 0., 0.))
    s2 = shape_mgr.createShape('sphere', radius=0.2, origin=(0.1, 0.2, 0.3))
    b = shape_mgr.createShape('box', origin=(0.1, 0.2, 0.3),
                              lengths=(1.1, 1.2, 1.3))
    c = shape_mgr.createShape('cylinder', radius=0.5, origin=(0.3, 0.4, 0.5),
                              lengths=(1.0, 0.0, 0.0))
    geom = c*b - s2 - s1
    shape_mgr.saveShape(shape=geom, file_name='geom.ply', file_type='binary')
    geom2 = shape_mgr.loadAsShape('geom.ply')
    shape_mgr.showShape(geom2)


def testShapeComposition():
    shape_mgr = ShapeManager(file_format='vtk', vtk_dataset_type='POLYDATA')
    s1 = shape_mgr.createShape('sphere', radius=1, origin=(0., 0., 0.))
    s2 = shape_mgr.createShape('sphere', radius=1.2, origin=(0.8, 0., 0.))
    s3 = shape_mgr.composeShapes([('s1', s1), ('s2', s2)], 's1 + s2')
    shape_mgr.saveShape(shape=s3, file_name='s3.vtk', file_type='ascii')


if __name__ == '__main__':
    testPrimitiveShapes()
    testSaveLoad()
    testConstructiveGeometry()
    testShapeComposition()
