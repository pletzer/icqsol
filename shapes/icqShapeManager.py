#!/usr/bin/env python
"""
@brief A base class for constructing shapes
"""
import os
import vtk
import numpy
# We need the following to handle expressions received from callers.
from numpy import linspace
from math import sin, cos, tan, log, exp, pi, asin, acos, atan, atan2, e
from csg.geom import Vector, Vertex, Polygon, BSPNode
from csg.core import CSG
from icqShape import Box, Cone, Cylinder, Sphere
from icqShape import DEFAULTS, CompositeShape
from icqsol.color.icqColorMap import ColorMap

VTK_DATASET_TYPES = ['STRUCTURED_GRID', 'POLYDATA', 'UNSTRUCTURED_GRID']
FILE_FORMATS = ['ply', 'vtk']


class ShapeManager(object):

    def __init__(self, file_format=None, vtk_dataset_type=None):
        """
        This class incorporates features from 2 primary classes: CSG and
        VTK.  It uses CSG to assemble shapes, which have no knowledge of
        fields.  VTK data objects are grids with fields typically attached
        to them, but not always.

        If file_format is 'vtk', this class's methods work with the POLYDATA
        vtk_dataset_type, so other types are converted as a first step.
        """
        self.file_format = file_format
        self.vtk_dataset_type = vtk_dataset_type
        if self.file_format is None:
            self.reader = None
            self.writer = None
            self.vtk_geometry_filter = None
        elif self.file_format == 'vtk':
            assert self.vtk_dataset_type in VTK_DATASET_TYPES, "Invalid vtk_dataset_type %s" % str( self.vtk_dataset_type )
            self.vtk_geometry_filter = self.setVtkGeometryFilter(self.vtk_dataset_type)
            self.setReader(self.file_format, self.vtk_dataset_type)
            self.setWriter(self.file_format, self.vtk_dataset_type)
        else:
            # We have a PLY file
            self.setReader(self.file_format)
            self.setWriter(self.file_format)
            self.vtk_geometry_filter = None

    def createShape(self, type, origin=None, lengths=None, radius=None,
                    angle=None, n_theta=None, n_phi=None):
        """
        Create a primitive shape which can be one of box, cone, cylinder or
        sphere.
        @param type, the type of shape: box, cone, cylinder or sphere
        @param origin, an (x,y,z) tuple consisting of float origin coordinates
        @param lengths, (optional) an (x,y,z) tuple consisting of float length coordinates
        @param radius, (optional) float radius
        @param angle, (optional) float angle
        @param n_theta, (optional) float that controls tessellation along longitude
        @param n_phi, (optional) float that controls tessellation along latitude
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

    def addSurfaceFieldFromExpression(self, shape, field_name, expression, time_points):
        """
        Add a surface field to a shape using an expression consisting of
        legal variables x,y,z (shape point coordinates) and t (time).
        @param shape
        @param field_name, the name of the surface field
        @param expression, expression consisting of legal variables x, y, z, and t
        @param time_points, list of floating point values defining
               snapshots in a time sequence
        @return pdata, VTKPolyData converted from shape with added surface field
        """
        # When we attach a surface field, we can no longer use CSG objects
        # so this method returns just the data.
        valid_field_name = field_name.replace( ' ', '_' )
        # Get the points from the shape.
        pdata = self.shapeToVTKPolyData(shape)
        points = pdata.GetPoints()
        num_points = points.GetNumberOfPoints()
        # Define the data.
        data = vtk.vtkDoubleArray()
        data.SetName(valid_field_name)
        # Handle time points.
        num_time_points = len(time_points)
        # Update the data.
        data.SetNumberOfComponents(num_time_points)
        data.SetNumberOfTuples(num_points)
        # Add the surface field.
        for i in range(num_points):
            x, y, z = points.GetPoint(i)
            for j in range(num_time_points):
                t = time_points[ j ]
                field_value = eval(expression)
                data.SetComponent(i, j, field_value)
        pdata.GetPointData().SetScalars(data)
        return pdata

    def colorSurfaceField(self, vtk_poly_data, color_map, field_name=None):
        """
        Color a selected surface field of a shape.
        @param vtk_poly_data, data defining the shape with surface field information
        @param field_name, the name of the surface field to color
        @param color_map, name of the color map to use
        @return colored_vtk_poly_data, color applied to vtk_poly_data
        """
        # Copy the received vtk_poly_data, creating another
        # vtkPolyData object with the same points and cells.
        vtk_poly_data_copy = vtk.vtkPolyData()
        vtk_poly_data_copy.SetPoints(vtk_poly_data.GetPoints())
        vtk_poly_data_copy.SetPolys(vtk_poly_data.GetPolys())
        # Get information from the data.
        pointData = vtk_poly_data.GetPointData()
        numComps = pointData.GetNumberOfComponents()
        numPoints = vtk_poly_data.GetPoints().GetNumberOfPoints()
        numArrays = pointData.GetNumberOfArrays()
        assert field_name is not None or numArrays == 1, "Field name must be specified since data has %d arrays" % numArrays
        # Get the min/max field values.
        array = pointData.GetScalars(field_name)
        fmin, fmax = array.GetRange()
        # Prepare for coloring the points.
        redArray = vtk.vtkUnsignedCharArray()
        greenArray = vtk.vtkUnsignedCharArray()
        blueArray = vtk.vtkUnsignedCharArray()
        redArray.SetNumberOfComponents(numComps)
        greenArray.SetNumberOfComponents(numComps)
        blueArray.SetNumberOfComponents(numComps)
        redArray.SetNumberOfTuples(numPoints)
        greenArray.SetNumberOfTuples(numPoints)
        blueArray.SetNumberOfTuples(numPoints)
        rs = numpy.zeros((numComps,), numpy.uint8)
        gs = numpy.zeros((numComps,), numpy.uint8)
        bs = numpy.zeros((numComps,), numpy.uint8)
        # Build the color map.
        colorMap = ColorMap(fmin, fmax)
        colorMethod = eval('colorMap.' + color_map)
        # Color the points.
        for i in range(numPoints):
            fs = array.GetTuple(i)
            for j in range(numComps):
                rs[j], gs[j], bs[j] = colorMethod(fs[j])
            redArray.SetTuple(i, rs)
            greenArray.SetTuple(i, gs)
            blueArray.SetTuple(i, bs)
        colored_vtk_poly_data = vtk_poly_data_copy.GetPointData()
        for arr in redArray, greenArray, blueArray:
            colored_vtk_poly_data.AddArray(arr)
        colored_vtk_poly_data.GetArray(0).SetName('red')
        colored_vtk_poly_data.GetArray(1).SetName('green')
        colored_vtk_poly_data.GetArray(2).SetName('blue')
        return vtk_poly_data_copy

    def cloneShape(self, shape):
        """
        Clone shape
        @param shape
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
            raise IOError, 'File {} not found'.format(file_name)
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
                vtk_poly_data = self.convertToPolyData( vtk_data )
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

    def rotateShape(self, shape, axis=(1., 0., 0.), angleDeg=0.0):
        """
        Rotate along axis
        @param shape
        @param axis rotation axis
        @param angleDeg angle in degrees
        """
        return shape.rotate(axis, angleDeg)

    def saveShape(self, shape, file_name, file_type):
        """
        Save the shape to a file
        @param file_name file name
        @param file_type either 'ascii' or 'binary'
        @param shape for saving
        """
        vtk_poly_data = self.shapeToVTKPolyData(shape)
        self.saveVtkPolyData(vtk_poly_data, file_name, file_type)

    def saveVtkPolyData(self, vtk_poly_data, file_name, file_type):
        """
        Save the vtk_poly_data to a file
        @param file_name file name
        @param file_type either 'ascii' or 'binary'
        @param vtk_poly_data for saving
        """
        self.writer.SetFileName(file_name)
        if file_type.lower() == 'ascii':
            self.writer.SetFileTypeToASCII()
        else:
            self.writer.SetFileTypeToBinary()
        if vtk.VTK_MAJOR_VERSION >= 6:
            self.writer.SetInputData(vtk_poly_data)
        else:
            self.writer.SetInput(vtk_poly_data)
        self.writer.Write()
        self.writer.Update()

    def shapeFromPolygons(self, polys):
        return CSG.fromPolygons(polys)

    def shapeToPolygons(self, shape):
        return shape.toPolygons()

    def shapeFromVTKPolyData(self, pdata):
        """
        Create a shape from a VTK PolyData object
        @param pdata vtkPolyData instance
        @return shape
        @note field data will get lost
        """
        # store the cell connectivity as CSG polygons
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
            polygons.append(Polygon(verts))
        # instantiate the shape
        return CSG.fromPolygons(polygons)

    def shapeToVTKPolyData(self, shape):
        """
        Convert shape to a VTK polydata object
        """
        verts, polys, count = shape.toVerticesAndPolygons()
        shape.points = vtk.vtkPoints()
        numPoints = len(verts)
        shape.points.SetNumberOfPoints(numPoints)
        for i in range(numPoints):
            shape.points.SetPoint(i, verts[i])
        pdata = vtk.vtkPolyData()
        pdata.SetPoints(shape.points)
        numCells = len(polys)
        pdata.Allocate(numCells, 1)
        ptIds = vtk.vtkIdList()
        for i in range(numCells):
            npts = len(polys[i])
            ptIds.SetNumberOfIds(npts)
            for j in range(npts):
                ptIds.SetId(j, polys[i][j])
            pdata.InsertNextCell(vtk.VTK_POLYGON, ptIds)
        return pdata

    def show(self, shape, windowSizeX=600, windowSizeY=400, filename=''):
        """
        Show the boundary surface or write image to file
        @param windowSizeX number of pixels in x
        @param windowSizeY number of pixels in y
        @param filename write to a file if this keyword
               is present and a non-empty string
        """
        pdata = self.shapeToVTKPolyData(shape)

        # create a rendering window and renderer
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
            print 'WARNING: Cannot call show method -- likely missing VTK'
            return
        renWin.AddRenderer(ren)
        renWin.SetSize(windowSizeX, windowSizeY)
        # create a renderwindowinteractor
        iren.SetRenderWindow(renWin)
        # camera
        xmin, xmax, ymin, ymax, zmin, zmax = pdata.GetBounds()
        lo = numpy.array([xmin, ymin, zmin])
        hi = numpy.array([xmax, ymax, zmax])
        camera.SetFocalPoint(hi)
        center = 0.5*(lo + hi)
        camera.SetPosition(center + hi - lo)
        camera.Zoom(1.0)
        ren.SetActiveCamera(camera)
        # mapper
        if vtk.VTK_MAJOR_VERSION >= 6:
            mapper.SetInputData(pdata)
        else:
            mapper.SetInput(pdata)
        mapper.ScalarVisibilityOff()
        # actor
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(1, 1, 1)
        # add axes
        axesColrs = [(1., 0., 0.,), (0., 1., 0.,), (0., 0., 1.,)]
        for a in axes:
            a.SetShaftRadius(0.01)
            a.SetTipLength(0.2)
            a.SetTipRadius(0.03)
        for at in axesTransf:
            at.PostMultiply()
        # rotate the y and z arrows (initially along x)
        axesTransf[1].RotateZ(90.0)
        axesTransf[2].RotateY(-90.0)
        # scale
        for i in range(3):
            factor = hi[i] - lo[i]
            scale = [1., 1., 1.]
            scale[i] = factor
            axesTransf[i].Scale(scale)
        # translate to loBounds
        for at in axesTransf:
            at.Translate(lo)
        for i in range(3):
            axesTPD[i].SetInputConnection(axes[i].GetOutputPort())
            axesTPD[i].SetTransform(axesTransf[i])
            axesMappers[i].SetInputConnection(axesTPD[i].GetOutputPort())
            axesActors[i].SetMapper(axesMappers[i])
            axesActors[i].GetProperty().SetColor(axesColrs[i])
            ren.AddActor(axesActors[i])
        # assign actor to the renderer
        ren.AddActor(actor)
        # write to file
        # FIX_ME: set this to self.writer.
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
            # fire up interactor
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
        @param vtk_data, data whose vtk_dataset_type one of the VTK_DATASET_TYPES
               besides POLYDATA
        """
        if vtk.VTK_MAJOR_VERSION >= 6:
            self.vtk_geometry_filter.SetInputData(vtk_data)
        else:
            self.vtk_geometry_filter.SetInput(vtk_poly_data)
        self.vtk_geometry_filter.Update()
        vtk_poly_data = self.vtk_geometry_filter.GetOutput()
        return vtk_poly_data

    def chooseVtkGeometryFilter(self, vtk_dataset_type):
        """
        Return a specific vtk geometry filter based on vtk_dataset_type.
        @param vtk_dataset_type, one of the VTK_DATASET_TYPES
        """
        if vtk_dataset_type == 'STRUCTURED_GRID':
            return vtk.vtkStructuredGridGeometryFilter()
        elif vtk_dataset_type == 'POLYDATA':
            return vtk.vtkGeometryFilter()
        elif vtk_dataset_type == 'UNSTRUCTURED_GRID':
            return vtk.vtkUnstructuredGridGeometryFilter()

    def getVtkGeometryFilter(self):
        return self.vtk_geometry_filter

    def setVtkGeometryFilter(self, vtk_dataset_type):
        """
        Set the vtk_geometry_filter.
        @param vtk_dataset_type, one of the VTK_VTK_DATASET_TYPES
        """
        self.vtk_geometry_filter = self.chooseVtkGeometryFilter(vtk_dataset_type)

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

    def setReader(self, file_format=None, vtk_dataset_type=None ):
        """
        Set the reader.
        @param file_format, one of vtk or ply - defauylt to self.file_format
        @param vtk_dataset_type, one of the VTK_VTK_DATASET_TYPES
        """
        if file_format is None:
            file_format = self.file_format
        assert file_format in FILE_FORMATS, "Invalid file_format %s" % str( file_format )
        if file_format == 'vtk':
            if vtk_dataset_type is None:
                vtk_dataset_type = self.vtk_dataset_type
            assert vtk_dataset_type in VTK_DATASET_TYPES, "Invalid vtk_dataset_type %s" % str( vtk_dataset_type )
        self.reader = self.chooseReader(file_format, vtk_dataset_type)

    def chooseWriter(self, file_format, vtk_dataset_type):
        """
        Return a writer based on file_format and possibly vtk_dataset_type.
        @param vtk_dataset_type, None or one of the VTK_DATASET_TYPES
        """
        if file_format == 'ply':
            return vtk.vtkPLYWriter()
        # For now we'll just return the POLYDATA writer since methods work
        # only with that vtk_dataset_type.  The rest is here for possible future
        # use.
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
        assert file_format in FILE_FORMATS, "Invalid file_format %s" % str( file_format )
        if file_format == 'vtk':
            if vtk_dataset_type is None:
                vtk_dataset_type = self.vtk_dataset_type
            assert vtk_dataset_type in VTK_DATASET_TYPES, "Invalid vtk_dataset_type %s" % str( vtk_dataset_type )
        self.writer = self.chooseWriter(file_format, vtk_dataset_type)


###############################################################################
def testPrimitiveShapes():
    shape_mgr = ShapeManager(file_format='vtk', vtk_dataset_type='POLYDATA')

    # Box
    box = shape_mgr.createShape('box', origin=(0.,  0.,  0.), lengths=(0.5,  1.,  2.))
    shape_mgr.saveShape(shape=box, file_name='box.vtk', file_type='ascii')
    shape_mgr.show(box)

    # Cone
    con = shape_mgr.createShape('cone', radius=1.0, origin=(0.,  0.,  0.), lengths=[1., 0., 0.], n_theta=8)
    shape_mgr.saveShape(shape=con, file_name='con.vtk', file_type='ascii')
    shape_mgr.show(con)

    # Cylinder
    cyl = shape_mgr.createShape('cylinder', radius=1.0, origin=(0., 0., 0.), lengths=(1., 0., 0.), n_theta=8)
    shape_mgr.saveShape(shape=cyl, file_name='cyl.vtk', file_type='ascii')
    shape_mgr.show(cyl)

    # Sphere
    sph = shape_mgr.createShape('shpere', radius=1.0, origin=(0., 0., 0.), n_theta=8, n_phi=4)
    shape_mgr.saveShape(shape=sph, file_name='sph.vtk', file_type='ascii')
    shape_mgr.show(sph)


def testSaveLoad():
    shape_mgr = ShapeManager(file_format='vtk', vtk_dataset_type='POLYDATA')
    s = shape_mgr.createShape('shpere', radius=0.7, origin=(0., 0., 0.), n_theta=8, n_phi=4)
    shape_mgr.saveShape(shape=s, file_name='t.vtk', file_type='ascii')
    s2 = shape_mgr.loadAsShape('t.vtk')
    s3 = shape_mgr.createShape('box', origin=(0.1, 0.2, 0.3), lengths=(1.1, 1.2, 1.3))
    s4 = s2 + s3
    s4.debug()


def testConstructiveGeometry():
    shape_mgr = ShapeManager(file_format='ply')
    s1 = shape_mgr.createShape('sphere', radius=0.7, origin=(0., 0., 0.))
    s2 = shape_mgr.createShape('sphere', radius=0.2, origin=(0.1, 0.2, 0.3))
    b = shape_mgr.createShape('box', origin=(0.1, 0.2, 0.3), lengths=(1.1, 1.2, 1.3))
    c = shape_mgr.createShape('cylinder', radius=0.5, origin=(0.3, 0.4, 0.5), lengths=(1.0, 0.0, 0.0))
    geom = c*b - s2 - s1
    shape_mgr.saveShape(shape=geom, file_name='geom.ply', file_type='binary')
    geom2 = shape_mgr.loadAsShape('geom.ply')
    shape_mgr.show(geom2)


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
