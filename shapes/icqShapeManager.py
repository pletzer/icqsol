#!/usr/bin/env python
"""
@brief A base class for constructing shapes
"""
import os
import vtk
import numpy
from csg.geom import Vector, Vertex, Polygon, BSPNode
from csg.core import CSG
from icqShape import Box, Cone, Cylinder, Sphere
from icqShape import DEFAULTS, CompositeShape, Shape


class ShapeManager(object):

    def createShape(self, type, origin=None, lengths=None, radius=None,
                    angle=None, n_theta=None, n_phi=None):
        origin = origin or DEFAULTS.get('origin', [0.0, 0.0, 0.0])
        lengths = lengths or DEFAULTS.get('lengths', [1.0, 1.0, 1.0])
        radius = radius or DEFAULTS.get('radius', [0.5, 0.5, 0.5])
        angle = angle or DEFAULTS.get('angle', 90.0)
        n_theta = n_theta or DEFAULTS.get('n_theta', 16)
        n_phi = n_phi or DEFAULTS.get('n_phi', 8)
        if type == 'box':
            return Box(origin, lengths)
        if type == 'cone':
            return Cone(radius, origin, lengths, n_theta)
        if type == 'cylinder':
            return Cylinder(radius, origin, lengths, n_theta)
        if type == 'sphere':
            return Sphere(radius, origin, n_theta, n_phi)
        return None

    def cloneShape(self, shape):
        """
        Clone shape
        @param shape
        @return new shape
        """
        return Shape(shape.csg.clone())

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
        a = BSPNode(shape.csg.clone().polygons)
        b = BSPNode(other.csg.clone().polygons)
        b.invert()
        a.clipTo(b)
        return Shape(CSG.fromPolygons(a.allPolygons()))

    def load(self, file_name):
        """
        Load geometry from file
        @param file_name file name, suffix should be .ply or .vtk
        @return a Shape object
        """
        if not os.path.exists(file_name):
            raise IOError, 'File {} not found'.format(file_name)
        reader = None
        # select the reader according to the suffix
        if file_name.lower().find('.ply') >= 0:
            reader = vtk.vtkPLYReader()
        else:
            reader = vtk.vtkPolyDataReader()
        reader.SetFileName(file_name)
        # read
        reader.Update()
        # vtkPolyData
        pdata = reader.GetOutput()
        return self.shapeFromVTKPolyData(pdata)

    def rotateShape(self, shape, axis=(1., 0., 0.), angleDeg=0.0):
        """
        Rotate along axis
        @param shape
        @param axis rotation axis
        @param angleDeg angle in degrees
        """
        return shape.csg.rotate(axis, angleDeg)

    def save(self, shape, file_name, file_format, file_type):
        """
        Save the shape in file
        @param file_name file name
        @param file_format file format, currently either VTK or PLY
        @param file_type either 'ascii' or 'binary'
        """
        writer = None
        if file_format.lower() == 'ply':
            writer = vtk.vtkPLYWriter()
        else:
            writer = vtk.vtkPolyDataWriter()
        writer.SetFileName(file_name)
        if file_type.lower() == 'ascii':
            writer.SetFileTypeToASCII()
        else:
            writer.SetFileTypeToBinary()
        pdata = self.shapeToVTKPolyData(shape)
        if vtk.VTK_MAJOR_VERSION >= 6:
            writer.SetInputData(pdata)
        else:
            writer.SetInput(pdata)
        writer.Write()
        writer.Update()

    def shapeFromPolygons(self, polys):
        return Shape(csg=CSG.fromPolygons(polys))

    def shapeToPolygons(self, shape):
        return shape.csg.toPolygons()

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
        csg = CSG.fromPolygons(polygons)
        return Shape(csg=csg)

    def shapeToVTKPolyData(self, shape):
        """
        Convert shape to a VTK polydata object
        """
        verts, polys, count = shape.csg.toVerticesAndPolygons()
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
        shape.csg.translate(disp)


###############################################################################
def testPrimitiveShapes():
    shape_mgr = ShapeManager()

    # Box
    box = shape_mgr.createShape('box', origin=(0.,  0.,  0.),
                                lengths=(0.5,  1.,  2.))
    shape_mgr.save(box, 'box.vtk', file_format='vtk', file_type='ascii')
    shape_mgr.show(box)

    # Cone
    con = shape_mgr.createShape('cone', radius=1.0, origin=(0.,  0.,  0.),
                                lengths=[1., 0., 0.], n_theta=8)
    shape_mgr.save(con, 'con.vtk', file_format='vtk', file_type='ascii')
    shape_mgr.show(con)

    # Cylinder
    cyl = shape_mgr.createShape('cylinder', radius=1.0, origin=(0., 0., 0.),
                                lengths=(1., 0., 0.), n_theta=8)
    shape_mgr.save(cyl, 'cyl.vtk', file_format='vtk', file_type='ascii')
    shape_mgr.show(cyl)

    # Sphere
    sph = shape_mgr.createShape('shpere', radius=1.0, origin=(0., 0., 0.),
                                n_theta=8, n_phi=4)
    shape_mgr.save(sph, 'sph.vtk', file_format='vtk', file_type='ascii')
    shape_mgr.show(sph)


def testSaveLoad():
    shape_mgr = ShapeManager()
    s = shape_mgr.createShape('shpere', radius=0.7, origin=(0., 0., 0.),
                              n_theta=8, n_phi=4)
    shape_mgr.save(s, file_name='t.vtk', file_format='vtk', file_type='ascii')
    s2 = shape_mgr.load('t.vtk')
    s3 = shape_mgr.createShape('box', origin=(0.1, 0.2, 0.3),
                               lengths=(1.1, 1.2, 1.3))
    s4 = s2 + s3
    s4.debug()


def testConstructiveGeometry():
    shape_mgr = ShapeManager()
    s1 = shape_mgr.createShape('sphere', radius=0.7, origin=(0., 0., 0.))
    s2 = shape_mgr.createShape('sphere', radius=0.2, origin=(0.1, 0.2, 0.3))
    b = shape_mgr.createShape('box', origin=(0.1, 0.2, 0.3),
                              lengths=(1.1, 1.2, 1.3))
    c = shape_mgr.createShape('cylinder', radius=0.5, origin=(0.3, 0.4, 0.5),
                              lengths=(1.0, 0.0, 0.0))
    geom = c*b - s2 - s1
    shape_mgr.save(geom, 'geom.ply', file_format='ply', file_type='binary')
    geom2 = shape_mgr.load('geom.ply')
    shape_mgr.show(geom2)


def testShapeComposition():
    shape_mgr = ShapeManager()
    s1 = shape_mgr.createShape('sphere', radius=1, origin=(0., 0., 0.))
    s2 = shape_mgr.createShape('sphere', radius=1.2, origin=(0.8, 0., 0.))
    s3 = shape_mgr.composeShapes([('s1', s1), ('s2', s2)], 's1 + s2')
    shape_mgr.save(s3, 's3.vtk', file_format='vtk', file_type='ascii')


if __name__ == '__main__':
    testPrimitiveShapes()
    testSaveLoad()
    testConstructiveGeometry()
    testShapeComposition()
