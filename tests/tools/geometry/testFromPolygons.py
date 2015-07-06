from icqsol.tools.geometry.icqShape import Shape

from csg.core import CSG

cube = CSG.cube()
shp = Shape()
shp.fromPolygons(cube.toPolygons())
