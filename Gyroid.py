import numpy as np
from numpy import sin, cos, pi
from skimage import measure

from OCC.Display.SimpleGui import init_display
from OCC.CoreBRepBuilderAPI import BRepBuilderAPI_MakePolygon
from OCC.CoreBRepFill import BRepFill_CurveConstraint
from OCC.CoreGeomAbs import GeomAbs_C0
from OCC.Coregp import gp_Pnt
from OCC.CoreBRepFill import BRepFill_Filling
from OCC.Extend.TopologyUtils import TopologyExplorer, WireExplorer
from OCC.Extend.ShapeFactory import make_face, make_vertex
from OCC.CoreGeomAbs import GeomAbs_Intersection
from OCC.CoreBRepOffset import BRepOffset_MakeOffset, BRepOffset_Skin
from OCCUtils.Construct import make_n_sided, make_n_sections
from OCCUtils.Construct import make_edge
from OCCUtils.Topology import Topo


# DEFINE GEOMETRY
def gyroid(x, y, z, t):
    # Gyroïd
    return cos(x) * sin(y) + cos(y) * sin(z) + cos(z) * sin(x) - t


size = 1.
strut_param = 1.
resolution = 6j  # >6 crashes
res = int(np.imag(resolution))
x, y, z = pi / 2 * np.mgrid[-1:1:resolution, -1:1:resolution, -1:1:resolution]
vol = gyroid(x, y, z, strut_param) * size
verts, faces, norm, val = measure.marching_cubes_lewiner(vol, level=None)
xv, yv, zv = verts[:, 0], verts[:, 1], verts[:, 2]
gxyz = np.array([[x - max(xv) / 2, y - max(yv) / 2, z - max(zv) / 2]
                 for x, y, z in zip(xv, yv, zv)])
xmax, ymax, zmax = max(gxyz[:, 0]), max(gxyz[:, 1]), max(gxyz[:, 2])
xmin, ymin, zmin = min(gxyz[:, 0]), min(gxyz[:, 1]), min(gxyz[:, 2])
# FIND EDGES
gx_min = gxyz[np.nonzero((gxyz[:, 0] == xmin) & (gxyz[:, 2] > 0))]
gy_min = gxyz[np.nonzero((gxyz[:, 1] == ymin) & (gxyz[:, 0] > 0))]
gz_min = gxyz[np.nonzero((gxyz[:, 2] == zmin) & (gxyz[:, 1] > 0))]
gx_max = gxyz[np.nonzero((gxyz[:, 0] == xmax) & (gxyz[:, 2] < 0))]
gy_max = gxyz[np.nonzero((gxyz[:, 1] == ymax) & (gxyz[:, 0] < 0))]
gz_max = gxyz[np.nonzero((gxyz[:, 2] == zmax) & (gxyz[:, 1] < 0))]
# ORDER EDGES POINTS
gx_min = list(gx_min[np.argsort(np.linalg.norm(gx_min - gx_min[0], axis=1))])
gy_min = list(gy_min[np.argsort(np.linalg.norm(gy_min - gy_min[-1], axis=1))])
gz_min = list(gz_min[np.argsort(np.linalg.norm(gz_min - gz_min[0], axis=1))])
gx_max = list(gx_max[np.argsort(np.linalg.norm(gx_max - gx_max[0], axis=1))])
gy_max = list(gy_max[np.argsort(np.linalg.norm(gy_max - gy_max[0], axis=1))])
gz_max = list(gz_max[np.argsort(np.linalg.norm(gz_max - gz_max[0], axis=1))])
# BUILD EXTERNAL EDGE
xyz_max = gx_max + gz_min[::-1] + gy_max + gx_min[::-1] + gz_max + gy_min


def make_closed_polygon(*args):
    poly = BRepBuilderAPI_MakePolygon()
    for pt in args:
        if isinstance(pt, list) or isinstance(pt, tuple):
            for i in pt:
                poly.Add(i)
        else:
            poly.Add(pt)
    poly.Build()
    poly.Close()
    result = poly.Wire()
    return result


def geom_plate(event=None):
    display.EraseAll()

    # EDGE CONSTRAINT
    e_array = []
    for e in xyz_max:
        x, y, z = e
        e = gp_Pnt(float(x), float(y), float(z))
        e_array.append(e)
    poly = make_closed_polygon(e_array)
    top = Topo(poly)

    edges = []
    for e in top.edges():
        edges.append(e)

    # POINTS CONSTRAINT
    p_array = []
    for pt in gxyz:
        x, y, z = pt
        if (x < xmax) and (x > xmin) and (y < ymax) and (y > ymin) and (z < zmax) and (z > zmin):
            p = gp_Pnt(float(x), float(y), float(z))
            p_array.append(p)
    
    print(len(edges), len(p_array))
    face = make_n_sided(edges, p_array)

    # THICKEN SURFACE
    thickness = 0.15
    solid = BRepOffset_MakeOffset()
    solid.Initialize(face, thickness, 1.0E-5, BRepOffset_Skin, False, False,
                     GeomAbs_Intersection, True)  # The last True is important to make solid
    # solid.MakeOffsetShape()
    aShape = solid.Shape()

    display.DisplayShape(poly)
    for p in p_array:
        display.DisplayShape(make_vertex(p))
    display.DisplayShape(face)
    #display.DisplayShape(aShape, update=True)


# DISPLAY ALL
display, start_display, add_menu, add_function_to_menu = init_display()
geom_plate()
display.FitAll()
start_display()
