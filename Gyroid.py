import numpy as np
from numpy import sin, cos, pi
from skimage import measure

from OCC.Display.SimpleGui import init_display
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakePolygon
from OCC.BRepFill import BRepFill_CurveConstraint
from OCC.GeomAbs import GeomAbs_C0
from OCC.gp import gp_Pnt
from OCC.BRepFill import BRepFill_Filling
from OCC.Extend.TopologyUtils import TopologyExplorer, WireExplorer
from OCC.Extend.ShapeFactory import make_face, make_vertex
from OCC.GeomAbs import GeomAbs_Intersection
from OCC.BRepOffset import BRepOffset_MakeOffset, BRepOffset_Skin
from OCCUtils.Construct import make_n_sided, make_n_sections
from OCCUtils.Construct import make_edge


# DEFINE GEOMETRY
def gyroid(x, y, z, t):
    return cos(x) * sin(y) + cos(y) * sin(z) + cos(z) * sin(x) - t  # Gyroïd


lattice_param = 1.
size = 1.
strut_param = 0.
resolution = 6j  # >6 crashes
res = int(np.imag(resolution))
x, y, z = pi / 2 * np.mgrid[-1:1:resolution, -
                            1:1:resolution, -1:1:resolution] * lattice_param
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

# OCC FUNCTIONS
# def make_n_sided(edges, points, continuity=GeomAbs_C0):
#    n_sided = BRepFill_Filling(2, len(points), 2, False, 0.00001, 0.0001, 0.01, 0.1, 8, 9)
#    for edg in edges:
#        n_sided.Add(edg, continuity)
#    for pt in points:
#        n_sided.Add(pt)
#    n_sided.Build()
#    face = n_sided.Face()
#    return face


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
    edges = [i for i in TopologyExplorer(poly).edges()]

    # POINTS CONSTRAINT
    p_array = []
    for pt in gxyz:
        x, y, z = pt
        if (x < xmax) and (x > xmin) and (y < ymax) and (y > ymin) and (z < zmax) and (z > zmin):
            p = gp_Pnt(float(x), float(y), float(z))
            p_array.append(p)

    face = make_n_sided(edges, p_array)

    # THICKEN SURFACE
    thickness = 0.15
    solid = BRepOffset_MakeOffset()
    solid.Initialize(face, thickness, 1.e-5, BRepOffset_Skin, False, False,
                     GeomAbs_Intersection, True)  # The last True is important to make solid
    solid.MakeOffsetShape()
    aShape = solid.Shape()

    display.DisplayShape(edges)
    for p in p_array:
        display.DisplayShape(make_vertex(p))
    display.DisplayShape(face, update=True)
    display.DisplayShape(aShape, update=True)


# DISPLAY ALL
display, start_display, add_menu, add_function_to_menu = init_display()
geom_plate()
display.FitAll()
start_display()
