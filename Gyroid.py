import numpy as np
from numpy import sin, cos, pi
from skimage import measure

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pnt
from OCC.Core.GeomAbs import GeomAbs_C0
from OCC.Core.GeomAbs import GeomAbs_Intersection, GeomAbs_Arc
from OCC.Core.BRepFill import BRepFill_Filling
from OCC.Core.BRepFill import BRepFill_CurveConstraint
from OCC.Core.BRepOffset import BRepOffset_MakeOffset, BRepOffset_Skin, BRepOffset_Interval
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakePolygon
from OCC.Extend.TopologyUtils import TopologyExplorer, WireExplorer
from OCC.Extend.ShapeFactory import make_face, make_vertex
from OCC.Extend.DataExchange import read_step_file, write_step_file
from OCCUtils.Construct import make_n_sided, make_n_sections
from OCCUtils.Construct import make_edge, make_polygon
from OCCUtils.Topology import Topo

from src.base_occ import dispocc

# https://www.opencascade.com/doc/occt-7.4.0/refman/html/class_b_rep_fill___filling.html
# N-Side Filling This algorithm avoids to build a face from:
#   a set of edges defining the bounds of the face and some constraints the surface support has to satisfy
#   a set of edges and points defining some constraints the support surface has to satisfy
#       an initial surface to deform for satisfying the constraints
#
#   a set of parameters to control the constraints.
#   The support surface of the face is computed by deformation of the initial surface in order to satisfy the given constraints.
#   The set of bounding edges defines the wire of the face.
#
#   If no initial surface is given, the algorithm computes it automatically.
#   If the set of edges is not connected (Free constraint) missing edges are automatically computed.
#
#   Limitations:
#
#   If some constraints are not compatible
#       The algorithm does not take them into account.
#       So the constraints will not be satisfyed in an area containing the incompatibilitries.
#       he constraints defining the bound of the face have to be entered in order to have a continuous wire.
#
#   Other Applications:
#   Deformation of a face to satisfy internal constraints
#   Deformation of a face to improve Gi continuity with connected faces
#
# BRepFill_Filling Example
#
# n_sided = BRepFill_Filling()
# n_sided.SetApproxParam()
# n_sided.SetResolParam()
# n_sided.SetConstrParam()
# for i, pnt in enumerate(pnts[:-1]):
#    i0, i1 = i, i + 1
#    n_sided.Add(pnt)
#    n_sided.Add(make_edge(pnts[i0], pnts[i1]), GeomAbs_C0)
# n_sided.Add(gp_Pnt(0, 0, 1))
# n_sided.Build()
# face = n_sided.Face()


def gyroid(x, y, z, t):
    # Gyroïd
    # DEFINE GEOMETRY
    return cos(x) * sin(y) + cos(y) * sin(z) + cos(z) * sin(x) - t


size = 1.
strut_param = 1.
resolution = 10j  # >6 crashes
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


class Gyroid (dispocc):

    def __init__(self):
        dispocc.__init__(self)

        print(gxyz.shape)
        for i, xyz in enumerate(gxyz):
            print(i, *xyz)
            self.display.DisplayShape(gp_Pnt(*xyz))

        e_array = []
        for e in xyz_max:
            x, y, z = e
            e = gp_Pnt(float(x), float(y), float(z))
            e_array.append(e)
        e_array.append(e_array[0])
        poly = make_polygon(e_array)

        n_sided = BRepFill_Filling()
        for e in Topo(poly).edges():
            n_sided.Add(e, GeomAbs_C0)
        for pt in gxyz:
            x, y, z = pt
            if (x < xmax) and (x > xmin) and (y < ymax) and (y > ymin) and (z < zmax) and (z > zmin):
                n_sided.Add(gp_Pnt(x, y, z))
        n_sided.Build()
        face = n_sided.Face()

        #face = make_n_sided(edges, p_array)

        # THICKEN SURFACE
        thickness = 0.15
        solid = BRepOffset_MakeOffset(
            face, thickness, 1.0E-5, BRepOffset_Skin, False, False, GeomAbs_Intersection, True)
        # The last True is important to make solid
        # solid.MakeOffsetShape()
        # solid.MakeThickSolid()
        #aShape = solid.Shape()

        self.display.DisplayShape(poly)
        self.display.DisplayShape(face)
        #display.DisplayShape(aShape, update=True)
        #write_step_file(aShape, "./tmp/gyroid.stp")

        self.export_stp(solid.Shape())


if __name__ == '__main__':
    obj = Gyroid()
    obj.show_axs_pln(scale=1.0)
    obj.show_occ()
