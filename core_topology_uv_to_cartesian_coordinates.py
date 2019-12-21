import numpy as np
import matplotlib.pyplot as plt
import sys

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Lin, gp_Sphere
from OCC.Core.BRep import BRep_Tool, BRep_Builder
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse, BRepAlgoAPI_Section
from OCC.Core.Geom import Geom_ToroidalSurface
from OCC.Core.Geom import Geom_Line
from OCC.Core.GeomAPI import GeomAPI_IntCS
from OCC.Core.GeomAPI import GeomAPI_PointsToBSplineSurface
from OCC.Core.GeomAbs import GeomAbs_C2
from OCC.Core.GeomLProp import GeomLProp_SurfaceTool
from OCC.Core.TopoDS import TopoDS_Compound
from OCC.Core.ShapeAnalysis import ShapeAnalysis_Surface, shapeanalysis_GetFaceUVBounds
from OCC.Extend.DataExchange import read_step_file, write_step_file
from OCCUtils.Construct import make_n_sided, make_n_sections
from OCCUtils.Construct import make_edge, make_polygon
from OCCUtils.Construct import vec_to_dir, dir_to_vec
from OCCUtils.Topology import Topo
from OCC.Core.TColgp import TColgp_Array2OfPnt

from base import plotocc, gen_ellipsoid


def build_surf():
    p1 = gp_Pnt(-15, 200, 10)
    p2 = gp_Pnt(5, 204, 0)
    p3 = gp_Pnt(15, 200, 0)
    p4 = gp_Pnt(-15, 20, 15)
    p5 = gp_Pnt(-5, 20, 0)
    p6 = gp_Pnt(15, 20, 35)

    array = TColgp_Array2OfPnt(1, 3, 1, 2)
    array.SetValue(1, 1, p1)
    array.SetValue(2, 1, p2)
    array.SetValue(3, 1, p3)
    array.SetValue(1, 2, p4)
    array.SetValue(2, 2, p5)
    array.SetValue(3, 2, p6)
    bspl_surf = GeomAPI_PointsToBSplineSurface(array, 3, 8, GeomAbs_C2,
                                               0.001).Surface()
    return bspl_surf


def build_points_network(bspl_srf):
    """ Creates a list of gp_Pnt points from a bspline surface
    """
    # first create a face
    face = BRepBuilderAPI_MakeFace(bspl_srf, 1e-6).Face()
    # get face uv bounds
    umin, umax, vmin, vmax = shapeanalysis_GetFaceUVBounds(face)
    print(umin, umax, vmin, vmax)

    pnts = []
    sas = ShapeAnalysis_Surface(bspl_srf)

    u = umin
    while u < umax:
        v = vmin
        while v < vmax:
            p = sas.Value(u, v)
            print("u=", u, " v=", v, "->X=", p.X(), " Y=", p.Y(), " Z=", p.Z())
            pnts.append(p)
            v += 0.1
        u += 0.1
    return pnts


# https://www.opencascade.com/doc/occt-7.4.0/refman/html/class_shape_analysis___surface.html
if __name__ == '__main__':
    surf = build_surf()

    obj = plotocc()
    obj.display.DisplayShape(surf, update=True)
    pts = build_points_network(surf)
    for pt in pts:
        obj.display.DisplayShape(pt, update=False)
    obj.display.Repaint()
    obj.show_axs_pln()
    obj.show()
