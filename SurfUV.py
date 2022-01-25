import numpy as np
import matplotlib.pyplot as plt
import sys

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Lin, gp_Sphere
from OCC.Core.BRep import BRep_Tool, BRep_Builder
from OCC.Core.BRepProj import BRepProj_Projection
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse, BRepAlgoAPI_Section
from OCC.Core.Geom import Geom_ToroidalSurface
from OCC.Core.Geom import Geom_Line
from OCC.Core.GeomAPI import GeomAPI_IntCS
from OCC.Core.GeomAPI import GeomAPI_PointsToBSplineSurface
from OCC.Core.GeomAbs import GeomAbs_C2
from OCC.Core.GeomLProp import GeomLProp_SurfaceTool
from OCC.Core.ProjLib import ProjLib_ProjectOnSurface
from OCC.Core.GeomProjLib import geomprojlib
from OCC.Core.TopoDS import TopoDS_Compound
from OCC.Core.ShapeAnalysis import ShapeAnalysis_Surface, shapeanalysis_GetFaceUVBounds
from OCC.Extend.DataExchange import read_step_file, write_step_file
from OCCUtils.Construct import make_n_sided, make_n_sections
from OCCUtils.Construct import make_edge, make_polygon
from OCCUtils.Construct import vec_to_dir, dir_to_vec
from OCCUtils.Topology import Topo
from OCC.Core.TColgp import TColgp_Array2OfPnt

from src.base_occ import dispocc, gen_ellipsoid


class SurfUV (dispocc):

    def __init__(self):
        dispocc.__init__(self)
        self.build_surf()
        write_step_file(self.bspl_face, self.tmpdir + "SurfUV.stp")

    def build_surf(self):
        p1 = gp_Pnt(-15, 200, 10)
        p2 = gp_Pnt(5, 204, 0)
        p3 = gp_Pnt(15, 200, 0)
        p4 = gp_Pnt(-15, 20, 15)
        p5 = gp_Pnt(-5, 20, 0)
        p6 = gp_Pnt(15, 20, 35)
        self.display.DisplayShape(p1, color="RED")
        self.display.DisplayShape(p2, color="RED")
        self.display.DisplayShape(p3, color="RED")
        self.display.DisplayShape(p4, color="RED")
        self.display.DisplayShape(p5, color="RED")
        self.display.DisplayShape(p6, color="RED")

        array = TColgp_Array2OfPnt(1, 3, 1, 2)
        array.SetValue(1, 1, p1)
        array.SetValue(2, 1, p2)
        array.SetValue(3, 1, p3)
        array.SetValue(1, 2, p4)
        array.SetValue(2, 2, p5)
        array.SetValue(3, 2, p6)
        self.bspl_surf = GeomAPI_PointsToBSplineSurface(array, 3, 8, GeomAbs_C2,
                                                        0.001).Surface()
        self.bspl_face = BRepBuilderAPI_MakeFace(self.bspl_surf, 1e-6).Face()

    def build_points_network(self):
        """ Creates a list of gp_Pnt points from a bspline surface
        """
        # get face uv bounds
        umin, umax, vmin, vmax = shapeanalysis_GetFaceUVBounds(self.bspl_face)
        print(umin, umax, vmin, vmax)

        pnts = []
        sas = ShapeAnalysis_Surface(self.bspl_surf)

        u = umin
        while u < umax:
            v = vmin
            while v < vmax:
                p = sas.Value(u, v)
                print("u=", u, " v=", v, "->X=", p.X(),
                      " Y=", p.Y(), " Z=", p.Z())
                self.display.DisplayShape(p, update=False)
                pnts.append(p)
                v += 0.1
            u += 0.1

    def show_object(self):
        self.display.DisplayShape(self.bspl_surf, update=True)
        self.display.Repaint()
        self.show_axs_pln()
        self.ShowOCC()


# https://www.opencascade.com/doc/occt-7.4.0/refman/html/class_shape_analysis___surface.html
# https://www.opencascade.com/doc/occt-7.4.0/refman/html/class_b_rep_proj___projection.html
if __name__ == '__main__':
    obj = SurfUV()
    obj.build_points_network()
    obj.show_object()
