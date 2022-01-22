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
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE
from OCC.Core.ShapeAnalysis import ShapeAnalysis_Surface, shapeanalysis_GetFaceUVBounds
from OCC.Extend.DataExchange import read_step_file, write_step_file
from OCCUtils.Construct import make_n_sided, make_n_sections
from OCCUtils.Construct import make_edge, make_polygon
from OCCUtils.Construct import vec_to_dir, dir_to_vec
from OCCUtils.Topology import Topo
from OCC.Core.TColgp import TColgp_Array2OfPnt

from src.base_occ import dispocc, gen_ellipsoid


class SurfBox (dispocc):

    def __init__(self):
        dispocc.__init__(self)
        self.shell = read_step_file(self.tmpdir + "SurfUV.stp")
        print(self.shell)
        top = TopExp_Explorer(self.shell, TopAbs_FACE)
        self.face = top.Current()
        print(top.Depth())
        print(self.face)
        self.surf = BRep_Tool.Surface(self.face)

        u0, u1, v0, v1 = shapeanalysis_GetFaceUVBounds(self.face)
        print(u0, u1, v0, v1)
        sas = ShapeAnalysis_Surface(self.surf)
        print(sas.Value(u0, v0))
        print(sas.Value(u0, v1))
        print(sas.Value(u1, v0))
        print(sas.Value(u1, v1))

        u = u0
        while u <= u1:
            v = v0
            while v <= v1:
                p = sas.Value(u, v)
                self.display.DisplayShape(p, update=False)
                v += 1 / 3
            u += 1 / 4

    def show_object(self):
        self.display.DisplayShape(self.face, update=True)
        self.display.Repaint()
        self.show_axs_pln()
        self.show_occ()


if __name__ == '__main__':
    #
    # https://www.opencascade.com/doc/occt-7.4.0/refman/html/class_shape_analysis___surface.html
    # https://www.opencascade.com/doc/occt-7.4.0/refman/html/class_b_rep_proj___projection.html
    #
    obj = SurfBox()
    obj.show_object()
