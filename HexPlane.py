import numpy as np
import matplotlib.pyplot as plt
import sys

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Lin, gp_Sphere
from OCC.Core.BRep import BRep_Tool, BRep_Builder
from OCC.Core.BRepFill import BRepFill_Filling
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse, BRepAlgoAPI_Section
from OCC.Core.Geom import Geom_ToroidalSurface
from OCC.Core.Geom import Geom_Line
from OCC.Core.GeomAbs import GeomAbs_G1, GeomAbs_G2
from OCC.Core.GeomAbs import GeomAbs_C0, GeomAbs_C1, GeomAbs_C2, GeomAbs_C3
from OCC.Core.GeomAPI import GeomAPI_IntCS
from OCC.Core.GeomLProp import GeomLProp_SurfaceTool
from OCC.Core.TopoDS import TopoDS_Compound
from OCC.Extend.DataExchange import read_step_file, write_step_file
from OCCUtils.Construct import make_n_sided, make_n_sections
from OCCUtils.Construct import make_edge, make_polygon
from OCCUtils.Construct import vec_to_dir, dir_to_vec
from OCCUtils.Topology import Topo

from base import plotocc, gen_ellipsoid, set_loc


class HexPlane (plotocc):

    def __init__(self):
        plotocc.__init__(self)

        axs1 = gp_Ax3(gp_Pnt(5, 0, 0), gp_Dir(0, 0, 1))
        # Fail (not display)
        # self.pln1 = self.make_PolyPlane(radi=2.0, num=5, shft=0.0, axs=axs1)
        # self.pln1 = self.make_PolyPlane(radi=2.0, num=5, shft=0.0, axs=axs1)
        # self.pln1 = self.make_PolyPlane(radi=2.0, num=5, shft=0.0, axs=axs1)
        # self.pln1 = self.make_PolyPlane(radi=2.0, num=5, shft=0.0, axs=axs1)
        # self.pln1 = self.make_PolyPlane(radi=2.0, num=5, shft=0.0, axs=axs1)
        #

        self.pln1 = self.make_PolyPlane(radi=1.0, num=16, shft=0.1, axs=axs1)
        self.pln2 = self.make_PolyPlane(radi=2.0, num=6, shft=30.0)

    def make_PolyPlane(self, num=6, radi=1.0, shft=0.0, axs=gp_Ax3()):
        pnts = []
        angl = 360 / num
        for i in range(num):
            thet = np.deg2rad(i * angl) + np.deg2rad(shft)
            x, y = radi * np.sin(thet), radi * np.cos(thet)
            pnts.append(gp_Pnt(x, y, 0))
        pnts.append(pnts[0])
        poly = make_polygon(pnts)

        n_sided = BRepFill_Filling()
        n_sided.SetResolParam()
        n_sided.SetConstrParam()
        for i, pnt in enumerate(pnts[:-1]):
            n_sided.Add(pnt)
            n_sided.Add(make_edge(pnts[i], pnts[i + 1]), GeomAbs_C0)
        n_sided.Add(gp_Pnt())
        n_sided.Build()
        face = n_sided.Face()
        face.Location(set_loc(gp_Ax3(), axs))
        return face

    def export_file(self):
        builder = BRep_Builder()
        compound = TopoDS_Compound()
        builder.MakeCompound(compound)
        builder.Add(compound, self.pln1)
        builder.Add(compound, self.pln2)
        #builder.Add(compound, make_polygon(self.pts))
        write_step_file(compound, "./tmp/HexPlane.stp")

    def display_ball(self):
        write_step_file(self.pln2, "./tmp/HexPlane.stp")
        self.display.DisplayShape(self.pln1)
        self.display.DisplayShape(self.pln2)
        self.show_axs_pln(scale=1.0)
        self.show()


if __name__ == '__main__':
    obj = HexPlane()
    obj.display_ball()
