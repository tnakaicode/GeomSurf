import numpy as np
import matplotlib.pyplot as plt
import sys

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Lin, gp_Sphere
from OCC.Core.BRep import BRep_Tool, BRep_Builder
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse, BRepAlgoAPI_Section
from OCC.Core.Geom import Geom_ToroidalSurface
from OCC.Core.Geom import Geom_Line
from OCC.Core.GeomAPI import GeomAPI_IntCS
from OCC.Core.GeomLProp import GeomLProp_SurfaceTool
from OCC.Core.TopoDS import TopoDS_Compound
from OCC.Extend.DataExchange import read_step_file, write_step_file
from OCCUtils.Construct import make_n_sided, make_n_sections
from OCCUtils.Construct import make_edge, make_polygon
from OCCUtils.Construct import vec_to_dir, dir_to_vec
from OCCUtils.Topology import Topo

from base import plotocc, gen_ellipsoid


class HexPlane (plotocc):

    def __init__(self):
        plotocc.__init__(self)

        self.pts = []
        num = 6
        radi = 1.0
        angl = 360 / num
        shft = 0.0
        for i in range(num):
            thet = np.deg2rad(i * angl) + shft
            x, y = radi * np.sin(thet), radi * np.cos(thet)
            self.pts.append(gp_Pnt(x, y, 0))
        self.pts.append(self.pts[0])

    def export_file(self):
        builder = BRep_Builder()
        compound = TopoDS_Compound()
        builder.MakeCompound(compound)
        builder.Add(compound, self.b1)
        builder.Add(compound, self.b2)
        builder.Add(compound, make_polygon(self.pts))
        write_step_file(compound, "./tmp/test.stp")

    def display_ball(self):
        self.display.DisplayShape(make_polygon(self.pts))
        self.show_axs_pln(scale=1.0)
        self.show()
        # self.export_file()


if __name__ == '__main__':
    obj = HexPlane()
    obj.display_ball()
