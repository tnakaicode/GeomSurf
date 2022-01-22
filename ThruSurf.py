import numpy as np
import matplotlib.pyplot as plt
import sys

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Lin, gp_Sphere, gp_Pln
from OCC.Core.gp import gp_Circ
from OCC.Core.BRep import BRep_Tool, BRep_Builder
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeWire
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_ThruSections, BRepOffsetAPI_MakeOffset, BRepOffsetAPI_MakeEvolved
from OCC.Core.BRepOffset import BRepOffset_MakeOffset, BRepOffset_Skin, BRepOffset_Interval
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Common
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Fuse
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Section
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Cut
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Splitter
from OCC.Core.BRepAlgo import BRepAlgo_Cut, BRepAlgo_Common, BRepAlgo_Tool
from OCC.Core.Geom import Geom_Circle, Geom_Curve
from OCC.Core.GeomAPI import GeomAPI_IntCS
from OCC.Core.GeomLProp import GeomLProp_CurveTool
from OCC.Core.GeomAbs import GeomAbs_C0
from OCC.Core.GeomAbs import GeomAbs_Intersection, GeomAbs_Arc
from OCC.Core.TopoDS import TopoDS_Compound
from OCC.Extend.DataExchange import read_step_file, write_step_file
from OCC.Extend.DataExchange import read_iges_file, write_iges_file
from OCC.Extend.DataExchange import read_stl_file, write_stl_file
from OCCUtils.Construct import make_polygon, make_circle
from OCCUtils.Construct import point_to_vector, vector_to_point
from OCCUtils.Construct import dir_to_vec, vec_to_dir
from OCCUtils.Topology import Topo

from src.base import dispocc, set_loc


class GenThruSurf (plotocc):

    def __init__(self):
        plotocc.__init__(self)
        self.axs = gp_Ax3()

        self.circle = Geom_Circle(gp_Circ(self.axs.Ax2(), 100))
        p0, v1, v2 = gp_Pnt(), gp_Vec(), gp_Vec()
        GeomLProp_CurveTool.D2(self.circle, 0, p0, v1, v2)
        self.poly_axs = gp_Ax3(p0, vec_to_dir(v1))
        for num in range(4, 9):
            self.poly = self.make_PolyWire(
                num=num, radi=20.0, axs=self.poly_axs)

            self.base = self.make_Thru(50)
            self.display.DisplayShape(
                self.base, transparency=0.7, color="BLUE1")
            write_step_file(self.base, self.tmpdir +
                            "ThruSurf_{:d}.stp".format(num))

    def make_Thru(self, num=50):
        api = BRepOffsetAPI_ThruSections()
        print(self.poly.Location().Transformation())
        for idx, phi in enumerate(np.linspace(0, 2 * np.pi, num)):
            ax = self.poly_axs.Rotated(self.axs.Axis(), phi)
            poly_i = self.poly.Located(set_loc(gp_Ax3(), ax))
            # print(poly_i.Location().Transformation())
            api.AddWire(poly_i)
            self.display.DisplayShape(poly_i)
        api.Build()
        return api.Shape()

    def make_PolyWire(self, num=6, radi=1.0, shft=0.0, axs=gp_Ax3()):
        lxy = radi - 0.1
        pnts = []
        angl = 360 / num
        for i in range(num):
            thet = np.deg2rad(i * angl) + np.deg2rad(shft)
            x, y = radi * np.sin(thet), radi * np.cos(thet)
            pnts.append(gp_Pnt(x, y, 0))
        pnts.append(pnts[0])
        poly = make_polygon(pnts)
        poly.Location(set_loc(gp_Ax3(), axs))
        return poly

    def export_file(self):
        builder = BRep_Builder()
        compound = TopoDS_Compound()
        builder.MakeCompound(compound)
        write_step_file(compound, self.tmpdir + "ThruSurf.stp")

    def display_object(self):
        self.display.DisplayShape(self.poly)
        self.display.DisplayShape(self.circle)
        self.show_axs_pln(self.axs, scale=100)
        self.show()


if __name__ == '__main__':
    obj = GenThruSurf()
    obj.display_object()
