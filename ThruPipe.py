import numpy as np
import matplotlib.pyplot as plt
import sys

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Lin, gp_Sphere, gp_Pln
from OCC.Core.gp import gp_Circ, gp_Elips
from OCC.Core.BRep import BRep_Tool, BRep_Builder
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeWire
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_ThruSections, BRepOffsetAPI_MakeOffset, BRepOffsetAPI_MakeEvolved, BRepOffsetAPI_MakePipe, BRepOffsetAPI_MakePipeShell
from OCC.Core.BRepOffset import BRepOffset_MakeOffset, BRepOffset_Skin, BRepOffset_Interval
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Common
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Fuse
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Section
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Cut
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Splitter
from OCC.Core.BRepAlgo import BRepAlgo_Cut, BRepAlgo_Common, BRepAlgo_Tool
from OCC.Core.Geom import Geom_Circle, Geom_Curve
from OCC.Core.GeomAPI import GeomAPI_IntCS
from OCC.Core.GeomAPI import GeomAPI_PointsToBSplineSurface
from OCC.Core.GeomAPI import GeomAPI_PointsToBSpline
from OCC.Core.GeomLProp import GeomLProp_CurveTool
from OCC.Core.GeomAbs import GeomAbs_C0
from OCC.Core.GeomAbs import GeomAbs_Intersection, GeomAbs_Arc
from OCC.Core.TopoDS import TopoDS_Compound
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Core.TColgp import TColgp_Array1OfPnt, TColgp_Array2OfPnt
from OCC.Core.TColgp import TColgp_HArray1OfPnt, TColgp_HArray2OfPnt
from OCC.Extend.DataExchange import read_step_file, write_step_file
from OCC.Extend.DataExchange import read_iges_file, write_iges_file
from OCC.Extend.DataExchange import read_stl_file, write_stl_file
from OCCUtils.Construct import make_polygon, make_circle, make_vertex, make_edge, make_wire
from OCCUtils.Construct import point_to_vector, vector_to_point
from OCCUtils.Construct import dir_to_vec, vec_to_dir
from OCCUtils.Topology import Topo

from src.base import plotocc, set_loc


class GenThruSurf (plotocc):

    def __init__(self):
        plotocc.__init__(self)
        self.axs = gp_Ax3()

        p_array = TColgp_Array1OfPnt(1, 100)
        for idx, pt in enumerate(np.linspace(-2.0, 2.0, 100)):
            pnt = gp_Pnt(300 * np.sin(pt), 100 * np.cos(3 * pt), 0)
            p_array.SetValue(idx + 1, pnt)
        api = GeomAPI_PointsToBSpline(p_array)
        self.curv = api.Curve()
        print(self.curv)

        api = BRepOffsetAPI_ThruSections()
        api.SetSmoothing(True)
        num_list = [
            3, 3, 3, 3, 3,
            6, 6, 6, 6, 6, 6, 6, 6, 6,
            7, 7, 7, 7, 7, 7, 7, 7,
            4, 4, 4
        ]
        p1, v1, v2 = gp_Pnt(), gp_Vec(), gp_Vec()
        for idx, pt in enumerate(np.linspace(0.0, 0.5, len(num_list))):
            GeomLProp_CurveTool.D2(self.curv, pt, p1, v1, v2)
            v3 = v1.Crossed(v2)
            axis = gp_Ax3(p1, vec_to_dir(v1), vec_to_dir(v2))
            poly = self.make_PolyWire(num=num_list[idx], radi=20, axs=axis)
            api.AddWire(poly)
            self.show_axs_pln(axis, scale=10)
            self.display.DisplayShape(poly)
        api.Build()
        self.display.DisplayShape(api.Shape())
        write_step_file(api.Shape(), "./tmp/ThruPipe_Hex.stp")

        api = BRepOffsetAPI_ThruSections()
        # api.SetSmoothing(True)
        p1, v1, v2 = gp_Pnt(), gp_Vec(), gp_Vec()
        for idx, pt in enumerate(np.linspace(0.55, 1.0, 20)):
            GeomLProp_CurveTool.D2(self.curv, pt, p1, v1, v2)
            v3 = v1.Crossed(v2)
            axis = gp_Ax3(p1, vec_to_dir(v1), vec_to_dir(v2))
            shft = 90 * pt
            poly = self.make_Ellip(rxy=[15, 10], shft=shft, axs=axis)
            api.AddWire(poly)
            self.display.DisplayShape(poly)
        api.Build()
        self.display.DisplayShape(api.Shape())
        write_step_file(api.Shape(), "./tmp/ThruPipe_Ellip.stp")

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

    def make_Ellip(self, rxy=[1.0, 1.0], shft=0.0, axs=gp_Ax3()):
        rx, ry = rxy
        if rx > ry:
            major_radi = rx
            minor_radi = ry
            axis = axs.Ax2()
            axis.SetXDirection(axs.XDirection())
        else:
            major_radi = ry
            minor_radi = rx
            axis = axs.Ax2()
            axis.SetXDirection(axs.YDirection())
        axis.Rotate(axs.Axis(), np.deg2rad(shft))
        elip = make_edge(gp_Elips(axis, major_radi, minor_radi))
        poly = make_wire(elip)
        return poly

    def export_file(self):
        builder = BRep_Builder()
        compound = TopoDS_Compound()
        builder.MakeCompound(compound)
        write_step_file(compound, "./tmp/ThruSurf.stp")

    def display_object(self):
        self.display.DisplayShape(self.curv)
        self.show_axs_pln(self.axs, scale=100)
        self.show()


if __name__ == '__main__':
    obj = GenThruSurf()
    obj.display_object()
