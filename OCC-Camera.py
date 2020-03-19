import numpy as np
import matplotlib.pyplot as plt
import sys

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Lin, gp_Sphere, gp_Pln
from OCC.Core.BRep import BRep_Tool, BRep_Builder, BRep_PointsOnSurface
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse, BRepAlgoAPI_Section
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCone
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeWire
from OCC.Core.Geom import Geom_Line
from OCC.Core.GeomAbs import GeomAbs_G1, GeomAbs_G2
from OCC.Core.GeomAbs import GeomAbs_C0, GeomAbs_C1, GeomAbs_C2, GeomAbs_C3
from OCC.Core.GeomAPI import GeomAPI_IntCS, GeomAPI_IntSS
from OCC.Core.GeomLProp import GeomLProp_SurfaceTool
from OCC.Core.TopoDS import TopoDS_Compound, TopoDS_Builder, TopoDS_Face, TopoDS_Shell
from OCC.Core.math import math_NewtonMinimum
from OCC.Extend.DataExchange import read_step_file, write_step_file
from OCCUtils.Construct import make_n_sided, make_n_sections
from OCCUtils.Construct import make_edge, make_polygon, make_vertex, make_box
from OCCUtils.Construct import vec_to_dir, dir_to_vec
from OCCUtils.Construct import point_to_vector, vector_to_point
from OCCUtils.Topology import Topo, dumpTopology

from base import plotocc, gen_ellipsoid, set_loc, spl_face


def line_from_axs(axs=gp_Ax3(), length=100):
    vec = point_to_vector(axs.Location()) + \
        dir_to_vec(axs.Direction()) * length
    return make_edge(axs.Location(), vector_to_point(vec))


class Cap (plotocc):

    def __init__(self):
        plotocc.__init__(self)
        self.compound = TopoDS_Compound()
        self.builder = BRep_Builder()
        self.builder.MakeCompound(self.compound)

        box = BRepPrimAPI_MakeBox(100, 200, 300).Shape()
        self.display.DisplayShape(box)

        con = BRepPrimAPI_MakeCone(
            gp_Ax2(gp_Pnt(500, 0, 0), gp_Dir(0, 0, 1)), 100, 200, 300).Shape()
        self.display.DisplayShape(con)

        self.camera = self.display.View.Camera()
        print(self.camera.Scale(), dir_to_vec(self.camera.OrthogonalizedUp()))
        print(self.camera.ViewDimensions())

    def export_capture(self):
        self.display.FitAll()
        self.display.View.Dump("./tmp/cap.png")

        self.camera.OrthogonalizeUp()
        self.camera.SetUp(gp_Dir(0, 0, 1))
        print(self.camera.Scale(), self.camera.Eye(), self.camera.Center())

        self.camera.SetDirection(gp_Dir(-1, 0, 0))
        print(self.camera.Scale(), self.camera.Eye(), self.camera.Center())
        self.display.View.SetCamera(self.camera)
        self.display.View.Dump("./tmp/cap_x.png")

        self.camera.SetDirection(gp_Dir(0, -1, 0))
        print(self.camera.Scale(), self.camera.Eye(), self.camera.Center())
        self.display.View.SetCamera(self.camera)
        self.display.View.Dump("./tmp/cap_y.png")

        self.camera.SetDirection(gp_Dir(0, 0, -1))
        print(self.camera.Scale(), self.camera.Eye(), self.camera.Center())
        self.display.View.SetCamera(self.camera)
        self.display.View.Dump("./tmp/cap_z.png")

    def display_Plane(self):
        self.show_axs_pln(scale=200.0)
        self.export_capture()
        self.start_display()


if __name__ == '__main__':
    obj = Cap()
    obj.display_Plane()
