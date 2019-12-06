import numpy as np
import matplotlib.pyplot as plt
import sys

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Lin, gp_Sphere, gp_Pln
from OCC.Core.BRep import BRep_Tool, BRep_Builder, BRep_PointsOnSurface
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse, BRepAlgoAPI_Section
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeWire
from OCC.Core.Geom import Geom_Line
from OCC.Core.GeomAbs import GeomAbs_G1, GeomAbs_G2
from OCC.Core.GeomAbs import GeomAbs_C0, GeomAbs_C1, GeomAbs_C2, GeomAbs_C3
from OCC.Core.GeomAPI import GeomAPI_IntCS, GeomAPI_IntSS
from OCC.Core.GeomLProp import GeomLProp_SurfaceTool
from OCC.Core.TopoDS import TopoDS_Compound, TopoDS_Builder, TopoDS_Face, TopoDS_Shell
from OCC.Core.math import math_Vector, math_Matrix, math_Jacobi
from OCC.Extend.DataExchange import read_step_file, write_step_file
from OCCUtils.Construct import make_n_sided, make_n_sections
from OCCUtils.Construct import make_edge, make_polygon, make_vertex
from OCCUtils.Construct import vec_to_dir, dir_to_vec
from OCCUtils.Construct import point_to_vector, vector_to_point
from OCCUtils.Topology import Topo, dumpTopology

from base import plotocc, gen_ellipsoid, set_loc, spl_face


def line_from_axs(axs=gp_Ax3(), length=100):
    vec = point_to_vector(axs.Location()) + \
        dir_to_vec(axs.Direction()) * length
    return make_edge(axs.Location(), vector_to_point(vec))


class ODE1 (plotocc):

    def __init__(self):
        plotocc.__init__(self)
        self.compound = TopoDS_Compound()
        self.builder = BRep_Builder()
        self.builder.MakeCompound(self.compound)

    def export_file(self):
        write_step_file(self.compound, "./tmp/ODE1.stp")

    def display_Plane(self):
        self.show_axs_pln(scale=1.0)
        self.show()


if __name__ == '__main__':
    obj = ODE1()
    obj.export_file()
    obj.display_Plane()
