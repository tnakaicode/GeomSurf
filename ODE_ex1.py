import numpy as np
import matplotlib.pyplot as plt
import sys

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_XYZ
from OCC.Core.gp import gp_Lin, gp_Sphere, gp_Pln
from OCC.Core.BRep import BRep_Tool, BRep_Builder, BRep_PointsOnSurface
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse, BRepAlgoAPI_Section
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeWire
from OCC.Core.Geom import Geom_Line
from OCC.Core.GeomAbs import GeomAbs_G1, GeomAbs_G2
from OCC.Core.GeomAbs import GeomAbs_C0, GeomAbs_C1, GeomAbs_C2, GeomAbs_C3
from OCC.Core.GeomAPI import GeomAPI_IntCS, GeomAPI_IntSS
from OCC.Core.GeomLProp import GeomLProp_SurfaceTool
from OCC.Core.GeomFill import GeomFill_Filling
from OCC.Core.TopoDS import TopoDS_Compound, TopoDS_Builder, TopoDS_Face, TopoDS_Shell
from OCC.Core.math import math_Vector, math_Matrix, math_Jacobi
from OCC.Core.math import math_NewtonMinimum, math_BracketMinimum
from OCC.Core.math import math_GaussMultipleIntegration
from OCC.Core.math import math_MultipleVarFunctionWithHessian
from OCC.Core.math import math_MultipleVarFunctionWithGradient
from OCC.Core.math import math_TrigonometricEquationFunction
from OCC.Core.FEmTool import FEmTool_Curve
from OCC.Core.FairCurve import FairCurve_Newton
from OCC.Extend.DataExchange import read_step_file, write_step_file
from OCCUtils.Construct import make_n_sided, make_n_sections
from OCCUtils.Construct import make_edge, make_polygon, make_vertex
from OCCUtils.Construct import vec_to_dir, dir_to_vec
from OCCUtils.Construct import point_to_vector, vector_to_point
from OCCUtils.Topology import Topo, dumpTopology

from base import plotocc


class ODE1 (plotocc):

    def __init__(self):
        plotocc.__init__(self)
        self.compound = TopoDS_Compound()
        self.builder = BRep_Builder()
        self.builder.MakeCompound(self.compound)

        mmat = math_Matrix(1, 4, 1, 3, 0.0)
        mmat.SetDiag(1.0)
        print(mmat.Determinant())
        print(mmat.DumpToString())
        mvec = math_Vector(1, 3)
        print(mvec.Norm())
        print(mvec.DumpToString())

    def export_file(self):
        write_step_file(self.compound, "./tmp/ODE1.stp")

    def display_Plane(self):
        self.show_axs_pln(scale=1.0)
        self.show()


if __name__ == '__main__':
    obj = ODE1()
    # obj.export_file()
    # obj.display_Plane()
