import numpy as np
import matplotlib.pyplot as plt
import sys

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Lin, gp_Sphere
from OCC.Core.gp import gp_Mat, gp_XYZ
from OCC.Core.gp import gp_Trsf, gp_GTrsf
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.Core.BRepProj import BRepProj_Projection
from OCC.Extend.DataExchange import write_step_file
from OCCUtils.Construct import make_edge

from base import plotocc, spl_face

# https://www.opencascade.com/doc/occt-7.4.0/refman/html/class_b_rep_builder_a_p_i___make_face.html
if __name__ == '__main__':
    obj = plotocc()
    obj.show_axs_pln(obj.base_axs, scale=1)

    px = np.linspace(-1, 1, 10) * 5
    py = np.linspace(-1, 1, 10) * 5
    mesh = np.meshgrid(px, py)
    surf = mesh[0]**2 / 10 + mesh[1]**2 / 20
    axis = gp_Ax3(gp_Pnt(0.5, 0.0, 0.0), gp_Dir(0, 0, 1))
    face = spl_face(*mesh, surf, axs=axis)
    poly = obj.make_EllipWire(rxy=[1.1, 1.0], axs=axis)
    poly = obj.make_PolyWire(num=6, axs=axis)
    poly = obj.make_PolyWire(num=10, axs=axis)
    proj = BRepProj_Projection(poly, face, axis.Direction())

    bound_poly = proj.Current()
    print(bound_poly)
    api = BRepBuilderAPI_MakeFace(face, bound_poly)
    bound_face = api.Face()

    print(face.Location().Transformation())
    obj.display.DisplayShape(bound_face)
    obj.display.DisplayShape(poly)
    obj.display.DisplayShape(bound_poly)

    write_step_file(face, "./tmp/CutSurf_org.stp")
    write_step_file(bound_face, "./tmp/CutSurf.stp")

    obj.show()
