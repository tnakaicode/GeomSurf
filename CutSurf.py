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

from src.base import plotocc, spl_face

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
    trsf = face.Location().Transformation()

    axis_0 = axis.Transformed(trsf)
    axis_0.Translate(gp_Pnt(0, 0, 0), gp_Pnt(2, 0, 0))
    poly_0 = obj.make_EllipWire(rxy=[1.1, 1.0], axs=axis_0)
    proj = BRepProj_Projection(poly_0, face, axis.Direction())
    bound_poly_0 = proj.Current()

    axis_1 = axis.Transformed(trsf)
    axis_1.Translate(gp_Pnt(0, 0, 0), gp_Pnt(-2, 0, 0))
    poly_1 = obj.make_PolyWire(num=6, axs=axis_1)
    proj = BRepProj_Projection(poly_1, face, axis.Direction())
    bound_poly_1 = proj.Current()

    axis_2 = axis.Transformed(trsf)
    axis_2.Translate(gp_Pnt(0, 0, 0), gp_Pnt(0, 2, 0))
    poly_2 = obj.make_PolyWire(num=10, axs=axis_2)
    proj = BRepProj_Projection(poly_2, face, axis.Direction())
    bound_poly_2 = proj.Current()

    proj = BRepProj_Projection(poly_0, face, axis.Direction())
    bound_poly = proj.Current()
    print(bound_poly)
    api = BRepBuilderAPI_MakeFace(face, bound_poly_0)
    bound_face = api.Face()
    api = BRepBuilderAPI_MakeFace(bound_face, bound_poly_1)
    bound_face = api.Face()

    print(face.Location().Transformation())
    obj.display.DisplayShape(bound_face)
    obj.display.DisplayShape(poly_0)
    obj.display.DisplayShape(poly_1)
    obj.display.DisplayShape(poly_2)
    obj.display.DisplayShape(bound_poly_0)
    obj.display.DisplayShape(bound_poly_1)
    obj.display.DisplayShape(bound_poly_2)
    obj.display.DisplayShape(bound_poly)

    write_step_file(face, obj.tempname + "_org.stp")
    write_step_file(bound_face, obj.tempname + ".stp")

    obj.show()
