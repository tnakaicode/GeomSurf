import numpy as np
import matplotlib.pyplot as plt
import sys

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Lin, gp_Sphere
from OCC.Core.gp import gp_Mat, gp_XYZ
from OCC.Core.gp import gp_Trsf, gp_GTrsf
from OCC.Core.TColgp import TColgp_Array1OfPnt, TColgp_Array2OfPnt
from OCC.Core.Geom import Geom_BezierSurface, Geom_BSplineSurface, Geom_RectangularTrimmedSurface
from OCC.Core.GeomAPI import GeomAPI_PointsToBSplineSurface
from OCC.Core.GeomAbs import GeomAbs_C0, GeomAbs_C1, GeomAbs_C2
from OCC.Core.GeomAbs import GeomAbs_G1, GeomAbs_G2
from OCC.Core.GeomProjLib import geomprojlib_ProjectOnPlane
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.Core.BRepProj import BRepProj_Projection
from OCC.Core.BRep import BRep_Tool
from OCC.Core.GeomPlate import GeomPlate_MakeApprox, GeomPlate_BuildPlateSurface
from OCC.Core.GeomToStep import GeomToStep_MakeBoundedSurface
from OCC.Core.BRepFill import BRepFill_CurveConstraint
from OCC.Core.BOPAlgo import BOPAlgo_Tools
from OCC.Core.BOPTools import BOPTools_AlgoTools
from OCC.Extend.DataExchange import write_step_file
from OCCUtils.Construct import make_edge

from src.base_occ import dispocc, spl_face, set_loc


def bez_face(px, py, pz, axs=gp_Ax3()):
    nx, ny = px.shape
    pnt_2d = TColgp_Array2OfPnt(1, nx, 1, ny)
    for row in range(pnt_2d.LowerRow(), pnt_2d.UpperRow() + 1):
        for col in range(pnt_2d.LowerCol(), pnt_2d.UpperCol() + 1):
            i, j = row - 1, col - 1
            pnt = gp_Pnt(px[i, j], py[i, j], pz[i, j])
            pnt_2d.SetValue(row, col, pnt)
            #print (i, j, px[i, j], py[i, j], pz[i, j])

    surf = Geom_BezierSurface(pnt_2d)
    # api.Interpolate(pnt_2d)
    #surface = BRepBuilderAPI_MakeFace(curve, 1e-6)
    # return surface.Face()
    face = BRepBuilderAPI_MakeFace(surf, 1e-6).Face()
    face.Location(set_loc(gp_Ax3(), axs))
    return face


# https://www.opencascade.com/doc/occt-7.4.0/refman/html/class_b_rep_builder_a_p_i___make_face.html
if __name__ == '__main__':
    obj = dispocc()
    obj.show_axs_pln(obj.base_axs, scale=1)

    px = np.linspace(-1, 1, 10) * 5
    py = np.linspace(-1, 1, 10) * 5
    mesh = np.meshgrid(px, py)
    data = mesh[0]**2 / 10 + mesh[1]**2 / 20
    axis = gp_Ax3(gp_Pnt(0.5, 0.0, 0.0), gp_Dir(0, 0, 1))
    face = spl_face(*mesh, data, axs=axis)
    #face = bez_face(*mesh, data, axs=axis)
    trsf = face.Location().Transformation()
    surf = BRep_Tool.Surface(face)

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
    print(surf)
    api = BRepBuilderAPI_MakeFace(surf, 1.0)
    api.Add(bound_poly)
    bound_face = api.Face()
    print(api.IsDone())
    #api = BRepBuilderAPI_MakeFace(bound_face, bound_poly_1)
    #bound_face = api.Face()

    print(face.Location().Transformation())
    obj.display.DisplayShape(face, color="BLUE1", transparency=0.9)
    obj.display.DisplayShape(bound_face, color="RED", transparency=0.9)
    obj.display.DisplayShape(poly_0)
    obj.display.DisplayShape(poly_1)
    obj.display.DisplayShape(poly_2)
    obj.display.DisplayShape(bound_poly_0)
    obj.display.DisplayShape(bound_poly_1)
    obj.display.DisplayShape(bound_poly_2)
    obj.display.DisplayShape(bound_poly)

    #write_step_file(face, obj.tempname + "_org.stp")
    #write_step_file(bound_face, obj.tempname + ".stp")

    obj.ShowOCC()

    # ハンドルされない例外が 0x00007FFA7A40E8AF (TKSTEP.dll) で発生しました(python.exe 内): 0xC0000005: 場所 0x0000000000000000 の読み取り中にアクセス違反が発生しました。
    # 00007FFA7A40E8AF  mov         rax,qword ptr [rcx]
