import numpy as np
import matplotlib.pyplot as plt
import sys
import time
import os

sys.path.append(os.path.join('../'))
from src.base_occ import dispocc, which

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Pln, gp_Circ, gp_Sphere, gp_Elips
from OCC.Core.Geom import Geom_Plane, Geom_Circle, Geom_SphericalSurface, Geom_Ellipse
from OCC.Core.BRep import BRep_PolygonOnSurface, BRep_PolygonOnTriangulation, BRep_PointOnSurface
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_NurbsConvert, BRepBuilderAPI_ModifyShape
from OCC.Core.BRepOffset import BRepOffset_MakeOffset, BRepOffset_Skin, BRepOffset_Interval
from OCCUtils.Construct import make_box, make_line, make_wire, make_edge, make_face


if __name__ == '__main__':
    obj = dispocc(touch=True)

    axs = gp_Ax3(gp_Pnt(0, 0, 1.0), gp_Dir(1, 0.5, 1))
    ply = obj.make_EllipWire(rxy=[10, 15], axs=axs)
    pln = Geom_Plane(axs)
    sph = gp_Sphere(axs, 30)
    box = make_box(gp_Pnt(-100, -100, -100), gp_Pnt(100, 100, 100))
    builder = BRepBuilderAPI_NurbsConvert(box)
    # builder.Build()
    # obj.display.DisplayShape(builder.Shape())

    plan_face = make_face(pln, ply)
    ball_face = make_face(sph)
    elip_proj = obj.proj_rim_pln(ply, ball_face, axs, idx=1)
    ball_face = make_face(ball_face, elip_proj)

    print(box)
    print(builder.Shape())
    print(axs.Location())
    obj.display.DisplayShape(plan_face)
    obj.display.DisplayShape(ball_face, transparency=0.5, color="BLUE1")
    obj.display.DisplayShape(elip_proj)
    obj.show_axs_pln(scale=10)
    obj.ShowOCC()
