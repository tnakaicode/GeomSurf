import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import argparse
from linecache import getline, clearcache

basename = os.path.dirname(__file__)

sys.path.append(os.path.join("./"))
from src.base_occ import dispocc

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Circ, gp_Parab, gp_Parab2d
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Compound
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_ThruSections, BRepOffsetAPI_MakeOffset
from OCC.Core.BRepOffset import BRepOffset_MakeOffset, BRepOffset_Interval
from OCC.Core.BRepOffset import BRepOffset_Skin, BRepOffset_Pipe, BRepOffset_RectoVerso
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakePrism
from OCC.Core.GeomAbs import GeomAbs_C0
from OCC.Core.GeomAbs import GeomAbs_Intersection, GeomAbs_Arc, GeomAbs_Tangent
from OCC.Core.Geom import Geom_Parabola, Geom_Surface, Geom_Curve, Geom_TrimmedCurve
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCCUtils.Common import vertex2pnt
from OCCUtils.Construct import make_box, make_line, make_wire, make_edge
from OCCUtils.Construct import make_plane, make_polygon
from OCCUtils.Construct import point_to_vector, vector_to_point
from OCCUtils.Construct import dir_to_vec, vec_to_dir
from OCCUtils.Construct import make_wire, make_edge, make_polygon, make_face


def make_parabola(axs=gp_Ax3(), dst=100, width=[-50, 50], hight=[-50, 50], xyz="z"):
    """make Parabola Surface

    Args:
        axs (gp_Ax3, optional): Defaults to gp_Ax3().
        dst (float, optional): Parabola focus length. Defaults to 100.
        width (list, optional): Defaults to [-50, 50].
        hight (list, optional): Defaults to [-50, 50].
        xyz (str, optional): Defaults to "z".

    Returns:
        TopoDS_Shape : Parabola surface
    """
    if xyz == "z":
        axis = gp_Ax3(axs.Location(), axs.Direction(), axs.XDirection())
    elif xyz == "y":
        axis = gp_Ax3(axs.Location(), axs.YDirection(), axs.Direction())
    elif xyz == "x":
        axis = gp_Ax3(axs.Location(), axs.XDirection(), axs.Direction())
    else:
        axis = gp_Ax3(axs.Ax2())
    ax1 = dispocc.prop_axs(None, axis, hight[0], "z")
    ax2 = dispocc.prop_axs(None, axis, hight[1], "z")
    crv1 = Geom_TrimmedCurve(Geom_Parabola(ax1.Ax2(), dst), width[0], width[1])
    crv2 = Geom_TrimmedCurve(Geom_Parabola(ax2.Ax2(), dst), width[0], width[1])

    api = BRepOffsetAPI_ThruSections(False, True, 1.0e-6)
    api.AddWire(make_wire(make_edge(crv1)))
    api.AddWire(make_wire(make_edge(crv2)))
    api.Build()
    return api.Shape()


if __name__ == '__main__':
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument("--pxyz", dest="pxyz",
                        default=[0.0, 0.0, 0.0], type=float, nargs=3)
    opt = parser.parse_args()
    print(opt)

    obj = dispocc(touch=True)

    axs = gp_Ax3()
    para = Geom_Parabola(axs.Ax2(), 100)
    curv = Geom_TrimmedCurve(para, -100, 100)

    focus = 200
    face = make_parabola(xyz="z", width=[10, 200], hight=[-100,100], dst=focus)
    obj.display.DisplayShape(face, color="BLUE1")
    obj.show_axs_pln(gp_Ax3(gp_Pnt(focus, 0, 0), gp_Dir(1, 0, 0)), scale=50)

    beams = [
        gp_Ax3(gp_Pnt(100, 120, -10), gp_Dir(-0.9, 0.0, 0.0)),
        gp_Ax3(gp_Pnt(100, 150, 10), gp_Dir(-0.9, -0.0, 0.0)),
        gp_Ax3(gp_Pnt(100, 120, -10), gp_Dir(-0.9, 0.1, 0.0)),
        gp_Ax3(gp_Pnt(100, 150, 10), gp_Dir(-0.9, -0.3, 0.0))
    ]
    for i, beam in enumerate(beams):
        beam1 = obj.reflect_beam(face, beam)
        beam2 = obj.prop_axs(beam1, scale=300)
        pts = [b.Location() for b in [beam, beam1, beam2]]
        lin = make_polygon(pts)
        obj.display.DisplayShape(lin, color=obj.colors[i])
        obj.show_axs_pln(beam, scale=10, name="beam" + str(i))

    face = make_parabola(xyz="y")
    obj.display.DisplayShape(face, color="GREEN")

    face = make_parabola(xyz="x")
    obj.display.DisplayShape(face, color="RED")

    obj.show_axs_pln()
    obj.display.DisplayShape(curv)
    obj.ShowOCC()
