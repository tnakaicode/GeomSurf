import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import argparse
from linecache import getline, clearcache, updatecache

basename = os.path.dirname(__file__)

sys.path.append(os.path.join("./"))
from base_occ import dispocc

import logging

logging.getLogger("matplotlib").setLevel(logging.ERROR)

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3, gp_Trsf
from OCC.Core.gp import gp_Circ, gp_Elips, gp_Parab, gp_Hypr, gp_Pln, gp_Lin
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_ThruSections
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Fuse
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeShell, BRepBuilderAPI_MakeSolid
from OCC.Core.BRep import BRep_Builder
from OCC.Core.TopoDS import TopoDS_Compound
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Compound, TopoDS_Shell
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCCUtils.Common import vertex2pnt
from OCCUtils.Construct import make_box, make_line, make_wire, make_edge
from OCCUtils.Construct import make_plane, make_polygon, make_face
from OCCUtils.Construct import point_to_vector, vector_to_point
from OCCUtils.Construct import dir_to_vec, vec_to_dir


def make_elps(axs=gp_Ax3(), r0=1, r1=1, r2=1):
    """楕円を作成する関数"""
    p0 = gp_Pnt(0, 0, 0)
    trf = gp_Trsf()
    trf.SetTransformation(axs)
    if r1 > r0:
        ax0 = gp_Ax2(p0, gp_Dir(1, 0, 0), gp_Dir(0, 0, 1))
        ax0.Transform(trf)
        ep0 = gp_Elips(ax0, r1, r0)
    else:
        ax0 = gp_Ax2(p0, gp_Dir(1, 0, 0), gp_Dir(0, 1, 0))
        ax0.Transform(trf)
        ep0 = gp_Elips(ax0, r0, r1)

    if r1 > r2:
        ax1 = gp_Ax2(p0, gp_Dir(0, 1, 0), gp_Dir(0, 0, 1))
        ax1.Transform(trf)
        ep1 = gp_Elips(ax1, r1, r2)
    else:
        ax1 = gp_Ax2(p0, gp_Dir(0, 1, 0), gp_Dir(1, 0, 0))
        ax1.Transform(trf)
        ep1 = gp_Elips(ax1, r2, r1)

    if r2 > r0:
        ax2 = gp_Ax2(p0, gp_Dir(0, 0, 1), gp_Dir(1, 0, 0))
        ax2.Transform(trf)
        ep2 = gp_Elips(ax2, r2, r0)
    else:
        ax2 = gp_Ax2(p0, gp_Dir(0, 0, 1), gp_Dir(0, 1, 0))
        ax2.Transform(trf)
        ep2 = gp_Elips(ax2, r0, r2)

    return ep0, ep1, ep2


if __name__ == "__main__":
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument(
        "--pxyz", dest="pxyz", default=[0.0, 0.0, 0.0], type=float, nargs=3
    )
    opt = parser.parse_args()
    print(opt)

    obj = dispocc(touch=True)

    axs = gp_Ax3()
    axs.Translate(gp_Vec(1, 2, 3))
    axs.Rotate(gp_Ax1(axs.Location(), axs.XDirection()), np.pi / 10)
    r0 = 16
    r1 = 15
    r2 = 12
    ep0, ep1, ep2 = make_elps(axs, r0, r1, r2)

    el0 = make_edge(ep0, -2 * np.pi / 3, 2 * np.pi / 3)
    el1 = make_edge(ep1, -2 * np.pi / 3, 2 * np.pi / 3)
    el2 = make_edge(ep2, -2 * np.pi / 3, 2 * np.pi / 3)

    # 楕円を基にした面を作成
    face0 = make_face(make_wire([el0]))
    face1 = make_face(make_wire([el1]))
    face2 = make_face(make_wire([el2]))

    # 面を結合して閉じた立体を作成
    fused_shape = BRepAlgoAPI_Fuse(face0, face1).Shape()
    fused_shape = BRepAlgoAPI_Fuse(fused_shape, face2).Shape()

    # 面を結合してシェルを作成
    builder = BRep_Builder()
    shell = TopoDS_Shell()
    builder.MakeShell(shell)
    builder.Add(shell, face0)
    builder.Add(shell, face1)
    builder.Add(shell, face2)

    # 楕円を結合して閉じた立体を作成
    builder = BRepOffsetAPI_ThruSections(True, False, 1e-6)
    builder.AddWire(make_wire([el0]))
    builder.AddWire(make_wire([el1]))
    builder.AddWire(make_wire([el2]))
    builder.Build()

    # 作成した形状を取得
    solid = builder.Shape()

    obj.display.DisplayShape(el0, color="RED")
    obj.display.DisplayShape(el1, color="GREEN")
    obj.display.DisplayShape(el2, color="BLUE1")
    # obj.display.DisplayShape(shell, transparency=0.7, color="ORANGE")
    # obj.display.DisplayShape(solid, transparency=0.7, color="ORANGE")
    # obj.display.DisplayShape(fused_shape, transparency=0.7, color="ORANGE")

    obj.ShowOCC()
