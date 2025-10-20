import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import argparse
from math import pi
from linecache import getline, clearcache, updatecache

basename = os.path.dirname(__file__)

sys.path.append(os.path.join("./"))
from base_occ import dispocc

import logging

logging.getLogger("matplotlib").setLevel(logging.ERROR)

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Pnt2d, gp_Ax2d, gp_Dir2d, gp_Circ2d
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Compound, TopoDS_Wire
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_WIRE
from OCC.Core.Geom import Geom_CylindricalSurface
from OCC.Core.GCE2d import GCE2d_MakeSegment, GCE2d_MakeCircle
from OCC.Core.Geom2dAPI import Geom2dAPI_Interpolate
from OCC.Core.TColgp import TColgp_Array1OfPnt2d, TColgp_HArray1OfPnt2d
from OCC.Core.BRepBuilderAPI import (
    BRepBuilderAPI_MakeEdge,
    BRepBuilderAPI_MakeWire,
    BRepBuilderAPI_MakeFace,
)
from OCC.Core.BRepAlgo import BRepAlgo_FaceRestrictor
from OCC.Core.ShapeAnalysis import ShapeAnalysis_Wire
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCCUtils.Common import vertex2pnt
from OCCUtils.Construct import make_box, make_line, make_wire, make_edge
from OCCUtils.Construct import make_plane, make_polygon
from OCCUtils.Construct import point_to_vector, vector_to_point
from OCCUtils.Construct import dir_to_vec, vec_to_dir


def uv_polygon2d_to_wire_on_surface(surface, uv_points, close=True):
    """
    UVパラメータ空間上の2D多角形を、指定したパラメトリックサーフェス上のWireに変換する。
    surface: Geom_Surface (例: Geom_CylindricalSurface, Geom_SphericalSurface, Geom_BSplineSurface)
    uv_points: [(u1, v1), (u2, v2), ...] のリスト
    close: Trueなら閉じる
    戻り値: TopoDS_Wire
    """
    edges = []
    n = len(uv_points)
    for i in range(n if close else n - 1):
        p1 = gp_Pnt2d(*uv_points[i])
        p2 = gp_Pnt2d(*uv_points[(i + 1) % n])
        seg2d = GCE2d_MakeSegment(p1, p2).Value()
        edge = BRepBuilderAPI_MakeEdge(
            seg2d, surface, seg2d.FirstParameter(), seg2d.LastParameter()
        ).Edge()
        edges.append(edge)
    wire_maker = BRepBuilderAPI_MakeWire()
    for edge in edges:
        wire_maker.Add(edge)
    return wire_maker.Wire()


def uv_circle_to_wire_on_surface(surface, center_uv, radius, num_segments=100):
    """
    UV空間上の円を指定したパラメトリックサーフェス上のWireに変換する。
    surface: Geom_Surface (例: Geom_CylindricalSurface)
    center_uv: (u, v) 円の中心のUV座標
    radius: 円の半径
    num_segments: 円を近似するセグメント数（デフォルト: 100）
    戻り値: TopoDS_Wire
    """
    # UV空間上の円を作成
    center = gp_Pnt2d(center_uv[0] % (2 * pi), center_uv[1])  # Uを周期的に扱う
    axis = gp_Ax2d(center, gp_Dir2d(1, 0))  # 円の中心と法線ベクトル
    circle = gp_Circ2d(axis, radius)
    circle_edge = GCE2d_MakeCircle(circle).Value()

    # 円をSurface上のエッジに変換
    edge = BRepBuilderAPI_MakeEdge(circle_edge, surface).Edge()

    # Wireを作成
    wire_maker = BRepBuilderAPI_MakeWire()
    wire_maker.Add(edge)
    return wire_maker.Wire()


def uv_spline_to_wire_on_surface(surface, uv_points):
    """
    UV空間上のSplineを指定したパラメトリックサーフェス上のWireに変換する。
    surface: Geom_Surface (例: Geom_CylindricalSurface)
    uv_points: [(u1, v1), (u2, v2), ...] のリスト
    戻り値: TopoDS_Wire
    """
    # UV空間上の点列をSplineに変換
    array = TColgp_Array1OfPnt2d(1, len(uv_points))
    for i, (u, v) in enumerate(uv_points, start=1):
        array.SetValue(i, gp_Pnt2d(u, v))  # Uを周期的に扱う
    harray = TColgp_HArray1OfPnt2d(array)
    curve_builder = Geom2dAPI_Interpolate(harray, True, 0.1e-6)
    curve_builder.Perform()
    spline_2d = curve_builder.Curve()

    # SplineをSurface上のエッジに変換
    edge = BRepBuilderAPI_MakeEdge(spline_2d, surface).Edge()

    # Wireを作成
    wire_maker = BRepBuilderAPI_MakeWire()
    wire_maker.Add(edge)
    return wire_maker.Wire()


def is_wire_closed_and_oriented(wire=TopoDS_Wire()):
    saw = ShapeAnalysis_Wire()
    wire = wire.Oriented(1)
    saw.Load(wire)
    api = saw.WireData().WireAPIMake()
    print("Wireの向き:", api.Orientation())
    print("Wireの向き付け可能性:", api.Orientable())
    print("Wireの閉じ:", api.Closed())


def extract_outer_and_inner_wires(face):
    outer_wire = None
    inner_wires = []
    exp = TopExp_Explorer(face, TopAbs_WIRE)
    while exp.More():
        wire = exp.Current()
        print("Wireの向き:", wire.Orientation())
        exp.Next()
    return outer_wire, inner_wires


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

    # 円筒面の作成
    cylinder = Geom_CylindricalSurface(gp_Ax3(), 5.0)
    face = BRepBuilderAPI_MakeFace(cylinder.Cylinder(), 0, 2 * pi, -10, 10).Face()

    # UV空間上の四角形
    uv_poly = [
        (2 * pi - pi / 3, 1),
        (2 * pi + pi / 4, 5),
        (2 * pi - pi / 4, -5),
    ]
    wire = uv_polygon2d_to_wire_on_surface(cylinder, uv_poly, True)

    # builder_face = BRepBuilderAPI_MakeFace(face)
    builder_face = BRepAlgo_FaceRestrictor()
    builder_face.Init(face, False, True)
    builder_face.Add(wire)
    builder_face.Perform()
    trimmed_face = builder_face.Current()

    obj.display.DisplayShape(face)
    obj.display.DisplayShape(wire)
    obj.display.DisplayShape(trimmed_face, color="RED")
    obj.ShowOCC()
