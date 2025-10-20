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
from base_occ import dispocc, spl_face

import logging

logging.getLogger("matplotlib").setLevel(logging.ERROR)

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Pnt2d, gp_Circ2d, gp_Ax2d, gp_Dir2d
from OCC.Core.gp import gp_Circ
from OCC.Core.GC import GC_MakeCircle, GC_MakeEllipse
from OCC.Core.TColgp import TColgp_Array2OfPnt
from OCC.Core.TColgp import TColgp_Array1OfPnt2d, TColgp_HArray1OfPnt2d
from OCC.Core.GCE2d import GCE2d_MakeSegment, GCE2d_MakeCircle
from OCC.Core.GeomAbs import GeomAbs_G2
from OCC.Core.Geom import Geom_BSplineSurface
from OCC.Core.Geom import Geom_SphericalSurface
from OCC.Core.Geom import Geom_BSplineSurface
from OCC.Core.Geom import Geom_CylindricalSurface
from OCC.Core.Geom2dAPI import Geom2dAPI_PointsToBSpline, Geom2dAPI_Interpolate
from OCC.Core.GeomAPI import GeomAPI_PointsToBSplineSurface
from OCC.Core.GeomAPI import GeomAPI_ProjectPointOnSurf, GeomAPI_PointsToBSpline
from OCC.Core.BRepBuilderAPI import (
    BRepBuilderAPI_MakeEdge,
    BRepBuilderAPI_MakeWire,
    BRepBuilderAPI_MakeFace,
    BRepBuilderAPI_Sewing,
)
from OCC.Core.BRepAlgo import BRepAlgo_FaceRestrictor
from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRep import BRep_Builder
from OCC.Core.BRepProj import BRepProj_Projection
from OCC.Core.BRepAdaptor import BRepAdaptor_CompCurve
from OCC.Core.TopAbs import TopAbs_WIRE, TopAbs_Orientation, TopAbs_EDGE
from OCC.Core.TopoDS import TopoDS_Wire, topods
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Compound
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.ShapeAnalysis import ShapeAnalysis_Wire
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCCUtils.Common import vertex2pnt
from OCCUtils.Construct import make_box, make_line, make_wire, make_edge
from OCCUtils.Construct import make_plane, make_polygon
from OCCUtils.Construct import point_to_vector, vector_to_point
from OCCUtils.Construct import dir_to_vec, vec_to_dir


def polygon2d_to_wire_on_surface(surface, uv_points, close=True):
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


def sample_points_on_wire(wire: TopoDS_Wire, n: int):
    # ワイヤをパラメータ化して等間隔サンプリング
    compcurve = BRepAdaptor_CompCurve(wire)
    first = compcurve.FirstParameter()
    last = compcurve.LastParameter()
    params = np.linspace(first, last, n, endpoint=False)
    points = []
    for u in params:
        pnt = compcurve.Value(u)
        points.append(pnt)
    return points


def close_projected_wire(wire, face, direction):
    # 投影方向を指定（例: Z方向）
    projector = BRepProj_Projection(wire, face, gp_Dir(0, 0, 1))
    projected_wire = projector.Current()
    # 2. Check if the projected wire is closed
    # analyzer = ShapeAnalysis_Wire(projected_wire, face, 1e-7)
    # if analyzer.IsClosed():
    return projected_wire  # Already closed

    # 3. Get start and end points of the open wire
    edges = []
    exp = TopExp_Explorer(projected_wire, TopAbs_EDGE)
    while exp.More():
        edges.append(topods.Edge(exp.Current()))
        exp.Next()
    if not edges:
        raise RuntimeError("No edges in projected wire")

    # Get start and end points
    from OCC.Core.BRep import BRep_Tool

    first_edge = edges[0]
    last_edge = edges[-1]
    first_vert = BRep_Tool().Pnt(BRep_Tool().Vertex(first_edge, True))
    last_vert = BRep_Tool().Pnt(BRep_Tool().Vertex(last_edge, False))

    # 4. Project these points onto the face (to ensure they are on the surface)
    surf = BRep_Tool().Surface(face)
    proj1 = GeomAPI_ProjectPointOnSurf(first_vert, surf)
    proj2 = GeomAPI_ProjectPointOnSurf(last_vert, surf)
    p1 = proj1.NearestPoint()
    p2 = proj2.NearestPoint()

    # 5. Create a BSpline curve on the face between these points
    interp = GeomAPI_PointsToBSpline([p1, p2])
    edge_interp = BRepBuilderAPI_MakeEdge(interp.Curve(), p1, p2).Edge()

    # 6. Build a new closed wire
    wire_builder = BRepBuilderAPI_MakeWire()
    for e in edges:
        wire_builder.Add(e)
    wire_builder.Add(edge_interp)
    closed_wire = wire_builder.Wire()
    return closed_wire


if __name__ == "__main__":
    print(
        "This module provides a function to close a projected wire onto a face in OpenCASCADE."
    )
    obj = dispocc(touch=True)
    axs = gp_Ax3()

    # 円筒面と球面を作成
    cylinder = Geom_CylindricalSurface(gp_Ax3(), 10.0)
    cylinder_face = BRepBuilderAPI_MakeFace(
        cylinder.Cylinder(), 0, 2 * pi, -10, 10
    ).Face()
    obj.display.DisplayShape(cylinder_face, transparency=0.8, color="BLUE1")

    sphere = Geom_SphericalSurface(gp_Ax3(), 8.0)
    sphere_face = BRepBuilderAPI_MakeFace(
        sphere.Sphere(), 0, 2 * pi, -pi / 2, pi / 2
    ).Face()
    obj.display.DisplayShape(sphere_face, transparency=0.8, color="RED")

    # 円筒面上の閉じたワイヤ（UV空間の多角形）
    uv_poly = [
        (0, pi / 3),
        (pi / 4, pi / 3),
        (-pi / 10, -pi / 10),
    ]
    wire = polygon2d_to_wire_on_surface(sphere, uv_poly, close=True)
    obj.display.DisplayShape(wire, color="GREEN")

    wire = BRepBuilderAPI_MakeWire(
        BRepBuilderAPI_MakeEdge(
            gp_Circ(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0)), 5)
        ).Edge()
    ).Wire()
    obj.display.DisplayShape(wire, color="GREEN")

    px = np.linspace(-1, 1, 100) * 10 + 0.50
    py = np.linspace(-1, 1, 200) * 10 - 0.50
    mesh = np.meshgrid(px, py)
    data = mesh[0] ** 2 / 500 + mesh[1] ** 2 / 200 + 8
    surf = spl_face(*mesh, data)
    uv_poly = [
        [-5, 5, 0],
        [5, 15, 0],
        [5, -5, 0],
    ]
    wire = make_polygon([gp_Pnt(*xyz) for xyz in uv_poly], closed=True)
    obj.display.DisplayShape(wire, color="GREEN")
    obj.display.DisplayShape(surf, transparency=0.5)

    # ワイヤを球面に投影し、open wireなら閉じる
    closed_wire = close_projected_wire(wire, surf, gp_Dir(0, 0, 1))
    obj.display.DisplayShape(closed_wire, color="YELLOW")

    # サンプリング点を可視化
    # pts = sample_points_on_wire(closed_wire, 50)
    # for p in pts:
    #    obj.display.DisplayShape(p, color="BLACK")

    obj.ShowOCC()
