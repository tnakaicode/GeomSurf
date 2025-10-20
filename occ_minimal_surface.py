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
from OCC.Core.gp import gp_Pnt2d, gp_Circ2d, gp_Ax2d, gp_Dir2d
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
from OCC.Core.BRepBuilderAPI import (
    BRepBuilderAPI_MakeEdge,
    BRepBuilderAPI_MakeWire,
    BRepBuilderAPI_MakeFace,
    BRepBuilderAPI_Sewing,
)
from OCC.Core.BRepAlgo import BRepAlgo_FaceRestrictor
from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRep import BRep_Builder
from OCC.Core.BRepAdaptor import BRepAdaptor_CompCurve
from OCC.Core.TopAbs import TopAbs_WIRE, TopAbs_Orientation
from OCC.Core.TopoDS import TopoDS_Wire
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


def minimal_surface_from_wires_minimal_like(
    wire1: TopoDS_Wire,
    wire2: TopoDS_Wire,
    n: int = 100,
    m: int = 30,
    shrink_ratio=0.5,
):
    """Create a minimal-like surface between two wires.
    This function generates a minimal-like surface by sampling points on two wires
    and creating a grid of points that shrink towards the center linearly.
    It uses the hyperbolic cosine function to create a smooth transition between the two wires.
    The resulting surface is approximated by triangular faces.
    It also calculates the area of the generated surface.
    This is useful for creating a minimal surface approximation between two curves in 3D space.

    Args:
        wire1 (TopoDS_Wire):
        wire2 (TopoDS_Wire):
        n (int, optional): Defaults to 100.
        m (int, optional): \Defaults to 30.
        shrink_ratio (float, optional): Defaults to 0.5.
    Returns:
        _type_: _description_
    """
    print("Wire1 Direction:", wire1.Orientation())
    print("Wire2 Direction:", wire2.Orientation())
    pts1 = sample_points_on_wire(wire1, n)
    pts2 = sample_points_on_wire(wire2, n)
    grid = []
    for k in range(m):
        t = k / (m - 1)
        s = (np.cosh((t - 0.5) * 2) - np.cosh(1)) / (np.cosh(0) - np.cosh(1))
        ring = []
        center = (1 - t) * np.mean(
            [[p.X(), p.Y(), p.Z()] for p in pts1], axis=0
        ) + t * np.mean([[p.X(), p.Y(), p.Z()] for p in pts2], axis=0)
        for i in range(n):
            p1 = np.array([pts1[i].X(), pts1[i].Y(), pts1[i].Z()])
            p2 = np.array([pts2[i].X(), pts2[i].Y(), pts2[i].Z()])
            p = (1 - t) * p1 + t * p2
            p = p + (center - p) * s * shrink_ratio
            ring.append(gp_Pnt(*p))
        grid.append(ring)
    edges = []
    area = 0.0
    for k in range(m - 1):
        for i in range(n):
            i_next = (i + 1) % n
            p1 = grid[k][i]
            p2 = grid[k + 1][i]
            p3 = grid[k + 1][i_next]
            p4 = grid[k][i_next]
            # 三角形1
            edges.append(
                BRepBuilderAPI_MakeWire(
                    BRepBuilderAPI_MakeEdge(p1, p2).Edge(),
                    BRepBuilderAPI_MakeEdge(p2, p3).Edge(),
                    BRepBuilderAPI_MakeEdge(p3, p1).Edge(),
                ).Wire()
            )
            # 三角形2
            edges.append(
                BRepBuilderAPI_MakeWire(
                    BRepBuilderAPI_MakeEdge(p1, p3).Edge(),
                    BRepBuilderAPI_MakeEdge(p3, p4).Edge(),
                    BRepBuilderAPI_MakeEdge(p4, p1).Edge(),
                ).Wire()
            )
            # 面積計算（2三角形分）
            # 三角形1
            v1 = np.array([p2.X() - p1.X(), p2.Y() - p1.Y(), p2.Z() - p1.Z()])
            v2 = np.array([p3.X() - p1.X(), p3.Y() - p1.Y(), p3.Z() - p1.Z()])
            area += 0.5 * np.linalg.norm(np.cross(v1, v2))
            # 三角形2
            v1 = np.array([p3.X() - p1.X(), p3.Y() - p1.Y(), p3.Z() - p1.Z()])
            v2 = np.array([p4.X() - p1.X(), p4.Y() - p1.Y(), p4.Z() - p1.Z()])
            area += 0.5 * np.linalg.norm(np.cross(v1, v2))
    print(f"Total area: {area:.4f}")
    sewing = BRepBuilderAPI_Sewing()
    for w in edges:
        face = BRepBuilderAPI_MakeFace(w)
        sewing.Add(face.Face())
    sewing.Perform()
    shell = sewing.SewedShape()
    return shell, area


def minimal_surface_on_single_wire(wire: TopoDS_Wire, n: int = 100):
    """
    1つのWireに張る極小曲面近似（三角形ファン or 平面）を生成
    平面上のWireなら平面、立体的なWireなら重心ファン
    Args:
        wire (TopoDS_Wire): 閉じたワイヤ
        n (int): サンプリング点数
        return_area (bool): Trueなら面積も返す
    Returns:
        shell (TopoDS_Shape): 曲面シェル
        area (float, optional): 面積
    """
    pts = sample_points_on_wire(wire, n)
    # 平面判定: 最小二乗平面を求め、全点の距離が小さければ平面とみなす
    coords = np.array([[p.X(), p.Y(), p.Z()] for p in pts])
    centroid = coords.mean(axis=0)
    uu, dd, vv = np.linalg.svd(coords - centroid)
    normal = vv[-1]
    # 各点から平面への距離
    dists = np.abs((coords - centroid) @ normal)
    max_dist = np.max(dists)
    tol = 1e-6  # 許容誤差
    if max_dist < tol:
        # 平面上なら、その平面に張る
        from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace

        face = BRepBuilderAPI_MakeFace(wire).Face()
        # 面積計算
        area = 0.0
        for i in range(n):
            i_next = (i + 1) % n
            p1 = coords[i]
            p2 = coords[i_next]
            v1 = p2 - centroid
            v2 = coords[i] - centroid
            # 三角形面積
            area += 0.5 * np.linalg.norm(
                np.cross(p2 - centroid, coords[i_next] - centroid)
            )
        return face, area
    else:
        # 立体的なWire→重心ファン
        center = centroid
        center_pnt = gp_Pnt(*center)
        edges = []
        area = 0.0
        for i in range(n):
            i_next = (i + 1) % n
            p1 = pts[i]
            p2 = pts[i_next]
            # 三角形ワイヤ
            wire_tri = BRepBuilderAPI_MakeWire(
                BRepBuilderAPI_MakeEdge(p1, p2).Edge(),
                BRepBuilderAPI_MakeEdge(p2, center_pnt).Edge(),
                BRepBuilderAPI_MakeEdge(center_pnt, p1).Edge(),
            ).Wire()
            edges.append(wire_tri)
            # 面積計算
            v1 = np.array([p2.X() - p1.X(), p2.Y() - p1.Y(), p2.Z() - p1.Z()])
            v2 = np.array(
                [
                    center_pnt.X() - p1.X(),
                    center_pnt.Y() - p1.Y(),
                    center_pnt.Z() - p1.Z(),
                ]
            )
            area += 0.5 * np.linalg.norm(np.cross(v1, v2))
        # Sewingでシェル化
        sewing = BRepBuilderAPI_Sewing()
        for w in edges:
            face = BRepBuilderAPI_MakeFace(w)
            sewing.Add(face.Face())
        sewing.Perform()
        shell = sewing.SewedShape()
        return shell, area


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

    # 円筒面
    cylinder = Geom_CylindricalSurface(gp_Ax3(), 6.0)
    cylinder_face = BRepBuilderAPI_MakeFace(
        cylinder.Cylinder(), 0, 2 * pi, -10, 10
    ).Face()
    obj.display.DisplayShape(cylinder_face, transparency=0.9, color="BLUE1")

    # 球面
    sphere = Geom_SphericalSurface(gp_Ax3(), 8.0)
    sphere_face = BRepBuilderAPI_MakeFace(
        sphere.Sphere(), 0, 2 * pi, -pi / 2, pi / 2
    ).Face()
    obj.display.DisplayShape(sphere_face, transparency=0.9, color="RED")

    # UV空間上の四角形
    uv_poly1 = [
        (pi, 5),
        (pi + pi / 4, 7),
        (pi + pi / 4 + pi / 11, 1),
        (pi - pi / 10, -1),
    ]
    wire1 = polygon2d_to_wire_on_surface(cylinder, uv_poly1, True)
    wire1_pts = sample_points_on_wire(wire1, 100)
    obj.display.DisplayShape(wire1)
    for p in wire1_pts:
        obj.display.DisplayShape(p)

    # UV空間上の四角形
    uv_poly2 = [
        (0, pi / 3),
        (pi / 4 + pi / 11, pi / 3 + pi / 8),
        (-pi / 10, -pi / 10),
    ]
    wire2 = polygon2d_to_wire_on_surface(sphere, uv_poly2, True)
    wire2_pts = sample_points_on_wire(wire2, 100)
    obj.display.DisplayShape(wire2)
    for p in wire2_pts:
        obj.display.DisplayShape(p)

    # 極小曲面近似
    shell, area = minimal_surface_from_wires_minimal_like(
        wire1, wire2.Reversed(), n=100, m=100, shrink_ratio=0.4
    )
    obj.display.DisplayShape(shell, transparency=0.5)

    # shell1, area = minimal_surface_on_single_wire(wire1, n=100)
    # print(f"面積: {area:.4f}")
    # obj.display.DisplayShape(shell1, update=True)

    obj.ShowOCC()
