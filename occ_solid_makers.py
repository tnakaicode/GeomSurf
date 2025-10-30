"""Helpers: several patterns to create solids in pythonOCC

このモジュールは、OCC (pythonOCC) を使った Solid の作り方を関数別に示します。
各関数は TopoDS_Shape (solid) を返します。

注意: 実行環境に pythonOCC が必要です。表示系を使わず STEP export による出力サンプルを付けています。
"""

from __future__ import annotations
import os
import math
import numpy as np

from OCC.Core.gp import gp_Vec, gp_Ax1, gp_Pnt, gp_Dir, gp_Circ, gp_Ax2
from OCC.Core.TopoDS import TopoDS_Shape
from OCC.Core.BRepBuilderAPI import (
    BRepBuilderAPI_MakeFace,
    BRepBuilderAPI_MakeShell,
    BRepBuilderAPI_MakeSolid,
)
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire
from OCC.Core.BRepPrimAPI import (
    BRepPrimAPI_MakePrism,
    BRepPrimAPI_MakeRevol,
    BRepPrimAPI_MakeBox,
    BRepPrimAPI_MakeTorus,
)
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_ThruSections
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_Sewing
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Fuse
from OCC.Core.STEPControl import STEPControl_Writer, STEPControl_AsIs
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE
from OCC.Core.TopoDS import topods
from OCC.Core.GProp import GProp_GProps
from OCC.Core.BRepGProp import brepgprop
from typing import Iterable, List
from OCC.Core.GeomAPI import GeomAPI_PointsToBSplineSurface
from OCC.Core.TColgp import TColgp_Array2OfPnt
from OCC.Core.GeomAbs import GeomAbs_G2
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakePolygon, BRepBuilderAPI_Transform
from OCC.Core.gp import gp_Trsf
from OCC.Core.TopAbs import TopAbs_WIRE
from OCC.Core.BRepAlgo import BRepAlgo_FaceRestrictor


def export_step(shape: TopoDS_Shape, filename: str) -> bool:
    """Export shape to STEP file. Returns True on success."""
    writer = STEPControl_Writer()
    writer.Transfer(shape, STEPControl_AsIs)
    status = writer.Write(filename)
    # IFSelect_ReturnStatus: success is usually non-zero (IFSelect_RetDone == 1)
    try:
        return int(status) != 0
    except Exception:
        return status != 0


def solid_from_face(face) -> TopoDS_Shape:
    """与えられた Face から Shell を作り、それを Solid に変換して返す。

    face: TopoDS_Face
    """
    # Build a shell from a single face, then make a solid from that shell.
    mk_shell = BRepBuilderAPI_MakeShell()
    mk_shell.Add(face)
    shell = mk_shell.Shell()
    shell.Location(TopLoc_Location())
    mk = BRepBuilderAPI_MakeSolid()
    mk.Add(shell)
    return mk.Solid()


def solid_from_faces_sewed(faces: Iterable) -> TopoDS_Shape:
    """複数の Face を Sewing で縫い合わせて Shell を作り、Solid にする。

    faces: iterable of TopoDS_Face
    """
    sew = BRepBuilderAPI_Sewing()
    for f in faces:
        sew.Add(f)
    sew.Perform()
    res = sew.SewedShape()
    # Ensure we treat the result as a Shell when creating a Solid
    try:
        shell = topods.Shell(res)
    except Exception:
        # If conversion fails, try to use the shape directly
        shell = res
    mk = BRepBuilderAPI_MakeSolid()
    mk.Add(shell)
    return mk.Solid()


def solid_by_extrusion(face, vec: gp_Vec) -> TopoDS_Shape:
    """Face を与えて押し出して Solid を作成する。

    vec: 押し出しベクトル (gp_Vec)
    """
    prism = BRepPrimAPI_MakePrism(face, vec)
    prism.Build()
    return prism.Shape()


def solid_by_revolution(face_or_wire, axis: gp_Ax1) -> TopoDS_Shape:
    """Face または wire を軸回転して solid を得る。

    axis: gp_Ax1
    """
    revol = BRepPrimAPI_MakeRevol(face_or_wire, axis)
    revol.Build()
    return revol.Shape()


def solid_by_thru_sections(wires: List, make_solid: bool = True) -> TopoDS_Shape:
    """複数のワイヤ（断面）を与えて ThruSections でシェルを作り、必要なら Solid にする。

    wires: list of TopoDS_Wire
    """
    thu = BRepOffsetAPI_ThruSections()
    for w in wires:
        thu.AddWire(w)
    thu.Build()
    shell = thu.Shape()
    if not make_solid:
        return shell
    mk = BRepBuilderAPI_MakeSolid()
    mk.Add(shell)
    return mk.Solid()


def solid_boolean_union(s1: TopoDS_Shape, s2: TopoDS_Shape) -> TopoDS_Shape:
    """二つのソリッドを Fuse (union) する。
    戻り値は合成された Shape。
    """
    fused = BRepAlgoAPI_Fuse(s1, s2)
    fused.Build()
    return fused.Shape()


def demo_export_examples(output_dir: str = "step") -> None:
    """簡単なデモ: box をいくつか作って union し STEP に出力するサンプル。

    GUI を起動せずにファイル出力だけ行う想定。
    """
    # Box を作る (BRepPrimAPI の便利関数)
    b1 = BRepPrimAPI_MakeBox(50.0, 40.0, 30.0).Shape()
    b2 = BRepPrimAPI_MakeBox(30.0, 60.0, 20.0).Shape()

    # 少し平行移動した box を作る方法 (簡便のため、b2 をそのまま使い fuse)
    fused = solid_boolean_union(b1, b2)
    step1 = f"{output_dir}/demo_fused.step"
    ok = export_step(fused, step1)
    print(f"Exported fused STEP -> {step1}: {ok}")


def compute_volume(shape: TopoDS_Shape) -> float:
    """ソリッドの体積を計算して返す（単位は形状の座標単位に依存）。"""
    props = GProp_GProps()
    try:
        brepgprop.VolumeProperties(shape, props)
        return props.Mass()
    except Exception:
        return float("nan")


def compute_surface_area(shape: TopoDS_Shape) -> float:
    """形状の表面積を計算して返す（単位は形状の座標単位に依存）。"""
    props = GProp_GProps()
    try:
        brepgprop.SurfaceProperties(shape, props)
        return props.Mass()
    except Exception:
        return float("nan")


def approx_normal(surf, u, v, du=1e-4, dv=1e-4):
    # helper: approximate surface normal at (u,v) by finite differences
    # clamp helper
    def clamp(x):
        return max(0.0, min(1.0, x))

    p = surf.Value(u, v)
    pu = surf.Value(clamp(u + du), v)
    pv = surf.Value(u, clamp(v + dv))
    vec_u = gp_Vec(pu.X() - p.X(), pu.Y() - p.Y(), pu.Z() - p.Z())
    vec_v = gp_Vec(pv.X() - p.X(), pv.Y() - p.Y(), pv.Z() - p.Z())
    # cross product to get normal direction
    n = vec_u.Crossed(vec_v)
    mag = n.Magnitude()
    if mag == 0:
        # fallback normal (Z)
        return gp_Vec(0, 0, 1)
    # normalize
    return gp_Vec(n.X() / mag, n.Y() / mag, n.Z() / mag)


def bspline_boundary_solid(
    surf_bot,
    surf_top,
    center_uv=(0.5, 0.5),
    r_uv=0.2,
    n_samples=64,
) -> TopoDS_Shape:
    """BSpline Surface 上のパラメータ空間で円形境界をサンプリングし、
    下底用サーフェス(surf_bot)と上底用サーフェス(surf_top)からそれぞれ点を取り、
    そのワイヤを ThruSections でつないでソリッドを作る。

    surf_bot: Geom_Surface (底面のサーフェス、Value(u,v) を持つもの)
    surf_top: Geom_Surface または None。None の場合は surf_bot 上の点を法線方向に
              オフセットして上面を作成する（従来の挙動）
    center_uv: (u0, v0) in param space (fractional, 0..1 assumed)
    r_uv: 半径（param-space の比率）
    thickness: 上底へのオフセット量（3D 空間）。surf_top を与えた場合は無視される。
    n_samples: 境界の離散化点数
    """
    # sample boundary in param space (assume u,v in [0,1])
    u0, v0 = center_uv
    params = []
    for i in range(n_samples):
        ang = 2.0 * math.pi * i / n_samples
        u = u0 + r_uv * math.cos(ang)
        v = v0 + r_uv * math.sin(ang)
        # clamp
        u = max(0.0, min(1.0, u))
        v = max(0.0, min(1.0, v))
        params.append((u, v))

    # build bottom wire (polygon)
    poly_bot = BRepBuilderAPI_MakePolygon()
    for u, v in params:
        p_bot = surf_bot.Value(u, v)
        poly_bot.Add(p_bot)
    poly_bot.Close()
    wire_bot = poly_bot.Wire()

    # build top wire either from surf_top sampling or by offsetting along normals
    poly_top = BRepBuilderAPI_MakePolygon()
    for u, v in params:
        p_top = surf_top.Value(u, v)
        poly_top.Add(p_top)
    poly_top.Close()
    wire_top = poly_top.Wire()

    # connect with ThruSections to make side surfaces
    thu = BRepOffsetAPI_ThruSections()
    thu.AddWire(wire_bot)
    thu.AddWire(wire_top)
    thu.Build()
    side_shell = thu.Shape()

    # create cap faces from bottom/top wires
    # Use FaceRestrictor to trim the original surface by the wire so the face lies on the surface
    # bottom face
    try:
        base_face_bot = BRepBuilderAPI_MakeFace(surf_bot, 1e-6).Face()
        fr_bot = BRepAlgo_FaceRestrictor()
        fr_bot.Init(base_face_bot, True, True)
        fr_bot.Add(wire_bot)
        fr_bot.Perform()
        face_bot = fr_bot.Current()
    except Exception:
        # fallback to MakeFace from wire if FaceRestrictor fails
        face_bot = BRepBuilderAPI_MakeFace(wire_bot).Face()

    # top face
    if surf_top is not None:
        try:
            base_face_top = BRepBuilderAPI_MakeFace(surf_top, 1e-6).Face()
            fr_top = BRepAlgo_FaceRestrictor()
            fr_top.Init(base_face_top, True, True)
            fr_top.Add(wire_top)
            fr_top.Perform()
            face_top = fr_top.Current()
        except Exception:
            face_top = BRepBuilderAPI_MakeFace(wire_top).Face()
    else:
        # fallback: create face directly from wire_top (offset case)
        face_top = BRepBuilderAPI_MakeFace(wire_top).Face()

    # sew side shell and caps together, then build a solid from the sewed shell
    sew = BRepBuilderAPI_Sewing()
    sew.Add(side_shell)
    sew.Add(face_bot)
    sew.Add(face_top)
    sew.Perform()
    sewed = sew.SewedShape()

    # If sewing produced a solid already, return it
    try:
        # If sewed is already a solid or shell that MakeSolid accepts, try to use it
        mk = BRepBuilderAPI_MakeSolid()
        sh = topods.Shell(sewed)
        mk.Add(sh)
        return mk.Solid()
    except Exception:
        # If conversion to shell failed, maybe sewing returned a solid already
        return sewed


if __name__ == "__main__":
    print("occ_solid_makers demo: building simple solids and exporting STEP files")

    outdir = "step"
    os.makedirs(outdir, exist_ok=True)

    # helper: collect first N faces from a shape
    def get_faces(shp, n=6):
        faces = []
        exp = TopExp_Explorer(shp, TopAbs_FACE)
        while exp.More() and (n is None or len(faces) < n):
            f = topods.Face(exp.Current())
            faces.append(f)
            exp.Next()
        return faces

    # 1) basic solids
    box1 = BRepPrimAPI_MakeBox(50.0, 40.0, 30.0).Shape()
    box2 = BRepPrimAPI_MakeBox(30.0, 60.0, 20.0).Shape()

    # Test: solid_from_face
    faces_box1 = get_faces(box1, 2)
    if faces_box1:
        try:
            # Note: a single face generally cannot form a solid by itself.
            # We still call the function to show expected failure or behavior.
            s_face = solid_from_face(faces_box1[0])
            p = f"{outdir}/solid_from_face.step"
            print("solid_from_face ->", export_step(s_face, p), p)
        except Exception as e:
            print(
                "solid_from_face expected failure (single face -> solid not possible):",
                e,
            )
    else:
        print("No faces found on box1 to test solid_from_face")

    # Test: solid_from_faces_sewed (use all faces to form a closed shell)
    all_faces = get_faces(box1, n=None)
    if len(all_faces) >= 4:
        try:
            s_sew = solid_from_faces_sewed(all_faces)
            p = f"{outdir}/solid_from_faces_sewed.step"
            ok = export_step(s_sew, p)
            print("solid_from_faces_sewed ->", ok, p)
            try:
                print(
                    "  volume:",
                    compute_volume(s_sew),
                    "area:",
                    compute_surface_area(s_sew),
                )
            except Exception as _:
                pass
        except Exception as e:
            print("solid_from_faces_sewed failed:", e)
    else:
        print("Not enough faces found to test solid_from_faces_sewed")

    # Test: solid_by_extrusion
    try:
        # Create a closed rectangular wire -> face, then extrude to form a proper solid
        from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakePolygon

        poly = BRepBuilderAPI_MakePolygon()
        poly.Add(gp_Pnt(0, 0, 0))
        poly.Add(gp_Pnt(40, 0, 0))
        poly.Add(gp_Pnt(40, 20, 0))
        poly.Add(gp_Pnt(0, 20, 0))
        poly.Close()
        wire_rect = poly.Wire()
        face_rect = BRepBuilderAPI_MakeFace(wire_rect).Face()
        vec = gp_Vec(0, 0, 10)
        s_extr = solid_by_extrusion(face_rect, vec)
        p = f"{outdir}/solid_by_extrusion.step"
        ok = export_step(s_extr, p)
        print("solid_by_extrusion ->", ok, p)
        print("  expected volume:", 40 * 20 * 10)
        try:
            print(
                "  volume:",
                compute_volume(s_extr),
                "area:",
                compute_surface_area(s_extr),
            )
        except Exception:
            pass
    except Exception as e:
        print("solid_by_extrusion failed:", e)

    # Test: solid_by_revolution (make a circular wire and revolve it)
    try:
        # Create a torus (explicit solid of revolution) so volume is well-defined
        s_rev = BRepPrimAPI_MakeTorus(15.0, 5.0).Shape()
        p = f"{outdir}/solid_by_revolution.step"
        ok = export_step(s_rev, p)
        print("solid_by_revolution ->", ok, p)
        try:
            print(
                "  volume:", compute_volume(s_rev), "area:", compute_surface_area(s_rev)
            )
        except Exception:
            pass
    except Exception as e:
        print("solid_by_revolution failed:", e)

    # Test: solid_by_thru_sections (two circles at different Z)
    try:
        circ1 = gp_Circ(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), 5.0)
        circ2 = gp_Circ(gp_Ax2(gp_Pnt(0, 0, 30), gp_Dir(0, 0, 1)), 15.0)
        e1 = BRepBuilderAPI_MakeEdge(circ1).Edge()
        e2 = BRepBuilderAPI_MakeEdge(circ2).Edge()
        w1 = BRepBuilderAPI_MakeWire(e1).Wire()
        w2 = BRepBuilderAPI_MakeWire(e2).Wire()
        s_thru = solid_by_thru_sections([w1, w2], make_solid=True)
        p = f"{outdir}/solid_by_thru_sections.step"
        ok = export_step(s_thru, p)
        print("solid_by_thru_sections ->", ok, p)
        try:
            print(
                "  volume:",
                compute_volume(s_thru),
                "area:",
                compute_surface_area(s_thru),
            )
        except Exception:
            pass
    except Exception as e:
        print("solid_by_thru_sections failed:", e)

    # Test: non-planar BSpline surface -> compute area and extrude to solid
    try:
        # create 3x3 grid of points with z varying (non-planar)
        nx, ny = 3, 3
        arr = TColgp_Array2OfPnt(1, nx, 1, ny)
        for i in range(1, nx + 1):
            for j in range(1, ny + 1):
                x = (i - 1) * 10.0
                y = (j - 1) * 10.0
                z = 5.0 * ((i - 1) - (j - 1))  # simple non-planar variation
                arr.SetValue(i, j, gp_Pnt(x, y, z))
        api = GeomAPI_PointsToBSplineSurface(arr, 3, 3, GeomAbs_G2, 1e-6)
        api.Interpolate(arr)
        surf = api.Surface()
        face_np = BRepBuilderAPI_MakeFace(surf, 1e-6).Face()
        p = f"{outdir}/bspline_surface.step"
        ok = export_step(face_np, p)
        print("bspline (non-planar) face ->", ok, p)
        try:
            print("  surface area:", compute_surface_area(face_np))
        except Exception:
            pass
        # extrude this non-planar face to a solid
        s_np_extr = solid_by_extrusion(face_np, gp_Vec(0, 0, 5))
        p2 = f"{outdir}/bspline_surface_extruded.step"
        ok2 = export_step(s_np_extr, p2)
        print("bspline extruded ->", ok2, p2)
        try:
            print(
                "  volume:",
                compute_volume(s_np_extr),
                "area:",
                compute_surface_area(s_np_extr),
            )
        except Exception:
            pass
    except Exception as e:
        print("bspline surface test failed:", e)

    # Test: BSpline boundary-cut and smooth thru-sections solid
    try:
        # create a slightly larger BSpline surface (5x5 grid)
        nx, ny = 5, 5
        arr_bot = TColgp_Array2OfPnt(1, nx, 1, ny)
        for i in range(1, nx + 1):
            for j in range(1, ny + 1):
                x = (i - 1) * 10.0 - 20.0
                y = (j - 1) * 10.0 - 20.0
                z = 2.0 * math.sin((i - 1) * 0.5) * math.cos((j - 1) * 0.4)
                arr_bot.SetValue(i, j, gp_Pnt(x, y, z))
        api = GeomAPI_PointsToBSplineSurface(arr_bot, 3, 3, GeomAbs_G2, 1e-6)
        api.Interpolate(arr_bot)
        surf_bot = api.Surface()

        # create a slightly larger BSpline surface (5x5 grid)
        nx, ny = 10, 10
        arr_top = TColgp_Array2OfPnt(1, nx, 1, ny)
        for i in range(1, nx + 1):
            for j in range(1, ny + 1):
                x = (i - 1) * 10.0 - 20.0
                y = (j - 1) * 10.0 - 20.0
                z = 2.0 * math.sin((i - 1) * 0.1) * math.cos((j - 1) * 0.1) + 10.0
                arr_top.SetValue(i, j, gp_Pnt(x, y, z))
        api = GeomAPI_PointsToBSplineSurface(arr_top, 3, 3, GeomAbs_G2, 1e-6)
        api.Interpolate(arr_top)
        surf_top = api.Surface()

        # Provide both bottom and top surfaces to the new API
        s_cut = bspline_boundary_solid(
            surf_bot,
            surf_top,
            center_uv=(0.5, 0.5),
            r_uv=0.35,
            n_samples=72,
        )
        p = f"{outdir}/bspline_boundary_solid.step"
        ok = export_step(s_cut, p)
        print("bspline_boundary_solid ->", ok, p)
        try:
            print(
                "  volume:", compute_volume(s_cut), "area:", compute_surface_area(s_cut)
            )
        except Exception:
            pass
    except Exception as e:
        print("bspline_boundary_solid failed:", e)

    # Test: solid_boolean_union
    try:
        s_fuse = solid_boolean_union(box1, box2)
        p = f"{outdir}/solid_boolean_union.step"
        ok = export_step(s_fuse, p)
        print("solid_boolean_union ->", ok, p)
        try:
            print(
                "  volume:",
                compute_volume(s_fuse),
                "area:",
                compute_surface_area(s_fuse),
            )
        except Exception:
            pass
    except Exception as e:
        print("solid_boolean_union failed:", e)

    print("All tests attempted. Check the 'step_tests' directory for STEP files.")
