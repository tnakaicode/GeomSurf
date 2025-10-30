"""Helpers: several patterns to create solids in pythonOCC

このモジュールは、OCC (pythonOCC) を使った Solid の作り方を関数別に示します。
各関数は TopoDS_Shape (solid) を返します。

注意: 実行環境に pythonOCC が必要です。表示系を使わず STEP export による出力サンプルを付けています。
"""

from __future__ import annotations
import os

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
