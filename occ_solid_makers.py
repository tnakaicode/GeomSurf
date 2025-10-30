"""Helpers: several patterns to create solids in pythonOCC

このモジュールは、OCC (pythonOCC) を使った Solid の作り方を関数別に示します。
各関数は TopoDS_Shape (solid) を返します。

注意: 実行環境に pythonOCC が必要です。表示系を使わず STEP export による出力サンプルを付けています。
"""

from __future__ import annotations

from OCC.Core.gp import gp_Vec, gp_Ax1, gp_Pnt, gp_Dir
from OCC.Core.TopoDS import TopoDS_Shape
from OCC.Core.BRepBuilderAPI import (
    BRepBuilderAPI_MakeFace,
    BRepBuilderAPI_MakeShell,
    BRepBuilderAPI_MakeSolid,
)
from OCC.Core.BRepPrimAPI import (
    BRepPrimAPI_MakePrism,
    BRepPrimAPI_MakeRevol,
    BRepPrimAPI_MakeBox,
)
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_ThruSections
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_Sewing
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Fuse
from OCC.Core.STEPControl import STEPControl_Writer, STEPControl_AsIs
from OCC.Core.TopLoc import TopLoc_Location
from typing import Iterable, List


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
    shell = BRepBuilderAPI_MakeShell(face).Shell()
    shell.Location(TopLoc_Location())
    mk = BRepBuilderAPI_MakeSolid()
    try:
        mk.Add(shell)
    except Exception:
        # MakeSolid に失敗した場合、呼び出し側で処理してもらう
        raise
    return mk.Solid()


def solid_from_faces_sewed(faces: Iterable) -> TopoDS_Shape:
    """複数の Face を Sewing で縫い合わせて Shell を作り、Solid にする。

    faces: iterable of TopoDS_Face
    """
    sew = BRepBuilderAPI_Sewing()
    for f in faces:
        sew.Add(f)
    sew.Perform()
    res_shell = sew.SewedShape()
    mk = BRepBuilderAPI_MakeSolid()
    # SewedShape は Shell になっているはず
    mk.Add(res_shell)
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


if __name__ == "__main__":
    print("occ_solid_makers demo: building simple solids and exporting STEP files")
    demo_export_examples()
