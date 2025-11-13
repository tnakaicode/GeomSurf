import numpy as np
import sys
import os
import argparse

basename = os.path.dirname(__file__)
sys.path.append(os.path.join("./"))
from src.base_occ import dispocc

import logging

logging.getLogger("matplotlib").setLevel(logging.ERROR)

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir, gp_Ax1, gp_Ax2, gp_Ax3, gp_Circ
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeCylinder
from OCC.Core.BRepBuilderAPI import (
    BRepBuilderAPI_MakeWire,
    BRepBuilderAPI_MakeFace,
    BRepBuilderAPI_MakeEdge,
)
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_ThruSections
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
from OCC.Core.BRepFilletAPI import BRepFilletAPI_MakeFillet
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Edge
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_EDGE
from OCC.Core.Geom import Geom_Circle
from OCC.Core.BRepAdaptor import BRepAdaptor_Curve
from OCC.Core.GeomAbs import GeomAbs_Circle
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCCUtils.Construct import make_wire, make_face


def create_cylinder_with_holes():
    """基本のシリンダーを作成し、複数の穴を開ける"""
    # メインシリンダーを作成
    cylinder_axis = gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1))
    main_cylinder = BRepPrimAPI_MakeCylinder(cylinder_axis, 50.0, 100.0).Shape()

    # 穴1: 上部の円形穴
    hole1_axis = gp_Ax2(gp_Pnt(20, 0, 80), gp_Dir(0, 0, 1))
    hole1_cylinder = BRepPrimAPI_MakeCylinder(hole1_axis, 10.0, 20.0).Shape()

    # 穴2: 下部の円形穴
    hole2_axis = gp_Ax2(gp_Pnt(-20, 0, 20), gp_Dir(0, 0, 1))
    hole2_cylinder = BRepPrimAPI_MakeCylinder(hole2_axis, 8.0, 20.0).Shape()

    # 穴3: 側面の穴1（X軸方向）- 位置を記録
    hole3_center = gp_Pnt(0, 0, 60)
    hole3_axis = gp_Ax2(hole3_center, gp_Dir(1, 0, 0))
    hole3_cylinder = BRepPrimAPI_MakeCylinder(hole3_axis, 6.0, 60.0).Shape()

    # 穴4: 側面の穴2（Y軸方向）- 位置を記録
    hole4_center = gp_Pnt(0, 0, 40)
    hole4_axis = gp_Ax2(hole4_center, gp_Dir(0, 1, 0))
    hole4_cylinder = BRepPrimAPI_MakeCylinder(hole4_axis, 5.0, 60.0).Shape()

    # 順次穴を開ける
    cut1 = BRepAlgoAPI_Cut(main_cylinder, hole1_cylinder)
    cut1.Build()
    cylinder_with_hole1 = cut1.Shape()

    cut2 = BRepAlgoAPI_Cut(cylinder_with_hole1, hole2_cylinder)
    cut2.Build()
    cylinder_with_hole2 = cut2.Shape()

    cut3 = BRepAlgoAPI_Cut(cylinder_with_hole2, hole3_cylinder)
    cut3.Build()
    cylinder_with_hole3 = cut3.Shape()

    cut4 = BRepAlgoAPI_Cut(cylinder_with_hole3, hole4_cylinder)
    cut4.Build()
    cylinder_with_holes = cut4.Shape()

    # 上下穴の中心と側面穴の情報を返す
    side_holes_info = [
        {"center": hole3_center, "radius": 6.0, "direction": gp_Dir(1, 0, 0)},
        {"center": hole4_center, "radius": 5.0, "direction": gp_Dir(0, 1, 0)},
    ]

    return cylinder_with_holes, gp_Pnt(20, 0, 80), gp_Pnt(-20, 0, 20), side_holes_info


def create_thru_section(pt1, pt2):
    """2つの点を結ぶThruSectionを作成"""
    # 上部の円
    circle1_axis = gp_Ax2(pt1, gp_Dir(0, 0, 1))
    circle1_geom = Geom_Circle(circle1_axis, 10.0)
    circle1_edge = BRepBuilderAPI_MakeEdge(circle1_geom).Edge()
    wire1 = BRepBuilderAPI_MakeWire(circle1_edge).Wire()

    # 下部の円
    circle2_axis = gp_Ax2(pt2, gp_Dir(0, 0, 1))
    circle2_geom = Geom_Circle(circle2_axis, 8.0)
    circle2_edge = BRepBuilderAPI_MakeEdge(circle2_geom).Edge()
    wire2 = BRepBuilderAPI_MakeWire(circle2_edge).Wire()

    # ThruSectionでつなぐ
    thru_sections = BRepOffsetAPI_ThruSections(True, True)
    thru_sections.AddWire(wire1)
    thru_sections.AddWire(wire2)
    thru_sections.Build()

    return thru_sections.Shape()


def apply_side_hole_fillet(cylinder_shape, side_holes_info, radius=3.0):
    """側面穴の周囲のエッジにFilletを適用"""
    from OCC.Core.BRep import BRep_Tool
    from OCC.Core.TopoDS import topods

    fillet = BRepFilletAPI_MakeFillet(cylinder_shape)

    # エッジを探して、側面穴の近くにあるエッジにFilletを適用
    explorer = TopExp_Explorer(cylinder_shape, TopAbs_EDGE)
    fillet_count = 0

    while explorer.More():
        edge = explorer.Current()

        # エッジの中点を取得
        curve = BRep_Tool.Curve(edge)
        if curve[0] is not None:
            u_start = curve[1]
            u_end = curve[2]
            u_mid = (u_start + u_end) / 2

            # エッジの中点座標を取得
            pt = curve[0].Value(u_mid)
            edge_center = gp_Pnt(pt.X(), pt.Y(), pt.Z())

            # 各側面穴の中心からの距離をチェック
            for hole_info in side_holes_info:
                hole_center = hole_info["center"]
                hole_radius = hole_info["radius"]

                # 穴の中心からエッジまでの距離
                distance = edge_center.Distance(hole_center)

                # 穴の半径の近くにあるエッジ（穴の境界部分）をFillet対象とする
                if abs(distance - hole_radius) < 8.0:  # 許容範囲内
                    fillet.Add(radius, edge)
                    fillet_count += 1
                    print(
                        f"Fillet追加: 穴中心({hole_center.X():.1f}, {hole_center.Y():.1f}, {hole_center.Z():.1f}), 距離: {distance:.2f}"
                    )
                    break  # 一つの穴に対してマッチしたら次のエッジへ

        explorer.Next()

    print(f"Fillet対象エッジ数: {fillet_count}")

    if fillet_count > 0:
        fillet.Build()
        if fillet.IsDone():
            print("Fillet適用成功")
            return fillet.Shape()
        else:
            print(f"Fillet作成失敗: 半径{radius}が大きすぎる可能性があります")
            return cylinder_shape
    else:
        print("Fillet対象エッジが見つかりませんでした")
        return cylinder_shape


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--radius", dest="radius", default=3.0, type=float)
    opt = parser.parse_args()

    obj = dispocc(touch=True)

    # ステップ1: シリンダーに穴を開ける
    cylinder_with_holes, hole1_center, hole2_center, side_holes_info = (
        create_cylinder_with_holes()
    )

    # ステップ2: 側面穴の部分だけにFilletを適用
    cylinder_filleted = apply_side_hole_fillet(
        cylinder_with_holes, side_holes_info, opt.radius
    )
    obj.display.DisplayShape(cylinder_filleted, color="BLUE", transparency=0.3)

    # ステップ3: ThruSectionを作成
    thru_section = create_thru_section(hole1_center, hole2_center)
    obj.display.DisplayShape(thru_section, color="RED", transparency=0.3)

    # ステップ4: 形状を結合（Filletは側面穴のみに適用済み）
    fuse = BRepAlgoAPI_Fuse(cylinder_filleted, thru_section)
    fuse.Build()

    if fuse.IsDone():
        final_shape = fuse.Shape()
        obj.display.DisplayShape(final_shape, color="GREEN", transparency=0.8)
        print("結合成功: Filletが保持されました")
    else:
        print("結合失敗: 個別に表示します")
        # 結合が失敗した場合は個別に表示
        obj.display.DisplayShape(cylinder_filleted, color="BLUE", transparency=0.3)
        obj.display.DisplayShape(thru_section, color="RED", transparency=0.3)

    print(f"シリンダーに4つの穴（上下2つ、側面2つ）を開け、ThruSectionで接続し、")
    print(f"側面穴の部分のみにFillet半径{opt.radius}を適用しました。")

    obj.ShowOCC()
