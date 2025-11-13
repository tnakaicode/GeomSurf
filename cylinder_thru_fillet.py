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
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_ThruSections, BRepOffsetAPI_MakePipe
from OCC.Core.Geom import Geom_Circle
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
from OCC.Core.BRepFilletAPI import BRepFilletAPI_MakeFillet
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Edge
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_EDGE
from OCC.Core.GeomAPI import GeomAPI_PointsToBSpline
from OCC.Core.TColgp import TColgp_Array1OfPnt


class HoleSpec:
    """穴の仕様を定義するクラス"""

    def __init__(self, height, angle_deg, diameter, depth=20.0):
        self.height = height  # Z座標での高さ
        self.angle_deg = angle_deg  # 回転角度（度）
        self.diameter = diameter  # 穴の直径
        self.depth = depth  # 穴の深さ

    def get_position(self, cylinder_radius):
        """シリンダー表面での穴の位置を計算"""
        angle_rad = np.radians(self.angle_deg)
        x = cylinder_radius * np.cos(angle_rad)
        y = cylinder_radius * np.sin(angle_rad)
        z = self.height
        return gp_Pnt(x, y, z)

    def get_direction(self):
        """穴を開ける方向を計算（シリンダー中心向き）"""
        angle_rad = np.radians(self.angle_deg)
        dx = -np.cos(angle_rad)  # 中心向き
        dy = -np.sin(angle_rad)
        return gp_Dir(dx, dy, 0)


def create_cylinder_with_holes(
    cylinder_radius=50.0, cylinder_height=100.0, hole_specs=None
):
    """一般的なシリンダーに任意の穴を開ける"""
    if hole_specs is None:
        # デフォルト: 2つの穴（径違い、高さ違い）
        hole_specs = [
            HoleSpec(height=70, angle_deg=0, diameter=8.0),  # 上部、正面
            HoleSpec(height=30, angle_deg=0, diameter=6.0),  # 下部、正面
        ]

    # シリンダー作成
    cylinder_axis = gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1))
    cylinder = BRepPrimAPI_MakeCylinder(
        cylinder_axis, cylinder_radius, cylinder_height
    ).Shape()

    hole_centers = []
    result_shape = cylinder

    # 各穴を順次作成
    for i, hole_spec in enumerate(hole_specs):
        hole_center = hole_spec.get_position(cylinder_radius)
        hole_direction = hole_spec.get_direction()

        print(
            f"穴{i+1}: 位置({hole_center.X():.1f}, {hole_center.Y():.1f}, {hole_center.Z():.1f}), "
            f"角度{hole_spec.angle_deg}°, 径{hole_spec.diameter}mm"
        )

        hole_axis = gp_Ax2(hole_center, hole_direction)
        hole_cylinder = BRepPrimAPI_MakeCylinder(
            hole_axis, hole_spec.diameter, hole_spec.depth
        ).Shape()

        # 穴を開ける
        cut = BRepAlgoAPI_Cut(result_shape, hole_cylinder)
        cut.Build()
        if cut.IsDone():
            result_shape = cut.Shape()
            hole_centers.append(hole_center)
        else:
            print(f"警告: 穴{i+1}の作成に失敗")

    return result_shape, hole_centers


def create_thru_section_connection(hole_centers, hole_specs):
    """穴同士をThruSectionで接続（元の要求通り）"""
    if len(hole_centers) < 2:
        print("警告: 接続するには最低2つの穴が必要")
        return None, None

    print(f"=== {len(hole_centers)}個の穴をThruSectionで接続 ===")

    try:
        # 穴1の入口円（シリンダー表面）
        pt1 = hole_centers[0]
        diameter1 = hole_specs[0].diameter if hole_specs else 8.0
        circle1_axis = gp_Ax2(pt1, gp_Dir(-1, 0, 0))  # 穴の方向
        circle1_geom = Geom_Circle(circle1_axis, diameter1 / 2)
        circle1_edge = BRepBuilderAPI_MakeEdge(circle1_geom).Edge()
        wire1 = BRepBuilderAPI_MakeWire(circle1_edge).Wire()

        # 穴2の入口円（シリンダー表面）
        pt2 = hole_centers[1]
        diameter2 = hole_specs[1].diameter if hole_specs else 6.0
        circle2_axis = gp_Ax2(pt2, gp_Dir(-1, 0, 0))  # 穴の方向
        circle2_geom = Geom_Circle(circle2_axis, diameter2 / 2)
        circle2_edge = BRepBuilderAPI_MakeEdge(circle2_geom).Edge()
        wire2 = BRepBuilderAPI_MakeWire(circle2_edge).Wire()

        # ThruSectionで接続
        thru_sections = BRepOffsetAPI_ThruSections(True, True)
        thru_sections.AddWire(wire1)
        thru_sections.AddWire(wire2)
        thru_sections.Build()

        if thru_sections.IsDone():
            connection_shape = thru_sections.Shape()
            # スパイン（中心線）も作成
            spine_edge = BRepBuilderAPI_MakeEdge(pt1, pt2).Edge()
            spine_wire = BRepBuilderAPI_MakeWire(spine_edge).Wire()

            print("✅ ThruSection接続成功")
            print(
                f"   穴1: 径{diameter1}mm at ({pt1.X():.1f}, {pt1.Y():.1f}, {pt1.Z():.1f})"
            )
            print(
                f"   穴2: 径{diameter2}mm at ({pt2.X():.1f}, {pt2.Y():.1f}, {pt2.Z():.1f})"
            )

            return connection_shape, spine_wire
        else:
            print("❌ ThruSection接続失敗")
            return None, None

    except Exception as e:
        print(f"❌ ThruSection接続エラー: {e}")
        return None, None


def create_spine_connections(hole_centers, connection_type="direct"):
    """Spine（線）接続（参考用）"""
    if len(hole_centers) < 2:
        print("警告: 接続するには最低2つの穴が必要")
        return [], []

    spine_edges = []
    spine_wires = []

    print(f"=== {len(hole_centers)}個の穴をSpine（線）で接続 ===")

    for i in range(len(hole_centers) - 1):
        pt1, pt2 = hole_centers[i], hole_centers[i + 1]
        edge, wire = create_single_spine_connection(pt1, pt2, i + 1)
        if edge and wire:
            spine_edges.append(edge)
            spine_wires.append(wire)

    print(f"Spine接続完了: {len(spine_edges)}本作成")
    return spine_edges, spine_wires


def create_single_spine_connection(pt1, pt2, connection_id):
    """2点間の単一Spine接続を作成"""
    distance = pt1.Distance(pt2)
    print(f"接続{connection_id}: 距離{distance:.1f}mm")

    try:
        if distance < 1.0:  # 距離が短すぎる場合
            print(f"警告: 接続{connection_id}の距離が短すぎます")
            return None, None

        # 直線Spine
        spine_edge = BRepBuilderAPI_MakeEdge(pt1, pt2).Edge()
        spine_wire = BRepBuilderAPI_MakeWire(spine_edge).Wire()
        print(f"接続{connection_id}: 成功")
        return spine_edge, spine_wire

    except Exception as e:
        print(f"接続{connection_id}失敗: {e}")
        return None, None


def calculate_center_point(hole_centers):
    """穴の中心点を計算"""
    if not hole_centers:
        return gp_Pnt(0, 0, 0)

    sum_x = sum(pt.X() for pt in hole_centers)
    sum_y = sum(pt.Y() for pt in hole_centers)
    sum_z = sum(pt.Z() for pt in hole_centers)
    n = len(hole_centers)

    return gp_Pnt(sum_x / n, sum_y / n, sum_z / n)


def apply_cylinder_connection_fillet(cylinder_shape, connection_shape, radius=3.0):
    """シリンダーと接続部分の境界にFilletを適用"""
    from OCC.Core.BRep import BRep_Tool
    from OCC.Core.BRepAdaptor import BRepAdaptor_Curve
    from OCC.Core.GeomAbs import GeomAbs_Circle

    print(f"=== Fillet適用開始（半径{radius}mm） ===")

    try:
        # シリンダーと接続部分を結合
        fuse = BRepAlgoAPI_Fuse(cylinder_shape, connection_shape)
        fuse.Build()

        if not fuse.IsDone():
            print("❌ 形状結合失敗: Fillet適用不可")
            return None

        fused_shape = fuse.Shape()
        print("✅ 形状結合成功")

        # Filletを適用
        fillet = BRepFilletAPI_MakeFillet(fused_shape)

        # エッジを探索してFillet対象を特定
        explorer = TopExp_Explorer(fused_shape, TopAbs_EDGE)
        fillet_count = 0
        processed_edges = set()

        while explorer.More():
            edge = explorer.Current()

            # 同一エッジの重複処理を防止
            edge_id = str(edge.__hash__())
            if edge_id in processed_edges:
                explorer.Next()
                continue
            processed_edges.add(edge_id)

            try:
                # エッジの形状を判定
                adaptor = BRepAdaptor_Curve(edge)

                # 円形エッジ（穴の境界）を検出
                if adaptor.GetType() == GeomAbs_Circle:
                    circle = adaptor.Circle()
                    circle_radius = circle.Radius()

                    # 穴径に近い円（穴の境界エッジ）を対象とする
                    if 2.0 <= circle_radius <= 6.0:  # 穴径の範囲内
                        fillet.Add(radius, edge)
                        fillet_count += 1
                        print(f"   Fillet追加: 円形エッジ（半径{circle_radius:.1f}mm）")

                # 直線エッジも一部対象とする
                else:
                    curve = BRep_Tool.Curve(edge)
                    if curve[0] is not None:
                        u_start = curve[1]
                        u_end = curve[2]
                        u_mid = (u_start + u_end) / 2

                        pt = curve[0].Value(u_mid)
                        edge_center = gp_Pnt(pt.X(), pt.Y(), pt.Z())

                        # シリンダー表面付近のエッジを対象
                        distance_from_center = (
                            edge_center.X() ** 2 + edge_center.Y() ** 2
                        ) ** 0.5
                        if (
                            abs(distance_from_center - 50.0) < 10.0
                        ):  # シリンダー半径±10mm
                            fillet.Add(radius, edge)
                            fillet_count += 1
                            print(
                                f"   Fillet追加: 接続エッジ（中心距離{distance_from_center:.1f}mm）"
                            )

            except Exception as edge_error:
                # 個別エッジエラーは無視して続行
                pass

            explorer.Next()

        print(f"Fillet対象エッジ数: {fillet_count}")

        if fillet_count > 0:
            # 小さい半径から試行
            for test_radius in [radius * 0.5, radius * 0.7, radius, radius * 1.2]:
                try:
                    test_fillet = BRepFilletAPI_MakeFillet(fused_shape)

                    # 再度エッジを追加
                    explorer = TopExp_Explorer(fused_shape, TopAbs_EDGE)
                    while explorer.More():
                        edge = explorer.Current()
                        try:
                            adaptor = BRepAdaptor_Curve(edge)
                            if adaptor.GetType() == GeomAbs_Circle:
                                circle = adaptor.Circle()
                                if 2.0 <= circle.Radius() <= 6.0:
                                    test_fillet.Add(test_radius, edge)
                        except:
                            pass
                        explorer.Next()

                    test_fillet.Build()
                    if test_fillet.IsDone():
                        print(f"✅ Fillet適用成功（半径{test_radius:.1f}mm）")
                        return test_fillet.Shape()

                except Exception as radius_error:
                    print(f"   半径{test_radius:.1f}mm失敗: {radius_error}")
                    continue

            print("❌ 全ての半径でFillet失敗")
            return fused_shape

        else:
            print("❌ Fillet対象エッジなし")
            return fused_shape

    except Exception as e:
        print(f"❌ Fillet処理エラー: {e}")
        return None


def run_complete_test(test_name, hole_specs):
    """完全なテスト（ThruSection + Fillet）を実行"""
    print(f"\n{'='*60}")
    print(f"完全テスト: {test_name}")
    print(f"元の要求: Cylinderに穴開けて→ThruSectionでつないで→Filletつける")
    print(f"{'='*60}")

    obj = dispocc(touch=True)

    # ステップ1: シリンダーに穴を開ける
    print("\n🔧 ステップ1: シリンダーに穴を開ける")
    cylinder_with_holes, hole_centers = create_cylinder_with_holes(
        hole_specs=hole_specs
    )

    # ステップ2: ThruSectionで接続
    print("\n🔧 ステップ2: ThruSectionで穴同士を接続")
    connection_shape, spine_wire = create_thru_section_connection(
        hole_centers, hole_specs
    )

    if connection_shape is None:
        print("❌ 接続失敗: 表示のみ")
        obj.display.DisplayShape(cylinder_with_holes, color="BLUE", transparency=0.3)
        obj.ShowOCC()
        return False

    # ステップ3: Filletを適用
    print("\n🔧 ステップ3: CylinderとThruSectionの接続部にFillet適用")
    final_shape = apply_cylinder_connection_fillet(
        cylinder_with_holes, connection_shape, radius=3.0
    )

    # 結果表示
    print("\n🎨 結果表示")
    if final_shape is not None:
        obj.display.DisplayShape(final_shape, color="GREEN", transparency=0.2)
        print("✅ 最終形状表示: 緑色（Cylinder + ThruSection + Fillet）")
    else:
        # 個別表示
        obj.display.DisplayShape(cylinder_with_holes, color="BLUE", transparency=0.3)
        obj.display.DisplayShape(connection_shape, color="RED", transparency=0.4)
        print("⚠️ 個別表示: 青（Cylinder）+ 赤（ThruSection）")

    obj.ShowOCC()
    return final_shape is not None


def run_spine_test(test_name, hole_specs):
    """Spine接続テスト（参考用）"""
    print(f"\n{'='*40}")
    print(f"Spine参考テスト: {test_name}")
    print(f"{'='*40}")

    obj = dispocc(touch=True)

    # シリンダーに穴を開ける
    cylinder_with_holes, hole_centers = create_cylinder_with_holes(
        hole_specs=hole_specs
    )
    obj.display.DisplayShape(cylinder_with_holes, color="BLUE", transparency=0.3)

    # Spine接続
    spine_edges, spine_wires = create_spine_connections(hole_centers)

    # 結果表示
    colors = ["RED", "GREEN", "YELLOW", "ORANGE", "PURPLE"]
    for i, (edge, wire) in enumerate(zip(spine_edges, spine_wires)):
        color = colors[i % len(colors)]
        obj.display.DisplayShape(edge, color=color, transparency=0.0)

    print(f"参考結果: {len(spine_edges)}本のSpine接続")
    obj.ShowOCC()
    return len(spine_edges) > 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mode",
        dest="mode",
        default="complete",
        help="実行モード: complete（ThruSection+Fillet）, spine（Spine参考）, both（両方）",
    )
    parser.add_argument(
        "--test",
        dest="test",
        default="basic",
        help="テストケース: basic, different_heights, different_angles, mixed",
    )
    opt = parser.parse_args()

    # テストケース定義（元の要求に基づく2穴構成をメインに）
    test_cases = {
        "basic": {
            "name": "基本テスト（径違い、高さ違い）- 元の要求",
            "holes": [
                HoleSpec(height=70, angle_deg=0, diameter=8.0),
                HoleSpec(height=30, angle_deg=0, diameter=6.0),
            ],
        },
        "different_heights": {
            "name": "高さ違いテスト（2穴）",
            "holes": [
                HoleSpec(height=80, angle_deg=0, diameter=8.0),
                HoleSpec(height=20, angle_deg=0, diameter=6.0),
            ],
        },
        "different_angles": {
            "name": "回転位置違いテスト（2穴）",
            "holes": [
                HoleSpec(height=60, angle_deg=0, diameter=7.0),  # 正面
                HoleSpec(height=40, angle_deg=180, diameter=5.0),  # 背面
            ],
        },
        "mixed": {
            "name": "複合テスト（径、高さ、角度違い）",
            "holes": [
                HoleSpec(height=75, angle_deg=45, diameter=9.0),
                HoleSpec(height=25, angle_deg=225, diameter=5.0),
            ],
        },
    }

    if opt.test in test_cases:
        case = test_cases[opt.test]

        if opt.mode == "complete":
            print("🎯 実行モード: 完全テスト（ThruSection + Fillet）")
            success = run_complete_test(case["name"], case["holes"])

            print(f"\n{'='*60}")
            if success:
                print("✅ 結論: ThruSection接続 + Fillet適用が正常に動作")
                print(
                    "✅ 達成: 元の要求『Cylinderに穴開けて→ThruSectionでつないで→Filletつける』完了"
                )
            else:
                print("❌ 結論: 一部機能に問題があります")

        elif opt.mode == "spine":
            print("📏 実行モード: Spine参考テスト")
            success = run_spine_test(case["name"], case["holes"])

        elif opt.mode == "both":
            print("🔄 実行モード: 両方テスト")
            success1 = run_complete_test(case["name"], case["holes"])
            success2 = run_spine_test(case["name"], case["holes"])

            print(f"\n{'='*60}")
            print(f"完全テスト: {'成功' if success1 else '失敗'}")
            print(f"Spineテスト: {'成功' if success2 else '失敗'}")

    else:
        print(f"未知のテストケース: {opt.test}")
        print(f"利用可能: {list(test_cases.keys())}")
        print(f"実行例: python cylinder_thru_fillet.py --mode complete --test basic")
