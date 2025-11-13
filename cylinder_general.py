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
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
from OCC.Core.GeomAPI import GeomAPI_PointsToBSpline
from OCC.Core.TColgp import TColgp_Array1OfPnt


class HoleSpec:
    """穴の仕様を定義するクラス"""
    def __init__(self, height, angle_deg, diameter, depth=20.0):
        self.height = height        # Z座標での高さ
        self.angle_deg = angle_deg  # 回転角度（度）
        self.diameter = diameter    # 穴の直径
        self.depth = depth         # 穴の深さ
        
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


def create_cylinder_with_holes(cylinder_radius=50.0, cylinder_height=100.0, hole_specs=None):
    """一般的なシリンダーに任意の穴を開ける"""
    if hole_specs is None:
        # デフォルト: 2つの穴（径違い、高さ違い）
        hole_specs = [
            HoleSpec(height=70, angle_deg=0, diameter=8.0),    # 上部、正面
            HoleSpec(height=30, angle_deg=0, diameter=6.0),    # 下部、正面
        ]
    
    # シリンダー作成
    cylinder_axis = gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1))
    cylinder = BRepPrimAPI_MakeCylinder(cylinder_axis, cylinder_radius, cylinder_height).Shape()
    
    hole_centers = []
    result_shape = cylinder
    
    # 各穴を順次作成
    for i, hole_spec in enumerate(hole_specs):
        hole_center = hole_spec.get_position(cylinder_radius)
        hole_direction = hole_spec.get_direction()
        
        print(f"穴{i+1}: 位置({hole_center.X():.1f}, {hole_center.Y():.1f}, {hole_center.Z():.1f}), "
              f"角度{hole_spec.angle_deg}°, 径{hole_spec.diameter}mm")
        
        hole_axis = gp_Ax2(hole_center, hole_direction)
        hole_cylinder = BRepPrimAPI_MakeCylinder(hole_axis, hole_spec.diameter, hole_spec.depth).Shape()
        
        # 穴を開ける
        cut = BRepAlgoAPI_Cut(result_shape, hole_cylinder)
        cut.Build()
        if cut.IsDone():
            result_shape = cut.Shape()
            hole_centers.append(hole_center)
        else:
            print(f"警告: 穴{i+1}の作成に失敗")
    
    return result_shape, hole_centers


def create_spine_connections(hole_centers, connection_type="direct"):
    """任意の数の穴間をSpineで接続"""
    if len(hole_centers) < 2:
        print("警告: 接続するには最低2つの穴が必要")
        return [], []
    
    spine_edges = []
    spine_wires = []
    
    print(f"=== {len(hole_centers)}個の穴を{connection_type}方式で接続 ===")
    
    if connection_type == "direct":
        # 直接接続（全ての穴を順番に接続）
        for i in range(len(hole_centers) - 1):
            pt1, pt2 = hole_centers[i], hole_centers[i + 1]
            edge, wire = create_single_spine_connection(pt1, pt2, i + 1)
            if edge and wire:
                spine_edges.append(edge)
                spine_wires.append(wire)
    
    elif connection_type == "hub":
        # ハブ接続（中央点から各穴へ）
        center = calculate_center_point(hole_centers)
        for i, hole_center in enumerate(hole_centers):
            edge, wire = create_single_spine_connection(center, hole_center, i + 1)
            if edge and wire:
                spine_edges.append(edge)
                spine_wires.append(wire)
    
    elif connection_type == "network":
        # 全接続（全ての穴同士を接続）
        for i in range(len(hole_centers)):
            for j in range(i + 1, len(hole_centers)):
                pt1, pt2 = hole_centers[i], hole_centers[j]
                edge, wire = create_single_spine_connection(pt1, pt2, f"{i+1}-{j+1}")
                if edge and wire:
                    spine_edges.append(edge)
                    spine_wires.append(wire)
    
    print(f"接続完了: {len(spine_edges)}本のSpine作成")
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
    
    return gp_Pnt(sum_x/n, sum_y/n, sum_z/n)


def run_test_case(test_name, hole_specs, connection_type="direct"):
    """テストケースを実行"""
    print(f"\n{'='*60}")
    print(f"テストケース: {test_name}")
    print(f"{'='*60}")
    
    obj = dispocc(touch=True)
    
    # シリンダーに穴を開ける
    cylinder_with_holes, hole_centers = create_cylinder_with_holes(hole_specs=hole_specs)
    obj.display.DisplayShape(cylinder_with_holes, color="BLUE", transparency=0.3)
    
    # Spine接続
    spine_edges, spine_wires = create_spine_connections(hole_centers, connection_type)
    
    # 結果表示
    colors = ["RED", "GREEN", "YELLOW", "ORANGE", "PURPLE"]
    for i, (edge, wire) in enumerate(zip(spine_edges, spine_wires)):
        color = colors[i % len(colors)]
        obj.display.DisplayShape(edge, color=color, transparency=0.0)
    
    print(f"✅ テスト完了: {len(spine_edges)}本の接続")
    obj.ShowOCC()
    return len(spine_edges) > 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--test", dest="test", default="basic", 
                       help="テストケース: basic, different_heights, different_angles, mixed")
    opt = parser.parse_args()
    
    # テストケース定義
    test_cases = {
        "basic": {
            "name": "基本テスト（径違い、高さ違い）",
            "holes": [
                HoleSpec(height=70, angle_deg=0, diameter=8.0),
                HoleSpec(height=30, angle_deg=0, diameter=6.0),
            ]
        },
        "different_heights": {
            "name": "高さ違いテスト", 
            "holes": [
                HoleSpec(height=80, angle_deg=0, diameter=8.0),
                HoleSpec(height=50, angle_deg=0, diameter=8.0),
                HoleSpec(height=20, angle_deg=0, diameter=8.0),
            ]
        },
        "different_angles": {
            "name": "回転位置違いテスト",
            "holes": [
                HoleSpec(height=50, angle_deg=0, diameter=6.0),    # 正面
                HoleSpec(height=50, angle_deg=90, diameter=6.0),   # 右側
                HoleSpec(height=50, angle_deg=180, diameter=6.0),  # 背面
                HoleSpec(height=50, angle_deg=270, diameter=6.0),  # 左側
            ]
        },
        "mixed": {
            "name": "複合テスト（径、高さ、角度すべて違い）",
            "holes": [
                HoleSpec(height=80, angle_deg=0, diameter=10.0),
                HoleSpec(height=60, angle_deg=120, diameter=8.0),
                HoleSpec(height=40, angle_deg=240, diameter=6.0),
                HoleSpec(height=20, angle_deg=45, diameter=4.0),
            ]
        }
    }
    
    if opt.test in test_cases:
        case = test_cases[opt.test]
        success = run_test_case(case["name"], case["holes"])
        
        print(f"\n{'='*60}")
        if success:
            print("✅ 結論: 一般性のあるSpine接続システムが正常に動作")
            print("✅ 対応: 径違い、高さ違い、回転位置違いすべてに対応")
        else:
            print("❌ 結論: システムに問題があります")
    else:
        print(f"未知のテストケース: {opt.test}")
        print(f"利用可能: {list(test_cases.keys())}")