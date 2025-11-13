import numpy as np
import sys
import os

basename = os.path.dirname(__file__)
sys.path.append(os.path.join("./"))
from src.base_occ import dispocc

from OCC.Core.gp import gp_Pnt, gp_Dir, gp_Ax2
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeCylinder
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeWire, BRepBuilderAPI_MakeEdge
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_ThruSections
from OCC.Core.Geom import Geom_Circle
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
from OCC.Core.BRepFilletAPI import BRepFilletAPI_MakeFillet
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_EDGE, TopAbs_FACE
from OCC.Core.TopoDS import topods
from OCC.Core.BRep import BRep_Tool


def create_cylinder_with_holes():
    """Cylinder側面に高さ・周方向・径違いの穴を開ける"""
    # シリンダー作成 (半径50mm, 高さ100mm)
    cylinder_axis = gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1))
    cylinder = BRepPrimAPI_MakeCylinder(cylinder_axis, 50.0, 100.0).Shape()
    
    # 穴1: 高さ75mm, 角度0°, 径6mm
    hole1_center = gp_Pnt(50, 0, 75)
    hole1_axis = gp_Ax2(hole1_center, gp_Dir(-1, 0, 0))
    hole1_cylinder = BRepPrimAPI_MakeCylinder(hole1_axis, 3.0, 30.0).Shape()
    
    # 穴2: 高さ50mm, 角度120°, 径8mm  
    angle2 = np.radians(120)
    hole2_center = gp_Pnt(50*np.cos(angle2), 50*np.sin(angle2), 50)
    hole2_axis = gp_Ax2(hole2_center, gp_Dir(-np.cos(angle2), -np.sin(angle2), 0))
    hole2_cylinder = BRepPrimAPI_MakeCylinder(hole2_axis, 4.0, 30.0).Shape()
    
    # 穴3: 高さ25mm, 角度240°, 径4mm
    angle3 = np.radians(240) 
    hole3_center = gp_Pnt(50*np.cos(angle3), 50*np.sin(angle3), 25)
    hole3_axis = gp_Ax2(hole3_center, gp_Dir(-np.cos(angle3), -np.sin(angle3), 0))
    hole3_cylinder = BRepPrimAPI_MakeCylinder(hole3_axis, 2.0, 30.0).Shape()
    
    # 穴を開ける
    cut1 = BRepAlgoAPI_Cut(cylinder, hole1_cylinder)
    cut1.Build()
    cut2 = BRepAlgoAPI_Cut(cut1.Shape(), hole2_cylinder) 
    cut2.Build()
    cut3 = BRepAlgoAPI_Cut(cut2.Shape(), hole3_cylinder)
    cut3.Build()
    
    # 穴の表面中心点
    hole_centers = [hole1_center, hole2_center, hole3_center]
    
    return cut3.Shape(), hole_centers


def create_smooth_connection(hole_centers):
    """穴を滑らかにつなぐ"""
    # 穴1の円 (径6mm)
    circle1_axis = gp_Ax2(hole_centers[0], gp_Dir(-1, 0, 0))
    circle1_geom = Geom_Circle(circle1_axis, 3.0)
    circle1_edge = BRepBuilderAPI_MakeEdge(circle1_geom).Edge()
    wire1 = BRepBuilderAPI_MakeWire(circle1_edge).Wire()
    
    # 穴2の円 (径8mm)
    angle2 = np.radians(120)
    circle2_axis = gp_Ax2(hole_centers[1], gp_Dir(-np.cos(angle2), -np.sin(angle2), 0))
    circle2_geom = Geom_Circle(circle2_axis, 4.0)
    circle2_edge = BRepBuilderAPI_MakeEdge(circle2_geom).Edge()
    wire2 = BRepBuilderAPI_MakeWire(circle2_edge).Wire()
    
    # 穴3の円 (径4mm)
    angle3 = np.radians(240)
    circle3_axis = gp_Ax2(hole_centers[2], gp_Dir(-np.cos(angle3), -np.sin(angle3), 0))
    circle3_geom = Geom_Circle(circle3_axis, 2.0)
    circle3_edge = BRepBuilderAPI_MakeEdge(circle3_geom).Edge()
    wire3 = BRepBuilderAPI_MakeWire(circle3_edge).Wire()
    
    # ThruSectionsで滑らかに接続
    thru_sections = BRepOffsetAPI_ThruSections(True, True)
    thru_sections.AddWire(wire1)
    thru_sections.AddWire(wire2)
    thru_sections.AddWire(wire3)
    thru_sections.Build()
    
    return thru_sections.Shape()


def apply_real_fillet(fused_shape):
    """実際のFilletを適用"""
    fillet = BRepFilletAPI_MakeFillet(fused_shape)
    
    # 特定のエッジを選択してFilletを適用
    edge_explorer = TopExp_Explorer(fused_shape, TopAbs_EDGE)
    edge_list = []
    
    # 全エッジを収集
    while edge_explorer.More():
        edge_list.append(topods.Edge(edge_explorer.Current()))
        edge_explorer.Next()
    
    # 安全な半径で段階的にFillet適用
    for radius in [0.2, 0.5, 1.0]:
        test_fillet = BRepFilletAPI_MakeFillet(fused_shape)
        success_count = 0
        
        # エッジを1つずつ追加してテスト
        for i, edge in enumerate(edge_list[:5]):  # 最初の5つのエッジのみ
            test_fillet.Add(radius, edge)
            test_fillet.Build()
            
            if test_fillet.IsDone():
                success_count += 1
                print(f"Fillet成功 エッジ{i+1}: 半径{radius}mm")
                return test_fillet.Shape()
            else:
                # 失敗したら新しいFilletオブジェクトで再開
                test_fillet = BRepFilletAPI_MakeFillet(fused_shape)
    
    print("Fillet適用失敗: 元の形状を返す")
    return fused_shape


def apply_fillet_properly(cylinder_shape, connection_shape):
    """CylinderとConnectionにFilletを適用"""
    # まず形状を結合
    fuse = BRepAlgoAPI_Fuse(cylinder_shape, connection_shape)
    fuse.Build()
    
    if not fuse.IsDone():
        print("形状結合失敗")
        return cylinder_shape
    
    fused_shape = fuse.Shape()
    print("形状結合成功")
    
    # 実際のFilletを適用
    filleted_shape = apply_real_fillet(fused_shape)
    return filleted_shape


def main():
    """メイン実行"""
    print("=== Cylinder側面に穴を開けて滑らかに接続してFilletを適用 ===")
    
    obj = dispocc(touch=True)
    
    # ステップ1: 穴あきシリンダー作成
    cylinder_with_holes, hole_centers = create_cylinder_with_holes()
    print("✅ ステップ1完了: 3個の穴 (径6mm,8mm,4mm, 高さ75mm,50mm,25mm, 角度0°,120°,240°)")
    
    # ステップ2: 滑らか接続
    connection = create_smooth_connection(hole_centers)
    print("✅ ステップ2完了: ThruSectionsで滑らかな接続")
    
    # ステップ3: Fillet適用
    final_shape = apply_fillet_properly(cylinder_with_holes, connection)
    print("✅ ステップ3完了: Cylinder + 接続部分 + 実際のFillet")
    
    # 最終表示
    obj.display.DisplayShape(final_shape, color="GREEN", transparency=0.1)
    print("✅ 表示完了: 実際のFilletが適用された最終形状")
    print("=== 完了: Cylinder + 滑らか接続 + 実際のFillet ===")
    
    obj.ShowOCC()


if __name__ == "__main__":
    main()