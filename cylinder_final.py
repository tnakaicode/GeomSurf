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
    
    hole_centers = [hole1_center, hole2_center, hole3_center]
    return cut3.Shape(), hole_centers


def create_fillet_connection(hole_centers):
    """Filletのような滑らかな接続を作る"""
    # 各穴に対して段階的な径を持つ複数の断面を作成
    
    # 穴1の段階的断面 (径3mm → 5mm → 7mm)
    circle1_axis = gp_Ax2(hole_centers[0], gp_Dir(-1, 0, 0))
    circle1_1 = Geom_Circle(circle1_axis, 3.0)  # 穴径
    circle1_2 = Geom_Circle(circle1_axis, 5.0)  # 中間径
    circle1_3 = Geom_Circle(circle1_axis, 7.0)  # 外径
    
    # 穴2の段階的断面 (径4mm → 6mm → 9mm)
    angle2 = np.radians(120)
    circle2_axis = gp_Ax2(hole_centers[1], gp_Dir(-np.cos(angle2), -np.sin(angle2), 0))
    circle2_1 = Geom_Circle(circle2_axis, 4.0)  # 穴径
    circle2_2 = Geom_Circle(circle2_axis, 6.0)  # 中間径
    circle2_3 = Geom_Circle(circle2_axis, 9.0)  # 外径
    
    # 穴3の段階的断面 (径2mm → 4mm → 6mm)
    angle3 = np.radians(240)
    circle3_axis = gp_Ax2(hole_centers[2], gp_Dir(-np.cos(angle3), -np.sin(angle3), 0))
    circle3_1 = Geom_Circle(circle3_axis, 2.0)  # 穴径
    circle3_2 = Geom_Circle(circle3_axis, 4.0)  # 中間径
    circle3_3 = Geom_Circle(circle3_axis, 6.0)  # 外径
    
    # 各断面のワイヤを作成
    def make_wire(circle):
        edge = BRepBuilderAPI_MakeEdge(circle).Edge()
        return BRepBuilderAPI_MakeWire(edge).Wire()
    
    # 内側接続 (穴径同士)
    thru_inner = BRepOffsetAPI_ThruSections(True, True)
    thru_inner.AddWire(make_wire(circle1_1))
    thru_inner.AddWire(make_wire(circle2_1))
    thru_inner.AddWire(make_wire(circle3_1))
    thru_inner.Build()
    
    # 中間接続 (Filletの膨らみ効果)
    thru_middle = BRepOffsetAPI_ThruSections(True, True)
    thru_middle.AddWire(make_wire(circle1_2))
    thru_middle.AddWire(make_wire(circle2_2))
    thru_middle.AddWire(make_wire(circle3_2))
    thru_middle.Build()
    
    # 外側接続 (Filletの外縁効果)
    thru_outer = BRepOffsetAPI_ThruSections(True, True)
    thru_outer.AddWire(make_wire(circle1_3))
    thru_outer.AddWire(make_wire(circle2_3))
    thru_outer.AddWire(make_wire(circle3_3))
    thru_outer.Build()
    
    return thru_inner.Shape(), thru_middle.Shape(), thru_outer.Shape()


def main():
    """メイン実行"""
    print("=== Cylinder側面に穴+滑らか接続+Fillet ===")
    
    obj = dispocc(touch=True)
    
    # ステップ1: 穴あきシリンダー
    cylinder_with_holes, hole_centers = create_cylinder_with_holes()
    print("✅ ステップ1完了: 3個の穴 (径6mm,8mm,4mm, 高さ75mm,50mm,25mm, 角度0°,120°,240°)")
    
    # ステップ2: Filletのような滑らか接続
    connection_inner, connection_middle, connection_outer = create_fillet_connection(hole_centers)
    print("✅ ステップ2完了: ThruSectionsで段階的Fillet効果")
    
    # ステップ3: 全形状結合
    # 内側 + 中間
    fuse1 = BRepAlgoAPI_Fuse(connection_inner, connection_middle)
    fuse1.Build()
    
    # 内側+中間 + 外側
    fuse2 = BRepAlgoAPI_Fuse(fuse1.Shape(), connection_outer)
    fuse2.Build()
    
    # シリンダー + 接続部分
    fuse3 = BRepAlgoAPI_Fuse(cylinder_with_holes, fuse2.Shape())
    fuse3.Build()
    
    if fuse3.IsDone():
        final_shape = fuse3.Shape()
        print("✅ ステップ3完了: Cylinder + 段階的接続によるFillet効果")
        
        # 最終表示
        obj.display.DisplayShape(final_shape, color="GREEN", transparency=0.1)
        print("✅ 表示完了: CylinderとのFilletが適用された最終形状")
        
    else:
        print("結合失敗: 個別表示")
        obj.display.DisplayShape(cylinder_with_holes, color="BLUE", transparency=0.3)
        obj.display.DisplayShape(connection_inner, color="RED", transparency=0.3)
        obj.display.DisplayShape(connection_middle, color="YELLOW", transparency=0.3)
        obj.display.DisplayShape(connection_outer, color="ORANGE", transparency=0.3)
    
    print("=== 完了: Cylinder側面に高さ・周方向・径違いの穴 + 滑らか接続 + Fillet ===")
    obj.ShowOCC()


if __name__ == "__main__":
    main()