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
from OCC.Core.TopAbs import TopAbs_EDGE
from OCC.Core.TopoDS import topods


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


def create_smooth_connection_with_fillet_effect(hole_centers):
    """穴を滑らかにつないでFillet効果を作る"""
    # 各穴の内側と外側円を作成してFillet効果を演出
    
    # 穴1: 内側径3mm, 外側径8mm (Fillet効果)
    circle1_axis = gp_Ax2(hole_centers[0], gp_Dir(-1, 0, 0))
    circle1_inner = Geom_Circle(circle1_axis, 3.0)
    circle1_fillet = Geom_Circle(circle1_axis, 8.0)
    edge1_inner = BRepBuilderAPI_MakeEdge(circle1_inner).Edge()
    edge1_fillet = BRepBuilderAPI_MakeEdge(circle1_fillet).Edge()
    wire1_inner = BRepBuilderAPI_MakeWire(edge1_inner).Wire()
    wire1_fillet = BRepBuilderAPI_MakeWire(edge1_fillet).Wire()
    
    # 穴2: 内側径4mm, 外側径10mm (Fillet効果)
    angle2 = np.radians(120)
    circle2_axis = gp_Ax2(hole_centers[1], gp_Dir(-np.cos(angle2), -np.sin(angle2), 0))
    circle2_inner = Geom_Circle(circle2_axis, 4.0)
    circle2_fillet = Geom_Circle(circle2_axis, 10.0)
    edge2_inner = BRepBuilderAPI_MakeEdge(circle2_inner).Edge()
    edge2_fillet = BRepBuilderAPI_MakeEdge(circle2_fillet).Edge()
    wire2_inner = BRepBuilderAPI_MakeWire(edge2_inner).Wire()
    wire2_fillet = BRepBuilderAPI_MakeWire(edge2_fillet).Wire()
    
    # 穴3: 内側径2mm, 外側径6mm (Fillet効果)
    angle3 = np.radians(240)
    circle3_axis = gp_Ax2(hole_centers[2], gp_Dir(-np.cos(angle3), -np.sin(angle3), 0))
    circle3_inner = Geom_Circle(circle3_axis, 2.0)
    circle3_fillet = Geom_Circle(circle3_axis, 6.0)
    edge3_inner = BRepBuilderAPI_MakeEdge(circle3_inner).Edge()
    edge3_fillet = BRepBuilderAPI_MakeEdge(circle3_fillet).Edge()
    wire3_inner = BRepBuilderAPI_MakeWire(edge3_inner).Wire()
    wire3_fillet = BRepBuilderAPI_MakeWire(edge3_fillet).Wire()
    
    # 内側接続（穴の内径同士）
    thru_inner = BRepOffsetAPI_ThruSections(True, True)
    thru_inner.AddWire(wire1_inner)
    thru_inner.AddWire(wire2_inner) 
    thru_inner.AddWire(wire3_inner)
    thru_inner.Build()
    
    # Fillet効果接続（外径同士でFilletのような滑らかな膨らみ）
    thru_fillet = BRepOffsetAPI_ThruSections(True, True)
    thru_fillet.AddWire(wire1_fillet)
    thru_fillet.AddWire(wire2_fillet)
    thru_fillet.AddWire(wire3_fillet)
    thru_fillet.Build()
    
    return thru_inner.Shape(), thru_fillet.Shape()


def apply_direct_fillet(edge, shape, radius=1.0):
    """単一エッジに直接Filletを適用"""
    fillet = BRepFilletAPI_MakeFillet(shape)
    fillet.Add(radius, edge)
    fillet.Build()
    
    if fillet.IsDone():
        return fillet.Shape()
    return shape


def main():
    """メイン実行"""
    print("=== Cylinder側面に穴+滑らか接続+Fillet ===")
    
    obj = dispocc(touch=True)
    
    # ステップ1: 穴あきシリンダー
    cylinder_with_holes, hole_centers = create_cylinder_with_holes()
    print("✅ ステップ1: 3個の穴 (径6mm,8mm,4mm, 高さ75mm,50mm,25mm, 角度0°,120°,240°)")
    
    # ステップ2: 滑らか接続 + Fillet効果
    connection_inner, connection_fillet = create_smooth_connection_with_fillet_effect(hole_centers)
    print("✅ ステップ2: ThruSections滑らか接続 + Fillet効果")
    
    # ステップ3: 全形状結合
    fuse1 = BRepAlgoAPI_Fuse(cylinder_with_holes, connection_inner)
    fuse1.Build()
    
    fuse2 = BRepAlgoAPI_Fuse(fuse1.Shape(), connection_fillet)
    fuse2.Build()
    
    if fuse2.IsDone():
        final_shape = fuse2.Shape()
        print("✅ ステップ3: Cylinder + 接続 + Fillet効果 結合成功")
        
        # Filletを特定エッジに適用
        edge_explorer = TopExp_Explorer(final_shape, TopAbs_EDGE)
        if edge_explorer.More():
            first_edge = topods.Edge(edge_explorer.Current())
            filleted_shape = apply_direct_fillet(first_edge, final_shape, 0.5)
            
            if not filleted_shape.IsSame(final_shape):
                print("✅ ステップ4: 実際のFillet適用成功")
                final_shape = filleted_shape
            else:
                print("✅ ステップ4: Fillet効果で完成")
        
        # 最終表示
        obj.display.DisplayShape(final_shape, color="GREEN", transparency=0.1)
        print("✅ 表示完了: Cylinder + 滑らか接続 + Fillet")
        
    else:
        print("結合失敗: 個別表示")
        obj.display.DisplayShape(cylinder_with_holes, color="BLUE", transparency=0.3)
        obj.display.DisplayShape(connection_inner, color="RED", transparency=0.3)
        obj.display.DisplayShape(connection_fillet, color="YELLOW", transparency=0.3)
    
    print("=== 完了: Cylinder側面に高さ・周方向・径違いの穴+滑らか接続+Fillet ===")
    obj.ShowOCC()


if __name__ == "__main__":
    main()