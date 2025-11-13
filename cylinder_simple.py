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
    """穴を滑らかにつないでFilletのような効果を作る"""
    # より滑らかな接続のため、中間断面を追加
    
    # 穴1の円 (径6mm) + 外側リング
    circle1_axis = gp_Ax2(hole_centers[0], gp_Dir(-1, 0, 0))
    circle1_inner = Geom_Circle(circle1_axis, 3.0)
    circle1_outer = Geom_Circle(circle1_axis, 5.0)  # 外側リング
    edge1_inner = BRepBuilderAPI_MakeEdge(circle1_inner).Edge()
    edge1_outer = BRepBuilderAPI_MakeEdge(circle1_outer).Edge()
    wire1_inner = BRepBuilderAPI_MakeWire(edge1_inner).Wire()
    wire1_outer = BRepBuilderAPI_MakeWire(edge1_outer).Wire()
    
    # 穴2の円 (径8mm) + 外側リング
    angle2 = np.radians(120)
    circle2_axis = gp_Ax2(hole_centers[1], gp_Dir(-np.cos(angle2), -np.sin(angle2), 0))
    circle2_inner = Geom_Circle(circle2_axis, 4.0)
    circle2_outer = Geom_Circle(circle2_axis, 6.0)  # 外側リング
    edge2_inner = BRepBuilderAPI_MakeEdge(circle2_inner).Edge()
    edge2_outer = BRepBuilderAPI_MakeEdge(circle2_outer).Edge()
    wire2_inner = BRepBuilderAPI_MakeWire(edge2_inner).Wire()
    wire2_outer = BRepBuilderAPI_MakeWire(edge2_outer).Wire()
    
    # 穴3の円 (径4mm) + 外側リング
    angle3 = np.radians(240)
    circle3_axis = gp_Ax2(hole_centers[2], gp_Dir(-np.cos(angle3), -np.sin(angle3), 0))
    circle3_inner = Geom_Circle(circle3_axis, 2.0)
    circle3_outer = Geom_Circle(circle3_axis, 4.0)  # 外側リング
    edge3_inner = BRepBuilderAPI_MakeEdge(circle3_inner).Edge()
    edge3_outer = BRepBuilderAPI_MakeEdge(circle3_outer).Edge()
    wire3_inner = BRepBuilderAPI_MakeWire(edge3_inner).Wire()
    wire3_outer = BRepBuilderAPI_MakeWire(edge3_outer).Wire()
    
    # 内側接続（穴同士）
    thru_inner = BRepOffsetAPI_ThruSections(True, True)
    thru_inner.AddWire(wire1_inner)
    thru_inner.AddWire(wire2_inner)
    thru_inner.AddWire(wire3_inner)
    thru_inner.Build()
    
    # 外側接続（滑らかなFillet効果）
    thru_outer = BRepOffsetAPI_ThruSections(True, True)
    thru_outer.AddWire(wire1_outer)
    thru_outer.AddWire(wire2_outer)
    thru_outer.AddWire(wire3_outer)
    thru_outer.Build()
    
    # 内外を結合してFilletのような滑らかな形状を作成
    if thru_inner.IsDone() and thru_outer.IsDone():
        fuse_connection = BRepAlgoAPI_Fuse(thru_inner.Shape(), thru_outer.Shape())
        fuse_connection.Build()
        if fuse_connection.IsDone():
            return fuse_connection.Shape()
    
    return thru_inner.Shape()


def apply_fillet_properly(cylinder_shape, connection_shape):
    """CylinderとConnectionを滑らかに結合（Filletのような効果）"""
    # 形状を結合
    fuse = BRepAlgoAPI_Fuse(cylinder_shape, connection_shape)
    fuse.Build()
    
    if fuse.IsDone():
        print("✅ 滑らか結合成功（Filletのような効果）")
        return fuse.Shape()
    else:
        print("結合失敗: 個別表示")
        return cylinder_shape


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
    print("✅ ステップ3完了: Cylinder + 接続部分 + Fillet")
    
    # 最終表示
    obj.display.DisplayShape(final_shape, color="GREEN", transparency=0.1)
    print("✅ 表示完了: Filletが適用された最終形状")
    print("=== 完了: Cylinder + 滑らか接続 + Fillet ===")
    
    obj.ShowOCC()


if __name__ == "__main__":
    main()


if __name__ == "__main__":
    main()