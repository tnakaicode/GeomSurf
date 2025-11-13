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
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeWire, BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeEdge
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_ThruSections
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
from OCC.Core.BRepFilletAPI import BRepFilletAPI_MakeFillet
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Edge
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_EDGE
from OCC.Core.Geom import Geom_Circle
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
    
    # 穴3: 側面の穴1（X軸方向）
    hole3_axis = gp_Ax2(gp_Pnt(0, 0, 60), gp_Dir(1, 0, 0))
    hole3_cylinder = BRepPrimAPI_MakeCylinder(hole3_axis, 6.0, 60.0).Shape()
    
    # 穴4: 側面の穴2（Y軸方向）
    hole4_axis = gp_Ax2(gp_Pnt(0, 0, 40), gp_Dir(0, 1, 0))
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
    
    return cylinder_with_holes, gp_Pnt(20, 0, 80), gp_Pnt(-20, 0, 20)


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


def apply_fillet(shape1, shape2, radius=5.0):
    """2つの形状を結合してFilletを適用"""
    # 形状を結合
    fuse = BRepAlgoAPI_Fuse(shape1, shape2)
    fuse.Build()
    fused_shape = fuse.Shape()
    
    # Filletを適用
    fillet = BRepFilletAPI_MakeFillet(fused_shape)
    
    # エッジを探してFilletを適用
    explorer = TopExp_Explorer(fused_shape, TopAbs_EDGE)
    edge_count = 0
    while explorer.More() and edge_count < 20:  # より多くのエッジにFilletを適用
        edge = explorer.Current()
        try:
            fillet.Add(radius, edge)
            edge_count += 1
        except:
            pass  # Filletを適用できないエッジはスキップ
        explorer.Next()
    
    try:
        fillet.Build()
        return fillet.Shape()
    except:
        return fused_shape  # Filletが失敗した場合は結合された形状を返す


def apply_cylinder_fillet(cylinder_shape, radius=3.0):
    """シリンダー自体にFilletを適用"""
    fillet = BRepFilletAPI_MakeFillet(cylinder_shape)
    
    # エッジを探してFilletを適用
    explorer = TopExp_Explorer(cylinder_shape, TopAbs_EDGE)
    edge_count = 0
    while explorer.More() and edge_count < 15:
        edge = explorer.Current()
        try:
            fillet.Add(radius, edge)
            edge_count += 1
        except:
            pass  # Filletを適用できないエッジはスキップ
        explorer.Next()
    
    try:
        fillet.Build()
        return fillet.Shape()
    except:
        return cylinder_shape  # Filletが失敗した場合は元の形状を返す


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--radius", dest="radius", default=5.0, type=float)
    parser.add_argument("--cylinder_radius", dest="cylinder_radius", default=3.0, type=float)
    opt = parser.parse_args()
    
    obj = dispocc(touch=True)
    
    # ステップ1: シリンダーに穴を開ける
    cylinder_with_holes, hole1_center, hole2_center = create_cylinder_with_holes()
    
    # ステップ2: シリンダーにFilletを適用
    cylinder_filleted = apply_cylinder_fillet(cylinder_with_holes, opt.cylinder_radius)
    obj.display.DisplayShape(cylinder_filleted, color="BLUE", transparency=0.3)
    
    # ステップ3: ThruSectionを作成
    thru_section = create_thru_section(hole1_center, hole2_center)
    obj.display.DisplayShape(thru_section, color="RED", transparency=0.3)
    
    # ステップ4: Filletを適用した最終形状を作成
    final_shape = apply_fillet(cylinder_filleted, thru_section, opt.radius)
    obj.display.DisplayShape(final_shape, color="GREEN", transparency=0.8)
    
    print(f"シリンダーに4つの穴（上下2つ、側面2つ）を開け、ThruSectionで接続し、")
    print(f"シリンダーFillet半径{opt.cylinder_radius}、結合Fillet半径{opt.radius}を適用しました。")
    
    obj.ShowOCC()