import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import argparse
from linecache import getline, clearcache, updatecache

basename = os.path.dirname(__file__)

sys.path.append(os.path.join("./"))
from base_occ import dispocc

# 分離した関数モジュールをインポート
from shape_generators import *
from analysis_utils import *
from shape_helpers import *
from wire_generators import *
from complex_wire_generators import *
from solid_generators import *

import logging

logging.getLogger("matplotlib").setLevel(logging.ERROR)

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir, gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Trsf
from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_WIRE
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Compound
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Core.TColgp import TColgp_Array1OfPnt
from OCC.Core.Geom import Geom_Circle, Geom_BSplineCurve
from OCC.Core.GProp import GProp_GProps
from OCC.Core.GeomAPI import GeomAPI_PointsToBSpline
from OCC.Core.BRepBuilderAPI import (
    BRepBuilderAPI_MakeEdge,
    BRepBuilderAPI_MakeWire,
    BRepBuilderAPI_MakeFace,
    BRepBuilderAPI_Transform,
)
from OCC.Core.BRepPrimAPI import (
    BRepPrimAPI_MakePrism,
    BRepPrimAPI_MakeCylinder,
    BRepPrimAPI_MakeSphere,
)
from OCC.Core.BRepGProp import brepgprop
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Cut
from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRepBndLib import brepbndlib
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_ThruSections
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCCUtils.Common import vertex2pnt
from OCCUtils.Construct import make_box, make_line, make_wire, make_edge
from OCCUtils.Construct import make_plane, make_polygon
from OCCUtils.Construct import point_to_vector, vector_to_point
from OCCUtils.Construct import dir_to_vec, vec_to_dir


class YZPlaneSweep(dispocc):
    def __init__(self):
        """YZ平面上のWireをX軸方向にSweepするクラス（描画特化）"""
        dispocc.__init__(self)
        self.axis = gp_Ax3()

    def create_single_cross_section_demo(self):
        """単一の断面作成デモ"""
        print("\n=== 単一断面作成デモ ===")

        # テスト用の円形Wire
        test_wire = create_circle_wire_on_yz_plane(
            center_y=0.0, center_z=0.0, radius=15.0
        )
        test_solid = create_lofted_solid_from_wire(test_wire, length=80.0)

        # テストSolidを表示
        transform_test = gp_Trsf()
        transform_test.SetTranslation(gp_Vec(0.0, 0.0, -100.0))
        test_solid_transformed = BRepBuilderAPI_Transform(
            test_solid, transform_test
        ).Shape()
        self.display.DisplayShape(
            test_solid_transformed, color="LIGHTBLUE", transparency=0.5
        )

        # 特定の位置での断面を作成
        section_positions = [20.0, 40.0, 60.0]
        section_colors = ["RED", "GREEN", "BLUE"]

        for i, x_pos in enumerate(section_positions):
            try:
                section = get_cross_section_at_x(test_solid_transformed, x_pos)
                color = section_colors[i % len(section_colors)]
                self.display.DisplayShape(section, color=color, transparency=0.0)
                print(f"断面をX={x_pos}で作成（色: {color}）")
            except Exception as e:
                print(f"断面作成失敗 (X={x_pos}): {e}")

    def run_sweep_demo(self):
        """Sweepデモを実行"""
        print("YZ平面上にWireを作成し、X軸方向にSweepします...")

        # 1. 円形のWire（中空構造用）
        print("1. 円形の中空構造を作成...")
        outer_circle, inner_circle = create_concentric_circles(
            outer_radius=25.0, inner_radius=15.0, center_y=0.0, center_z=0.0
        )

        # 中空の円形立体（青の中に黄色の空洞）
        circle_hollow_sweep = create_hollow_solid(
            outer_circle, inner_circle, length=100.0, method="sweep"
        )
        circle_hollow_loft = create_hollow_solid(
            outer_circle, inner_circle, length=100.0, method="loft"
        )

        # 2. 多角形のWire（中空構造用）
        print("2. 多角形の中空構造を作成...")
        outer_polygon, inner_polygon = create_scaled_polygon_pair(scale_factor=0.6)

        # 中空の多角形立体（白の中に緑の空洞）
        polygon_hollow_sweep = create_hollow_solid(
            outer_polygon, inner_polygon, length=100.0, method="sweep"
        )
        polygon_hollow_loft = create_hollow_solid(
            outer_polygon, inner_polygon, length=100.0, method="loft"
        )

        # 3. スムーズ曲線のWire（中空構造用）
        print("3. スムーズ曲線の中空構造を作成...")
        outer_curve = create_smooth_curve_wire_on_yz_plane()
        inner_curve = scale_wire(outer_curve, 0.5, 0.0, 0.0)

        curve_hollow_sweep = create_hollow_solid(
            outer_curve, inner_curve, length=100.0, method="sweep"
        )

        # 形状を表示
        print("形状を表示中...")

        # 元のWireを表示（赤色）
        self.display.DisplayShape(outer_circle, color="RED", transparency=0.0)
        self.display.DisplayShape(inner_circle, color="ORANGE", transparency=0.0)
        self.display.DisplayShape(outer_polygon, color="RED", transparency=0.0)
        self.display.DisplayShape(inner_polygon, color="ORANGE", transparency=0.0)

        # 中空立体を表示（異なる位置に配置）
        # 円形中空立体（Sweep版）- 青色
        self.display.DisplayShape(circle_hollow_sweep, color="BLUE", transparency=0.3)

        # 円形中空立体（Loft版）- 水色、Y軸方向に80mmオフセット
        transform1 = gp_Trsf()
        transform1.SetTranslation(gp_Vec(0.0, 80.0, 0.0))
        circle_hollow_loft_transformed = BRepBuilderAPI_Transform(
            circle_hollow_loft, transform1
        ).Shape()
        self.display.DisplayShape(
            circle_hollow_loft_transformed, color="CYAN", transparency=0.3
        )

        # 多角形中空立体（Sweep版）- 白色、Y軸方向に-80mmオフセット
        transform2 = gp_Trsf()
        transform2.SetTranslation(gp_Vec(0.0, -80.0, 0.0))
        polygon_hollow_sweep_transformed = BRepBuilderAPI_Transform(
            polygon_hollow_sweep, transform2
        ).Shape()
        self.display.DisplayShape(
            polygon_hollow_sweep_transformed, color="WHITE", transparency=0.3
        )

        # 多角形中空立体（Loft版）- 薄緑色、Z軸方向に70mmオフセット
        transform3 = gp_Trsf()
        transform3.SetTranslation(gp_Vec(0.0, 0.0, 70.0))
        polygon_hollow_loft_transformed = BRepBuilderAPI_Transform(
            polygon_hollow_loft, transform3
        ).Shape()
        self.display.DisplayShape(
            polygon_hollow_loft_transformed, color="GREEN", transparency=0.3
        )

        # スムーズ曲線中空立体 - 紫色、Z軸方向に-70mmオフセット
        transform4 = gp_Trsf()
        transform4.SetTranslation(gp_Vec(0.0, 0.0, -70.0))
        curve_hollow_sweep_transformed = BRepBuilderAPI_Transform(
            curve_hollow_sweep, transform4
        ).Shape()
        self.display.DisplayShape(
            curve_hollow_sweep_transformed, color="PURPLE", transparency=0.3
        )

        # 断面解析のデモ（中空構造の断面を確認）
        print("4. 中空構造の断面解析を実行...")
        print("円形中空立体の断面解析:")
        circle_sections = analyze_solid_cross_sections(
            circle_hollow_sweep, num_sections=3
        )
        for i, (pos, section) in enumerate(circle_sections):
            if section:
                colors = ["ORANGE", "PURPLE", "CYAN"]
                color = colors[i % len(colors)]
                self.display.DisplayShape(section, color=color, transparency=0.0)

        print("多角形中空立体の断面解析:")
        polygon_sections = analyze_solid_cross_sections(
            polygon_hollow_loft_transformed, num_sections=3
        )
        for i, (pos, section) in enumerate(polygon_sections):
            if section:
                colors = ["MAGENTA", "BROWN", "YELLOW"]
                color = colors[i % len(colors)]
                self.display.DisplayShape(section, color=color, transparency=0.0)

        print("作成完了!")
        print("=== 中空構造 ===")
        print("- 青: 円形中空立体（Sweep版）")
        print("- シアン: 円形中空立体（Loft版、Y軸+80mmオフセット）")
        print("- 白: 多角形中空立体（Sweep版、Y軸-80mmオフセット）")
        print("- 緑: 多角形中空立体（Loft版、Z軸+70mmオフセット）")
        print("- 紫: スムーズ曲線中空立体（Z軸-70mmオフセット）")
        print("- その他の色: 各中空立体の断面形状")
        print("- 赤: 外側Wire（YZ平面上）")
        print("- オレンジ: 内側Wire（YZ平面上）")

        # 座標軸と平面を表示
        self.show_axs_pln()

        # ディスプレイを開始
        self.ShowOCC()

    def run_hollow_structure_demo(self):
        """中空構造専用のデモ"""
        print("\n=== 面白い中空構造専用デモ ===")

        # 複雑で面白い中空構造を作成
        complex_structures = create_complex_hollow_structures()

        # 配置位置の設定
        positions = [
            gp_Vec(0.0, 0.0, 0.0),  # 中央
            gp_Vec(0.0, 120.0, 0.0),  # Y軸 +120mm
            gp_Vec(0.0, -120.0, 0.0),  # Y軸 -120mm
            gp_Vec(0.0, 0.0, 100.0),  # Z軸 +100mm
            gp_Vec(0.0, 0.0, -100.0),  # Z軸 -100mm
        ]

        print("3. 複雑中空構造を配置...")

        for i, (name, structure, color) in enumerate(complex_structures):
            if i < len(positions):
                # 構造体を指定位置に配置
                transform = gp_Trsf()
                transform.SetTranslation(positions[i])
                transformed_structure = BRepBuilderAPI_Transform(
                    structure, transform
                ).Shape()

                # 表示
                self.display.DisplayShape(
                    transformed_structure, color=color, transparency=0.3
                )
                print(f"  {name}: {color}色で配置")

        # 断面解析のデモ（いくつかの構造で）
        print("4. 面白い中空構造の断面解析...")
        if len(complex_structures) > 0:
            # 最初の構造の断面を表示
            first_structure = complex_structures[0][1]
            sections = analyze_solid_cross_sections(first_structure, num_sections=2)
            for i, (pos, section) in enumerate(sections):
                if section:
                    colors = ["ORANGE", "PURPLE"]
                    color = colors[i % len(colors)]
                    self.display.DisplayShape(section, color=color, transparency=0.0)

            # 3番目の構造の断面も表示（存在する場合）
            if len(complex_structures) > 2:
                transform = gp_Trsf()
                transform.SetTranslation(positions[2])
                third_structure_transformed = BRepBuilderAPI_Transform(
                    complex_structures[2][1], transform
                ).Shape()
                sections2 = analyze_solid_cross_sections(
                    third_structure_transformed, num_sections=2
                )
                for i, (pos, section) in enumerate(sections2):
                    if section:
                        colors = ["CYAN", "MAGENTA"]
                        color = colors[i % len(colors)]
                        self.display.DisplayShape(section, color=color, transparency=0.0)

        # 元のWireも少し表示（参考用）
        star_wire = create_star_wire_on_yz_plane(15.0, 8.0, 6)
        heart_wire = create_heart_wire_on_yz_plane(12.0)
        flower_wire = create_flower_wire_on_yz_plane(15.0, 8, 5.0)

        # 参考Wireを表示（小さくして右上に配置）
        reference_transform = gp_Trsf()
        reference_transform.SetTranslation(gp_Vec(0.0, 150.0, 80.0))
        reference_transform.SetScaleFactor(0.5)  # 半分のサイズ

        star_ref = BRepBuilderAPI_Transform(star_wire, reference_transform).Shape()
        heart_ref = BRepBuilderAPI_Transform(heart_wire, reference_transform).Shape()

        self.display.DisplayShape(star_ref, color="RED", transparency=0.0)
        self.display.DisplayShape(heart_ref, color="RED", transparency=0.0)

        print("面白い中空構造デモ完了!")
        print("=== 面白い中空構造 ===")
        for i, (name, _, color) in enumerate(complex_structures):
            print(f"- {color}: {name}")
        print("- その他: 断面形状と参考Wire")

    def run_artistic_shapes_demo(self):
        """芸術的な形状のデモ"""
        print("\n=== 芸術的形状デモ ===")

        # 各種面白い形状を個別に表示
        print("1. 個別形状の作成と表示...")

        # 星形
        star_wire = create_star_wire_on_yz_plane(20.0, 10.0, 8)
        star_solid = sweep_wire_along_x_axis(star_wire, 60.0)

        # ハート形
        heart_wire = create_heart_wire_on_yz_plane(15.0)
        heart_solid = create_lofted_solid_from_wire(heart_wire, 70.0)

        # 花形
        flower_wire = create_flower_wire_on_yz_plane(18.0, 12, 6.0)
        flower_solid = sweep_wire_along_x_axis(flower_wire, 65.0)

        # 渦巻き形（新）
        spiral_wire = create_spiral_wire_on_yz_plane(22.0, 6.0, 2.5)
        spiral_solid = create_lofted_solid_from_wire(spiral_wire, 55.0)

        # 歯車形（新）
        gear_wire = create_gear_wire_on_yz_plane(20.0, 15.0, 16)
        gear_solid = sweep_wire_along_x_axis(gear_wire, 50.0)

        # 波形円（新）
        wave_wire = create_wave_wire_on_yz_plane(8.0, 6.0, 18.0)
        wave_solid = create_lofted_solid_from_wire(wave_wire, 75.0)

        # 形状を配置
        # 星形（金色）- 中央
        self.display.DisplayShape(star_solid, color="GOLD", transparency=0.4)

        # ハート形（ピンク）- Y軸+80mm
        transform1 = gp_Trsf()
        transform1.SetTranslation(gp_Vec(0.0, 80.0, 0.0))
        heart_transformed = BRepBuilderAPI_Transform(heart_solid, transform1).Shape()
        self.display.DisplayShape(heart_transformed, color="PINK", transparency=0.4)

        # 花形（オレンジ）- Y軸-80mm
        transform2 = gp_Trsf()
        transform2.SetTranslation(gp_Vec(0.0, -80.0, 0.0))
        flower_transformed = BRepBuilderAPI_Transform(flower_solid, transform2).Shape()
        self.display.DisplayShape(flower_transformed, color="ORANGE", transparency=0.4)

        # 渦巻き形（シアン）- Z軸+80mm
        transform3 = gp_Trsf()
        transform3.SetTranslation(gp_Vec(0.0, 0.0, 80.0))
        spiral_transformed = BRepBuilderAPI_Transform(spiral_solid, transform3).Shape()
        self.display.DisplayShape(spiral_transformed, color="CYAN", transparency=0.4)

        # 歯車形（青）- Z軸-80mm
        transform4 = gp_Trsf()
        transform4.SetTranslation(gp_Vec(0.0, 0.0, -80.0))
        gear_transformed = BRepBuilderAPI_Transform(gear_solid, transform4).Shape()
        self.display.DisplayShape(gear_transformed, color="BLUE2", transparency=0.4)

        # 波形円（緑）- Y軸+40mm, Z軸+40mm
        transform5 = gp_Trsf()
        transform5.SetTranslation(gp_Vec(0.0, 40.0, 40.0))
        wave_transformed = BRepBuilderAPI_Transform(wave_solid, transform5).Shape()
        self.display.DisplayShape(wave_transformed, color="GREEN", transparency=0.4)

        # 元のWireも表示
        self.display.DisplayShape(star_wire, color="RED", transparency=0.0)
        self.display.DisplayShape(heart_wire, color="RED", transparency=0.0)
        self.display.DisplayShape(flower_wire, color="RED", transparency=0.0)
        self.display.DisplayShape(spiral_wire, color="RED", transparency=0.0)
        self.display.DisplayShape(gear_wire, color="RED", transparency=0.0)
        self.display.DisplayShape(wave_wire, color="RED", transparency=0.0)

        # 断面解析
        print("2. 芸術的形状の断面解析...")
        star_sections = analyze_solid_cross_sections(star_solid, num_sections=2)
        heart_sections = analyze_solid_cross_sections(heart_transformed, num_sections=2)
        
        for i, (pos, section) in enumerate(star_sections):
            if section:
                colors = ["ORANGE", "PURPLE"]
                color = colors[i % len(colors)]
                self.display.DisplayShape(section, color=color, transparency=0.0)
        
        for i, (pos, section) in enumerate(heart_sections):
            if section:
                colors = ["CYAN", "MAGENTA"]
                color = colors[i % len(colors)]
                self.display.DisplayShape(section, color=color, transparency=0.0)

        print("芸術的形状デモ完了!")
        print("- 金: 8頂点星形立体（中央）")
        print("- ピンク: ハート形立体（Y軸+80mm）")
        print("- オレンジ: 12花びら花形立体（Y軸-80mm）")
        print("- シアン: 渦巻き形立体（Z軸+80mm）")
        print("- 青: 歯車形立体（Z軸-80mm）")
        print("- 緑: 波形円立体（Y軸+40mm, Z軸+40mm）")
        print("- 赤: 元のWire（YZ平面上）")
        print("- その他: 断面形状")

    def create_cylinder_wire_on_yz_plane(self, radius=30.0, center_y=0.0, center_z=0.0):
        """
        YZ平面上に円柱用の円形Wireを作成

        Parameters:
        radius: float - 円の半径
        center_y: float - 中心のY座標
        center_z: float - 中心のZ座標

        Returns:
        TopoDS_Wire: 作成された円形Wire
        """
        return self.create_circle_wire_on_yz_plane(center_y, center_z, radius)

    def create_rectangular_wire_on_yz_plane(
        self, width=60.0, height=40.0, center_y=0.0, center_z=0.0
    ):
        """
        YZ平面上に四角柱用の矩形Wireを作成

        Parameters:
        width: float - Y軸方向の幅
        height: float - Z軸方向の高さ
        center_y: float - 中心のY座標
        center_z: float - 中心のZ座標

        Returns:
        TopoDS_Wire: 作成された矩形Wire
        """
        # 矩形の頂点を定義 (YZ平面上, X=0)
        half_width = width / 2.0
        half_height = height / 2.0

        points = [
            gp_Pnt(0.0, center_y - half_width, center_z - half_height),  # 左下
            gp_Pnt(0.0, center_y + half_width, center_z - half_height),  # 右下
            gp_Pnt(0.0, center_y + half_width, center_z + half_height),  # 右上
            gp_Pnt(0.0, center_y - half_width, center_z + half_height),  # 左上
        ]

        # エッジを作成してWireに追加
        wire_builder = BRepBuilderAPI_MakeWire()

        for i in range(len(points)):
            start_point = points[i]
            end_point = points[(i + 1) % len(points)]  # 最後の点は最初の点に戻る

            edge = BRepBuilderAPI_MakeEdge(start_point, end_point).Edge()
            wire_builder.Add(edge)

        return wire_builder.Wire()

    def create_hollow_solid(self, outer_wire, inner_wire, length=100.0, method="loft"):
        """
        外側と内側のWireから中空のSolidを作成

        Parameters:
        outer_wire: TopoDS_Wire - 外側のワイヤー
        inner_wire: TopoDS_Wire - 内側のワイヤー
        length: float - X軸方向の距離
        method: str - "loft" or "sweep" 作成方法

        Returns:
        TopoDS_Shape: 中空のSolid
        """

        # 外側のSolidを作成
        if method == "loft":
            outer_solid = self.create_lofted_solid_from_wire(outer_wire, length)
        else:
            outer_solid = self.sweep_wire_along_x_axis(outer_wire, length)

        # 内側のSolidを作成
        if method == "loft":
            inner_solid = self.create_lofted_solid_from_wire(inner_wire, length)
        else:
            inner_solid = self.sweep_wire_along_x_axis(inner_wire, length)

        # ブール演算（引き算）で中空構造を作成
        cut_builder = BRepAlgoAPI_Cut(outer_solid, inner_solid)
        cut_builder.Build()

        if not cut_builder.IsDone():
            raise RuntimeError("Failed to create hollow solid")

        return cut_builder.Shape()

    def create_scaled_wire(
        self, original_wire, scale_factor, center_y=0.0, center_z=0.0
    ):
        """
        Wireを指定した中心点でスケールする

        Parameters:
        original_wire: TopoDS_Wire - 元のワイヤー
        scale_factor: float - スケール倍率
        center_y: float - スケール中心のY座標
        center_z: float - スケール中心のZ座標

        Returns:
        TopoDS_Wire: スケールされたワイヤー
        """
        # スケール変換を作成
        scale_transform = gp_Trsf()
        scale_center = gp_Pnt(0.0, center_y, center_z)
        scale_transform.SetScale(scale_center, scale_factor)

        # 変換を適用
        scaled_wire = BRepBuilderAPI_Transform(original_wire, scale_transform).Shape()
        return scaled_wire

    def create_concentric_circles(
        self, outer_radius=20.0, inner_radius=12.0, center_y=0.0, center_z=0.0
    ):
        """
        同心円のWireペアを作成

        Parameters:
        outer_radius: float - 外側の円の半径
        inner_radius: float - 内側の円の半径
        center_y: float - 中心のY座標
        center_z: float - 中心のZ座標

        Returns:
        tuple: (外側Wire, 内側Wire)
        """
        outer_wire = self.create_circle_wire_on_yz_plane(
            center_y, center_z, outer_radius
        )
        inner_wire = self.create_circle_wire_on_yz_plane(
            center_y, center_z, inner_radius
        )
        return outer_wire, inner_wire

    def create_scaled_polygon_pair(self, scale_factor=0.6):
        """
        多角形のWireペア（外側と内側）を作成

        Parameters:
        scale_factor: float - 内側の多角形のスケール倍率

        Returns:
        tuple: (外側Wire, 内側Wire)
        """
        outer_wire = self.create_polygon_wire_on_yz_plane()

        # 多角形の重心を計算
        props = self.get_wire_properties(outer_wire)
        center_y = props["center"][1]
        center_z = props["center"][2]

        # 内側の多角形を作成（スケール）
        inner_wire = self.create_scaled_wire(
            outer_wire, scale_factor, center_y, center_z
        )

        return outer_wire, inner_wire

    def create_star_wire_on_yz_plane(
        self, outer_radius=20.0, inner_radius=10.0, num_points=5
    ):
        """
        YZ平面上に星形のWireを作成

        Parameters:
        outer_radius: float - 外側の頂点までの半径
        inner_radius: float - 内側の谷までの半径
        num_points: int - 星の頂点数

        Returns:
        TopoDS_Wire: 作成された星形Wire
        """
        import math

        points = []
        for i in range(num_points * 2):
            angle = i * math.pi / num_points
            if i % 2 == 0:  # 外側の頂点
                radius = outer_radius
            else:  # 内側の谷
                radius = inner_radius

            y = radius * math.cos(angle)
            z = radius * math.sin(angle)
            points.append(gp_Pnt(0.0, y, z))

        # 最初の点に戻って閉じる
        points.append(points[0])

        # エッジを作成してWireに追加
        wire_builder = BRepBuilderAPI_MakeWire()

        for i in range(len(points) - 1):
            edge = BRepBuilderAPI_MakeEdge(points[i], points[i + 1]).Edge()
            wire_builder.Add(edge)

        return wire_builder.Wire()

    def create_heart_wire_on_yz_plane(self, scale=15.0):
        """
        YZ平面上にハート形のWireを作成

        Parameters:
        scale: float - スケール倍率

        Returns:
        TopoDS_Wire: 作成されたハート形Wire
        """
        import math

        points = []
        num_points = 100

        for i in range(num_points + 1):
            t = 2 * math.pi * i / num_points

            # ハート形の数式
            y = scale * (16 * math.sin(t) ** 3)
            z = scale * (
                13 * math.cos(t)
                - 5 * math.cos(2 * t)
                - 2 * math.cos(3 * t)
                - math.cos(4 * t)
            )

            # スケールを調整
            y = y / 16
            z = z / 16

            points.append(gp_Pnt(0.0, y, z))

        # B-スプライン曲線でスムーズなハートを作成
        points_array = TColgp_Array1OfPnt(1, len(points))
        for i, point in enumerate(points):
            points_array.SetValue(i + 1, point)

        bspline = GeomAPI_PointsToBSpline(points_array).Curve()
        edge = BRepBuilderAPI_MakeEdge(bspline).Edge()
        wire_builder = BRepBuilderAPI_MakeWire(edge)

        return wire_builder.Wire()

    def create_flower_wire_on_yz_plane(
        self, petal_radius=20.0, num_petals=6, center_radius=5.0
    ):
        """
        YZ平面上に花形のWireを作成

        Parameters:
        petal_radius: float - 花びらの半径
        num_petals: int - 花びらの数
        center_radius: float - 中心部の半径

        Returns:
        TopoDS_Wire: 作成された花形Wire
        """
        import math

        points = []
        num_points_per_petal = 20

        for petal in range(num_petals):
            petal_angle_start = petal * 2 * math.pi / num_petals
            petal_angle_end = (petal + 1) * 2 * math.pi / num_petals

            for i in range(num_points_per_petal):
                t = i / (num_points_per_petal - 1)
                angle = petal_angle_start + t * (petal_angle_end - petal_angle_start)

                # 花びらの形状：中心から外に向かって膨らむ
                radius_variation = center_radius + (
                    petal_radius - center_radius
                ) * math.sin(math.pi * t)

                y = radius_variation * math.cos(angle)
                z = radius_variation * math.sin(angle)
                points.append(gp_Pnt(0.0, y, z))

        # 最初の点に戻って閉じる
        points.append(points[0])

        # B-スプライン曲線でスムーズな花を作成
        points_array = TColgp_Array1OfPnt(1, len(points))
        for i, point in enumerate(points):
            points_array.SetValue(i + 1, point)

        bspline = GeomAPI_PointsToBSpline(points_array).Curve()
        edge = BRepBuilderAPI_MakeEdge(bspline).Edge()
        wire_builder = BRepBuilderAPI_MakeWire(edge)

        return wire_builder.Wire()

    def create_eccentric_pair(
        self, outer_shape="cylinder", inner_shape="circle", offset_y=8.0, offset_z=3.0
    ):
        """
        偏心した外側と内側のWireペアを作成

        Parameters:
        outer_shape: str - "cylinder", "rectangle", "star", "heart", "flower", "polygon", "circle", "spiral", "gear", "wave"
        inner_shape: str - "circle", "star", "heart", "polygon", "spiral", "gear", "wave"
        offset_y: float - 内側形状のY軸オフセット
        offset_z: float - 内側形状のZ軸オフセット

        Returns:
        tuple: (外側Wire, 内側Wire)
        """
        # 外側の形状を作成
        if outer_shape == "cylinder":
            outer_wire = self.create_cylinder_wire_on_yz_plane(30.0)
        elif outer_shape == "rectangle":
            outer_wire = self.create_rectangular_wire_on_yz_plane(60.0, 40.0)
        elif outer_shape == "star":
            outer_wire = self.create_star_wire_on_yz_plane(25.0, 12.0, 6)
        elif outer_shape == "heart":
            outer_wire = self.create_heart_wire_on_yz_plane(20.0)
        elif outer_shape == "flower":
            outer_wire = self.create_flower_wire_on_yz_plane(25.0, 8, 8.0)
        elif outer_shape == "polygon":
            outer_wire = self.create_polygon_wire_on_yz_plane()
        elif outer_shape == "spiral":
            outer_wire = self.create_spiral_wire_on_yz_plane(25.0, 10.0, 2.0)
        elif outer_shape == "gear":
            outer_wire = self.create_gear_wire_on_yz_plane(25.0, 20.0, 16)
        elif outer_shape == "wave":
            outer_wire = self.create_wave_wire_on_yz_plane(8.0, 6.0, 20.0)
        else:  # circle
            outer_wire = self.create_circle_wire_on_yz_plane(0.0, 0.0, 25.0)

        # 内側の形状を作成（偏心）
        if inner_shape == "star":
            inner_wire = self.create_star_wire_on_yz_plane(8.0, 4.0, 5)
        elif inner_shape == "heart":
            inner_wire = self.create_heart_wire_on_yz_plane(8.0)
        elif inner_shape == "polygon":
            # 小さな多角形を作成
            small_points = [
                gp_Pnt(0.0, -5.0, -3.0),
                gp_Pnt(0.0, 5.0, -3.0),
                gp_Pnt(0.0, 7.0, 0.0),
                gp_Pnt(0.0, 5.0, 3.0),
                gp_Pnt(0.0, -5.0, 3.0),
                gp_Pnt(0.0, -7.0, 0.0),
            ]

            wire_builder = BRepBuilderAPI_MakeWire()
            for i in range(len(small_points)):
                start_point = small_points[i]
                end_point = small_points[(i + 1) % len(small_points)]
                edge = BRepBuilderAPI_MakeEdge(start_point, end_point).Edge()
                wire_builder.Add(edge)

            inner_wire = wire_builder.Wire()
        elif inner_shape == "spiral":
            inner_wire = self.create_spiral_wire_on_yz_plane(10.0, 3.0, 1.5)
        elif inner_shape == "gear":
            inner_wire = self.create_gear_wire_on_yz_plane(8.0, 6.0, 8)
        elif inner_shape == "wave":
            inner_wire = self.create_wave_wire_on_yz_plane(3.0, 4.0, 7.0)
        else:  # circle
            inner_wire = self.create_circle_wire_on_yz_plane(0.0, 0.0, 8.0)

        # 内側形状を偏心移動
        transform = gp_Trsf()
        transform.SetTranslation(gp_Vec(0.0, offset_y, offset_z))
        inner_wire_eccentric = BRepBuilderAPI_Transform(inner_wire, transform).Shape()

        return outer_wire, inner_wire_eccentric

    def create_spiral_wire_on_yz_plane(
        self, outer_radius=25.0, inner_radius=5.0, turns=3.0
    ):
        """
        YZ平面上に渦巻き形のWireを作成

        Parameters:
        outer_radius: float - 外側半径
        inner_radius: float - 内側半径
        turns: float - 回転数

        Returns:
        TopoDS_Wire: 作成された渦巻きWire
        """
        import math

        points = []
        num_points = 200

        for i in range(num_points):
            t = turns * 2 * math.pi * i / (num_points - 1)

            # 半径を徐々に変化させる
            radius = inner_radius + (outer_radius - inner_radius) * (
                1 - i / (num_points - 1)
            )

            y = radius * math.cos(t)
            z = radius * math.sin(t)
            points.append(gp_Pnt(0.0, y, z))

        # B-スプライン曲線で滑らかな渦巻きを作成
        points_array = TColgp_Array1OfPnt(1, len(points))
        for i, point in enumerate(points):
            points_array.SetValue(i + 1, point)

        bspline = GeomAPI_PointsToBSpline(points_array).Curve()
        edge = BRepBuilderAPI_MakeEdge(bspline).Edge()
        wire_builder = BRepBuilderAPI_MakeWire(edge)

        return wire_builder.Wire()

    def create_gear_wire_on_yz_plane(
        self, outer_radius=25.0, inner_radius=20.0, num_teeth=12
    ):
        """
        YZ平面上に歯車形のWireを作成

        Parameters:
        outer_radius: float - 外側半径（歯の先端）
        inner_radius: float - 内側半径（歯の根元）
        num_teeth: int - 歯の数

        Returns:
        TopoDS_Wire: 作成された歯車Wire
        """
        import math

        points = []

        for i in range(num_teeth):
            # 各歯に対して4つの点を作成（根元->先端->先端->根元）
            base_angle = 2 * math.pi * i / num_teeth
            tooth_width = 2 * math.pi / num_teeth / 4  # 歯の幅の半分

            # 歯の根元（左側）
            angle1 = base_angle - tooth_width
            y1 = inner_radius * math.cos(angle1)
            z1 = inner_radius * math.sin(angle1)
            points.append(gp_Pnt(0.0, y1, z1))

            # 歯の先端（左側）
            angle2 = base_angle - tooth_width / 2
            y2 = outer_radius * math.cos(angle2)
            z2 = outer_radius * math.sin(angle2)
            points.append(gp_Pnt(0.0, y2, z2))

            # 歯の先端（右側）
            angle3 = base_angle + tooth_width / 2
            y3 = outer_radius * math.cos(angle3)
            z3 = outer_radius * math.sin(angle3)
            points.append(gp_Pnt(0.0, y3, z3))

            # 歯の根元（右側）
            angle4 = base_angle + tooth_width
            y4 = inner_radius * math.cos(angle4)
            z4 = inner_radius * math.sin(angle4)
            points.append(gp_Pnt(0.0, y4, z4))

        # 最初の点に戻って閉じる
        points.append(points[0])

        # 各線分をエッジとして接続
        wire_builder = BRepBuilderAPI_MakeWire()
        for i in range(len(points) - 1):
            edge = BRepBuilderAPI_MakeEdge(points[i], points[i + 1]).Edge()
            wire_builder.Add(edge)

        return wire_builder.Wire()

    def create_wave_wire_on_yz_plane(self, amplitude=15.0, frequency=4.0, radius=20.0):
        """
        YZ平面上に波形円のWireを作成

        Parameters:
        amplitude: float - 波の振幅
        frequency: float - 波の周波数
        radius: float - 基本半径

        Returns:
        TopoDS_Wire: 作成された波形円Wire
        """
        import math

        points = []
        num_points = 200

        for i in range(num_points + 1):
            angle = 2 * math.pi * i / num_points

            # 基本円に波形を追加
            wave_radius = radius + amplitude * math.sin(frequency * angle)

            y = wave_radius * math.cos(angle)
            z = wave_radius * math.sin(angle)
            points.append(gp_Pnt(0.0, y, z))

        # B-スプライン曲線で滑らかな波形を作成
        points_array = TColgp_Array1OfPnt(1, len(points))
        for i, point in enumerate(points):
            points_array.SetValue(i + 1, point)

        bspline = GeomAPI_PointsToBSpline(points_array).Curve()
        edge = BRepBuilderAPI_MakeEdge(bspline).Edge()
        wire_builder = BRepBuilderAPI_MakeWire(edge)

        return wire_builder.Wire()

    def create_complex_hollow_structures(self):
        """複雑で面白い中空構造を作成（実用的な外形 + 芸術的内形）"""
        print("\n=== 実用的外形 + 芸術的内形の中空構造 ===")

        structures = []

        # 1. 円柱の中に偏心した星
        print("1. 円柱 + 偏心星...")
        cylinder_outer, star_inner = self.create_eccentric_pair(
            "cylinder", "star", 5.0, 3.0
        )
        cylinder_star_hollow = self.create_hollow_solid(
            cylinder_outer, star_inner, 60.0, "loft"
        )
        structures.append(("円柱+偏心星", cylinder_star_hollow, "BLUE"))

        # 2. 四角柱の中に偏心ハート
        print("2. 四角柱 + 偏心ハート...")
        rectangle_outer, heart_inner = self.create_eccentric_pair(
            "rectangle", "heart", -4.0, 6.0
        )
        rectangle_heart_hollow = self.create_hollow_solid(
            rectangle_outer, heart_inner, 50.0, "sweep"
        )
        structures.append(("四角柱+偏心ハート", rectangle_heart_hollow, "RED"))

        # 3. 円柱の中に偏心花形
        print("3. 円柱 + 偏心花...")
        cylinder_outer2, flower_inner = self.create_eccentric_pair(
            "cylinder", "flower", 10.0, -4.0
        )
        cylinder_flower_hollow = self.create_hollow_solid(
            cylinder_outer2, flower_inner, 90.0, "loft"
        )
        structures.append(("円柱+偏心花", cylinder_flower_hollow, "GREEN"))

        # 4. 四角柱の中に偏心渦巻き
        print("4. 四角柱 + 偏心渦巻き...")
        rectangle_outer2, spiral_inner = self.create_eccentric_pair(
            "rectangle", "spiral", 12.0, 6.0
        )
        rectangle_spiral_hollow = self.create_hollow_solid(
            rectangle_outer2, spiral_inner, 85.0, "sweep"
        )
        structures.append(("四角柱+偏心渦巻き", rectangle_spiral_hollow, "YELLOW"))

        # 5. 円柱の中に偏心歯車
        print("5. 円柱 + 偏心歯車...")
        cylinder_outer3, gear_inner = self.create_eccentric_pair(
            "cylinder", "gear", -8.0, 10.0
        )
        cylinder_gear_hollow = self.create_hollow_solid(
            cylinder_outer3, gear_inner, 75.0, "loft"
        )
        structures.append(("円柱+偏心歯車", cylinder_gear_hollow, "PURPLE"))

        # 6. 四角柱の中に偏心波形
        print("6. 四角柱 + 偏心波形...")
        rectangle_outer3, wave_inner = self.create_eccentric_pair(
            "rectangle", "wave", 7.0, -9.0
        )
        rectangle_wave_hollow = self.create_hollow_solid(
            rectangle_outer3, wave_inner, 95.0, "sweep"
        )
        structures.append(("四角柱+偏心波形", rectangle_wave_hollow, "ORANGE"))

        return structures

    # ...existing code...


if __name__ == "__main__":
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument(
        "--pxyz", dest="pxyz", default=[0.0, 0.0, 0.0], type=float, nargs=3
    )
    parser.add_argument(
        "--demo",
        dest="demo",
        default="hollow",
        help="Demo mode: 'all', 'sweep', 'section', 'hollow', 'artistic', 'practical', 'morphing', 'area' (default: hollow)",
    )
    opt = parser.parse_args()
    print(opt)

    # YZPlaneSweepクラスのインスタンスを作成
    obj = YZPlaneSweep()

    # Wireの情報を表示
    print("\n=== Wire情報 ===")
    test_circle = create_circle_wire_on_yz_plane(radius=20.0)
    test_polygon = create_polygon_wire_on_yz_plane()
    test_curve = create_smooth_curve_wire_on_yz_plane()

    print_wire_info(test_circle, "円形Wire")
    print_wire_info(test_polygon, "多角形Wire")
    print_wire_info(test_curve, "スムーズ曲線Wire")

    # デモモードに応じて実行
    if opt.demo == "all":
        # 基本デモを実行
        obj.run_sweep_demo()
        obj.create_single_cross_section_demo()
        obj.run_hollow_structure_demo()
        obj.run_artistic_shapes_demo()
        
        # 新しい関数を使用したデモも実行
        print("\n=== 新しい関数を使用した追加デモ ===")
        try:
            # 基本形状デモ
            basic_shapes = demonstrate_basic_shapes()
            for i, shape in enumerate(basic_shapes[:2]):  # 最初の2つのみ表示
                obj.display.DisplayShape(shape, transparency=0.2)
            
            # 複雑形状デモ
            complex_shapes = demonstrate_complex_shapes()
            for i, shape in enumerate(complex_shapes[:2]):  # 最初の2つのみ表示
                obj.display.DisplayShape(shape, transparency=0.3)
            
            print("追加デモ完了")
        except Exception as e:
            print(f"追加デモでエラー: {e}")
    elif opt.demo == "sweep":
        # Sweepデモのみ
        obj.run_sweep_demo()
    elif opt.demo == "section":
        # 断面デモのみ
        obj.create_single_cross_section_demo()
        obj.show_axs_pln()
        obj.ShowOCC()
    elif opt.demo == "hollow":
        # 中空構造デモのみ
        obj.run_hollow_structure_demo()
        
        # 新しい中空構造デモも追加
        print("\n=== 追加中空構造デモ ===")
        try:
            hollow_shapes = demonstrate_hollow_structures()
            for shape in hollow_shapes[:1]:  # 最初の1つのみ表示
                obj.display.DisplayShape(shape, transparency=0.4)
        except Exception as e:
            print(f"追加中空構造デモでエラー: {e}")
        
        obj.show_axs_pln()
        obj.ShowOCC()
    elif opt.demo == "artistic":
        # 芸術的形状デモのみ
        obj.run_artistic_shapes_demo()
        
        # 新しい複雑形状デモも追加
        print("\n=== 追加芸術形状デモ ===")
        try:
            complex_shapes = demonstrate_complex_shapes()
            for shape in complex_shapes[:3]:  # 最初の3つのみ表示
                obj.display.DisplayShape(shape, transparency=0.3)
        except Exception as e:
            print(f"追加芸術形状デモでエラー: {e}")
        
        obj.show_axs_pln()
        obj.ShowOCC()
    elif opt.demo == "practical":
        # 実用的外形デモ
        print("=== 実用的外形デモ ===")
        try:
            # 歯車や波形など実用的な形状を表示
            gear = create_gear_wire_on_yz_plane(25.0, 20.0, 12)
            gear_solid = sweep_wire_along_x_axis(gear, 50.0)
            obj.display.DisplayShape(gear_solid, transparency=0.2)
            
            wave = create_wave_wire_on_yz_plane(8.0, 6.0, 20.0)
            wave_solid = create_lofted_solid_from_wire(wave, 80.0)
            obj.display.DisplayShape(wave_solid, transparency=0.3)
            
            obj.show_axs_pln()
            obj.ShowOCC()
        except Exception as e:
            print(f"practical デモでエラー: {e}")
    elif opt.demo == "morphing":
        # モーフィングデモ
        print("=== モーフィングデモ ===")
        try:
            morphing_shapes = demonstrate_morphing_shapes()
            for shape in morphing_shapes:
                obj.display.DisplayShape(shape, transparency=0.3)
            obj.show_axs_pln()
            obj.ShowOCC()
        except Exception as e:
            print(f"morphing デモでエラー: {e}")
    elif opt.demo == "area":
        # 断面面積解析デモ
        print("=== 断面面積解析デモ ===")
        try:
            # テストSolidを作成
            star = create_star_wire_on_yz_plane(20.0, 10.0, 6)
            star_solid = sweep_wire_along_x_axis(star, 100.0)
            
            # 断面解析
            sections = analyze_solid_cross_sections(star_solid, 5)
            obj.display.DisplayShape(star_solid, transparency=0.4)
            
            # 断面を表示
            for x_pos, section in sections:
                obj.display.DisplayShape(section, transparency=0.1)
            
            obj.show_axs_pln()
            obj.ShowOCC()
        except Exception as e:
            print(f"area デモでエラー: {e}")
    elif opt.demo == "visualization":
        # 切断面可視化デモ
        print("=== 切断面可視化デモ ===")
        try:
            # 複雑形状を作成
            flower = create_flower_wire_on_yz_plane(25.0, 8, 8.0)
            flower_solid = create_lofted_solid_from_wire(flower, 120.0)
            obj.display.DisplayShape(flower_solid, transparency=0.5)
            
            # 複数の断面を作成して可視化
            for i in range(0, 121, 20):
                try:
                    section = get_cross_section_at_x(flower_solid, i)
                    obj.display.DisplayShape(section, transparency=0.2)
                except:
                    pass
                    
            obj.show_axs_pln()
            obj.ShowOCC()
        except Exception as e:
            print(f"visualization デモでエラー: {e}")
    else:
        print(f"未知のデモモード: {opt.demo}")
        print(
            "利用可能なモード: 'all', 'sweep', 'section', 'hollow', 'artistic'"
        )
