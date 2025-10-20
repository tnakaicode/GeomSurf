from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge
from OCC.Core.BRep import BRep_Tool
from OCC.Core.Geom import Geom_BezierCurve
from OCC.Core.TColgp import TColgp_Array1OfPnt
from OCC.Core.gp import gp_Pnt, gp_Vec
from OCC.Display.SimpleGui import init_display

# ビューアーを初期化
display, start_display, add_menu, add_function_to_menu = init_display()

# 制御点を定義
points = [
    gp_Pnt(0, 0, 0),
    gp_Pnt(1, 2, 0),
    gp_Pnt(4, 0, 0),
]

# 制御点を表示
for p in points:
    display.DisplayShape(p, update=True)

# 始点・終点の接線ベクトルを定義
start_vec = gp_Vec(1, 0.9, 0.1).Normalized()
end_vec = gp_Vec(1, -1, 0).Normalized()

# 始点と終点の接線ベクトルに基づいて新しい制御点を計算
control_point1 = points[0].Translated(start_vec.Scaled(2.0))
control_point2 = points[-1].Translated(end_vec.Reversed().Scaled(2.0))

# 全ての制御点を使用して配列を作成
all_points = [points[0], control_point1, control_point2, points[-1]]
points_array = TColgp_Array1OfPnt(1, len(all_points))
for i, point in enumerate(all_points):
    points_array.SetValue(i + 1, point)

# ベジェ曲線を作成
curve = Geom_BezierCurve(points_array)

# ベジェ曲線の形状を作成
curve_shape = BRepBuilderAPI_MakeEdge(curve).Edge()

# 始点と終点の接線ベクトルを表示
display.DisplayVector(start_vec, points[0], update=True)
display.DisplayVector(end_vec, points[-1], update=True)

# 曲線の始点・終点の接線ベクトルを計算
first_parameter = curve.FirstParameter()
last_parameter = curve.LastParameter()
tangent_start = gp_Vec()
tangent_end = gp_Vec()

curve.D1(first_parameter, points[0], tangent_start)
curve.D1(last_parameter, points[-1], tangent_end)

# 接線ベクトルを正規化
tangent_start.Normalize()
tangent_end.Normalize()

# ベクトルの成分を表示
print(
    "Calculated start tangent vector: ({}, {}, {})".format(
        tangent_start.X(), tangent_start.Y(), tangent_start.Z()
    )
)
print(
    "Given start tangent vector: ({}, {}, {})".format(
        start_vec.X(), start_vec.Y(), start_vec.Z()
    )
)
print(
    "Calculated end tangent vector: ({}, {}, {})".format(
        tangent_end.X(), tangent_end.Y(), tangent_end.Z()
    )
)
print(
    "Given end tangent vector: ({}, {}, {})".format(
        end_vec.X(), end_vec.Y(), end_vec.Z()
    )
)

# 始点・終点の接線ベクトルが一致しているかの判定
start_vec.Normalize()
end_vec.Normalize()

# if tangent_start.IsEqual(start_vec, 1e-6) and tangent_end.IsEqual(end_vec, 1e-6):
#    print("The tangents match the given vectors.")
# else:
#    print("The tangents do not match the given vectors.")

# ベジェ曲線を表示
display.DisplayShape(curve_shape, update=True)
display.FitAll()
start_display()


# ウェイト（Weight）：
#
# ウェイトは、各制御点が曲線に与える影響の大きさを調整します。ウェイトが大きいほど、曲線はその制御点に強く引かれるようになります。
#
# 接線ベクトル（Tangent Vector）：
#
# ベジェ曲線の接線ベクトルは、曲線の始点や終点での方向を示します。これは、曲線の滑らかさや形状に影響を与えます。
#
# 具体的には、始点や終点での接線ベクトルは、隣接する制御点の位置とその重みに依存します。次のような関係があります：
#
# 始点での接線ベクトル：
#
# 始点での接線ベクトルは、最初の制御点とその次の制御点の位置、およびそれらのウェイトによって決まります。
#
# 接線ベクトルは、通常、最初の制御点から次の制御点へのベクトルに比例します。ウェイトが大きいほど、接線ベクトルはその制御点に引かれるようになります。
#
# 終点での接線ベクトル：
#
# 終点での接線ベクトルは、最後の制御点とその前の制御点の位置、およびそれらのウェイトによって決まります。
#
# 接線ベクトルは、通常、最後の制御点からその前の制御点へのベクトルに比例します。ウェイトが大きいほど、接線ベクトルはその制御点に引かれるようになります。
#
