import matplotlib.pyplot as plt
import numpy as np
from OCC.Core.GeomAPI import GeomAPI_PointsToBSpline
from OCC.Core.Geom import Geom_BezierCurve
from OCC.Core.TColgp import TColgp_Array1OfPnt
from OCC.Core.gp import gp_Pnt

# 制御点を定義
points = [gp_Pnt(0, 0, 0), gp_Pnt(1, 2, 0), gp_Pnt(4, 0, 0)]


# Bezier曲線を生成
def create_bezier_curve(points, weights):
    array = TColgp_Array1OfPnt(1, len(points))
    for i, point in enumerate(points):
        array.SetValue(i + 1, point)

    bezier = Geom_BezierCurve(array)
    for i, weight in enumerate(weights):
        bezier.SetPole(i + 1, points[i], weight)
    return bezier


# 接線ベクトルを計算
def calculate_tangent(bezier, u):
    tangent_vec = bezier.DN(u, 1)
    return np.array([tangent_vec.X(), tangent_vec.Y(), tangent_vec.Z()])


# 始点と終点のWeightを変更してBezier曲線を生成
weights_list = [
    # [0.1, 1, 0.1],
    [1, 1, 1],
    [2, 1, 2],
    [3, 1, 3],
]

# グラフ化
plt.figure(figsize=(10, 6))
for weights in weights_list:
    bezier = create_bezier_curve(points, weights)
    u_values = np.linspace(bezier.FirstParameter(), bezier.LastParameter(), 100)
    tangents = np.array([calculate_tangent(bezier, u) for u in u_values])
    plt.plot(tangents[:, 0], tangents[:, 1], label=f"Weights: {weights}")

plt.title("Bezier Curve Tangent Vectors")
plt.xlabel("u")
plt.ylabel("Tangent Vector X")
plt.legend()
plt.show()
