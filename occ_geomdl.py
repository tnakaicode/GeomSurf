import numpy as np

from geomdl import fitting, BSpline
from geomdl.visualization import VisMPL
from geomdl import exchange

from OCC.Core.gp import gp_Pnt
from OCC.Core.TColgp import TColgp_Array2OfPnt
from OCC.Core.GeomAPI import GeomAPI_PointsToBSplineSurface
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.Display.SimpleGui import init_display


def generate_point_cloud(center, size, num_points):
    # 点群の生成（例としてランダムな点群を生成）
    return center + size * (np.random.rand(num_points, 3) - 0.5)


def create_nurbs_surface_from_point_cloud(point_cloud, size_u, size_v):
    # 点群からNURBS曲面を生成（非整列な点群に対応）
    # 点群をリスト形式に変換
    points = point_cloud.tolist()

    # NURBS曲面を補間
    surf = fitting.interpolate_surface(
        points, size_u=size_u, size_v=size_v, degree_u=3, degree_v=3
    )
    return surf


def debug_nurbs_surface(surface):
    print("=== NURBS Surface Debug Info ===")
    print(f"Degree U: {surface.degree_u}")
    print(f"Degree V: {surface.degree_v}")
    print(f"Number of Control Points (U): {len(surface.ctrlpts2d)}")
    print(f"Number of Control Points (V): {len(surface.ctrlpts2d[0])}")
    print(f"Knot Vector U: {surface.knotvector_u}")
    print(f"Knot Vector V: {surface.knotvector_v}")
    print(f"Control Points (First Row): {surface.ctrlpts2d[0]}")
    print("================================")


def convert_to_occ_surface(surface):
    # geomdl の NURBS 曲面を pythonocc-core の形式に変換
    # geomdl の制御点を取得
    ctrlpts = surface.ctrlpts2d
    size_u = len(ctrlpts)
    size_v = len(ctrlpts[0])

    # 制御点を OCC の形式に変換
    pnt_2d = TColgp_Array2OfPnt(1, size_u, 1, size_v)
    for i in range(size_u):
        for j in range(size_v):
            x, y, z = ctrlpts[i][j]
            pnt_2d.SetValue(i + 1, j + 1, gp_Pnt(x, y, z))

    # B-Spline 曲面を生成
    api = GeomAPI_PointsToBSplineSurface(pnt_2d)
    occ_surface = BRepBuilderAPI_MakeFace(api.Surface(), 1e-6).Face()
    return occ_surface


# メイン処理
if __name__ == "__main__":
    # グリッドサイズを指定（例: 10x10）
    size_u = 10
    size_v = 10

    # 規則的な点群を生成
    x = np.linspace(0, 10, size_u)
    y = np.linspace(0, 10, size_v)
    xv, yv = np.meshgrid(x, y)
    zv = np.sin(xv) + np.cos(yv)  # Z値を計算
    point_cloud = np.dstack([xv, yv, zv]).reshape(-1, 3)
    point_cloud = generate_point_cloud(np.array([0, 0, 0]), 0.5, size_u * size_v)

    # 点群をリスト形式に変換
    points = point_cloud.tolist()

    # NURBS曲面を補間
    nurbs_surface = fitting.interpolate_surface(
        points, size_u=size_u, size_v=size_v, degree_u=3, degree_v=3
    )
    exchange.export_stl(nurbs_surface, "occ_geomdl.stl")

    display, start_display, add_menu, add_function_to_menu = init_display()

    # NURBS 曲面を pythonocc-core の形式に変換
    occ_surface = convert_to_occ_surface(nurbs_surface)

    # 描画
    display.DisplayShape(occ_surface, update=True)
    start_display()
