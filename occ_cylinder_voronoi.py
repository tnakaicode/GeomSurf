import numpy as np
from scipy.spatial import Voronoi
from OCC.Core.gp import gp_Pnt, gp_Ax2, gp_Dir, gp_Ax3
from OCC.Core.TopAbs import TopAbs_FACE
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.Geom import Geom_Surface
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
from OCC.Core.BRepLProp import BRepLProp_SLProps
from OCC.Core.GeomLProp import GeomLProp_SLProps
from OCC.Core.GCE2d import GCE2d_MakeSegment
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge
from OCC.Core.gp import gp_Pnt2d
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeCone
from OCCUtils.Construct import make_polygon
from base_occ import dispocc


def make_cylinder(axis, radius, height, taper_angle=0.0):
    """
    円筒またはテーパー付き円筒を生成し、表面 (Face) を返す関数
    :param axis: 円筒の軸 (gp_Ax2)
    :param radius: 円筒の底面半径
    :param height: 円筒の高さ
    :param taper_angle: テーパー角度 (ラジアン, デフォルトは0で通常の円筒)
    :return: TopoDS_Face (円筒の表面)
    """
    top_radius = radius - height * \
        np.tan(taper_angle) if taper_angle != 0.0 else radius
    if top_radius <= 0:
        raise ValueError("テーパー角度が大きすぎて上面半径が負になります。")

    # 円筒の形状を作成
    cone = BRepPrimAPI_MakeCone(axis, radius, top_radius, height).Shape()

    # 円筒の表面 (Face) を取得
    explorer = TopExp_Explorer(cone, TopAbs_FACE)
    face = None
    while explorer.More():
        face = explorer.Current()
        break  # 最初のFaceを取得 (外側の表面)

    if face is None:
        raise ValueError("円筒の表面を取得できませんでした。")

    return face


def generate_uv_points_random(num_points):
    """
    UVパラメータ空間上でランダムな点群を生成
    :param num_points: 点の数
    :return: UVパラメータの点群 (numpy.ndarray)
    """
    u = np.random.uniform(0, 2 * np.pi, num_points)  # U方向 (0 ～ 2π)
    v = np.random.uniform(0, 1, num_points)         # V方向 (0 ～ 1, 正規化された高さ)
    return np.column_stack((u, v))


def map_uv_to_cylinder(uv_points, cylinder, height):
    """
    UVパラメータ空間の点を円筒表面の3D座標にマッピング
    :param uv_points: UVパラメータの点群 (numpy.ndarray)
    :param cylinder: Taper Cylinder (TopoDS_Face)
    :param height: 円筒の高さ
    :return: 3D座標の点群 (list of gp_Pnt)
    """
    points = []
    for u, v in uv_points:
        point = uv_to_point_on_cylinder(cylinder, u, v, height)
        points.append(point)
    return points


def voronoi_on_cylinder(points, radius, height):
    """
    円筒表面上の点群をVoronoi分割 (U方向の周期性を考慮)
    :param points: 円筒表面上の3D点群 (numpy.ndarray)
    :param radius: 円筒の底面半径
    :param height: 円筒の高さ
    :return: Voronoi分割結果 (scipy.spatial.Voronoi)
    """
    # 点群を2D平面に投影 (角度と高さ)
    angles = np.arctan2(points[:, 1], points[:, 0])  # U方向 (角度)
    heights = points[:, 2]  # V方向 (高さ)
    projected_points = np.column_stack((angles, heights))

    # U方向の周期性を考慮するために、点群を複製
    extended_points = np.vstack([
        projected_points,
        np.column_stack(
            (projected_points[:, 0] + 2 * np.pi, projected_points[:, 1])),
        np.column_stack(
            (projected_points[:, 0] - 2 * np.pi, projected_points[:, 1]))
    ])

    # Voronoi分割を計算
    vor = Voronoi(extended_points)

    # 元の点群に対応する領域のみを返す
    valid_regions = []
    for region_idx in vor.point_region[:len(projected_points)]:
        region = vor.regions[region_idx]
        if not region or -1 in region:
            continue  # 無効な領域をスキップ
        valid_regions.append(region)

    # 元の点群に対応するVoronoiオブジェクトを再構築
    vor.filtered_regions = valid_regions
    vor.filtered_vertices = vor.vertices

    return vor


def voronoi_on_uv(uv_points):
    """
    UVパラメータ空間上でVoronoi分割を実行 (U方向の周期性を考慮)
    :param uv_points: UVパラメータ空間上の点群 (numpy.ndarray)
    :return: Voronoi分割結果 (scipy.spatial.Voronoi)
    """
    # U方向の周期性を考慮するために、点群を複製
    extended_points = np.vstack([
        uv_points,
        np.column_stack((uv_points[:, 0] + 2 * np.pi, uv_points[:, 1])),
        np.column_stack((uv_points[:, 0] - 2 * np.pi, uv_points[:, 1]))
    ])

    # Voronoi分割を計算
    vor = Voronoi(extended_points)

    # 元の点群に対応する領域のみを返す
    valid_regions = []
    for region_idx in vor.point_region[:len(uv_points)]:
        region = vor.regions[region_idx]
        if not region or -1 in region:
            continue  # 無効な領域をスキップ
        valid_regions.append(region)

    # 元の点群に対応するVoronoiオブジェクトを再構築
    vor.filtered_regions = valid_regions
    vor.filtered_vertices = vor.vertices

    return vor


def create_curve_from_uv(uv_start, uv_end, cylinder):
    """
    UVパラメータの2点を元にTaper Cylinder上の曲線を構成
    :param uv_start: UVパラメータの始点 (u0, v0)
    :param uv_end: UVパラメータの終点 (u1, v1)
    :param cylinder: Taper Cylinder (TopoDS_Shape)
    :return: 曲線 (TopoDS_Edge)
    """
    # 円筒のサーフェスを取得
    surface_adaptor = BRepAdaptor_Surface(cylinder)
    surface = surface_adaptor.Surface()

    # UVパラメータを2D線分として定義
    u0, v0 = uv_start
    u1, v1 = uv_end
    line_2d = GCE2d_MakeSegment(gp_Pnt2d(u0, v0), gp_Pnt2d(u1, v1)).Value()

    # 曲線を生成
    edge = BRepBuilderAPI_MakeEdge(line_2d, surface).Edge()
    return edge


def uv_to_point_on_cylinder(cylinder, u, v, height):
    """
    UVパラメータを元にTaper Cylinder上の3D座標を取得
    :param cylinder: Taper Cylinder (TopoDS_Face)
    :param u: Uパラメータ (角度)
    :param v: Vパラメータ (高さの正規化, 0 <= v <= 1)
    :param height: 円筒の高さ
    :return: 円筒表面上の3D座標 (gp_Pnt)
    """
    # 円筒のサーフェスを取得
    surface_adaptor = BRepAdaptor_Surface(cylinder)

    # Vパラメータを高さにスケール
    scaled_v = v * height

    # UVパラメータを元に3D座標を取得
    props = BRepLProp_SLProps(surface_adaptor, u, scaled_v, 2, 1.0e-9)
    return props.Value()


if __name__ == "__main__":
    # 描画の初期化
    obj = dispocc(touch=True)
    axs = gp_Ax3()

    # 円筒のパラメータ
    radius = 50.0
    height = 100.0
    num_points = 200  # 点の数
    taper_angle = np.radians(5)  # テーパー角度 (5度)

    # UVパラメータ空間でランダムな点群を生成
    uv_points = generate_uv_points_random(num_points)

    # Voronoi分割をUV空間で計算 (U方向の周期性を考慮)
    vor = voronoi_on_uv(uv_points)

    # 基準となる円筒を描画 (テーパー付き)
    cylinder = make_cylinder(axs.Ax2(), radius, height, taper_angle)
    obj.display.DisplayShape(cylinder, transparency=0.8, update=False)

    # UVパラメータを3D座標にマッピング
    points = map_uv_to_cylinder(uv_points, cylinder, height)

    # 点群を描画
    for point in points:
        obj.display.DisplayShape(point, update=False)

    # Voronoi分割結果を円筒表面に描画
    for region in vor.filtered_regions:
        vertices = []
        for idx in region:
            u, v = vor.filtered_vertices[idx]
            angle = u % (2 * np.pi)  # U方向の周期性を戻す
            z = v * height  # V方向を高さに変換

            # 半径を計算 (テーパーを考慮)
            r = radius - z * \
                np.tan(taper_angle) if taper_angle != 0.0 else radius

            # 円筒範囲外の頂点を除外
            if z < 0 or z > height or r <= 0:
                continue

            x = r * np.cos(angle)
            y = r * np.sin(angle)
            vertices.append(gp_Pnt(x, y, z))

        # 領域をポリゴンとして描画
        if len(vertices) > 2:
            polygon = make_polygon(vertices, closed=True)
            obj.display.DisplayShape(polygon)

    # 座標軸を表示
    obj.show_axs_pln()

    # 描画ループを開始
    obj.ShowOCC()
