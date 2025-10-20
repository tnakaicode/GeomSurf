import numpy as np
import open3d as o3d
from sklearn.cluster import KMeans
from sklearn.cluster import DBSCAN
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeVertex, BRepBuilderAPI_MakeFace
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
from OCC.Core.BRepTools import breptools
from OCC.Core.TopAbs import TopAbs_FACE
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.gp import gp_Pnt, gp_Pln, gp_Dir
from OCC.Display.SimpleGui import init_display


# 表面上に点群を生成
def generate_surface_point_cloud(shape, num_points):
    points = []
    explorer = TopExp_Explorer(shape, TopAbs_FACE)
    while explorer.More():
        face = explorer.Current()
        adaptor = BRepAdaptor_Surface(face)
        u_min, u_max, v_min, v_max = breptools.UVBounds(face)
        for _ in range(num_points // 6):  # キューブの6面に分けて点を生成
            u = np.random.uniform(u_min, u_max)
            v = np.random.uniform(v_min, v_max)
            point = adaptor.Value(u, v)
            points.append([point.X(), point.Y(), point.Z()])
        explorer.Next()
    return np.array(points)


# クラスタ分け
def cluster_point_cloud(points, n_clusters):
    kmeans = KMeans(n_clusters=n_clusters)
    kmeans.fit(points)
    return kmeans.labels_


# クラスタ分け
def cluster_point_cloud1(points):
    # 法線ベクトルを使用せず、点群データのみを使用して同一面上にある点をクラスタリングする方法として、DBSCAN (Density-Based Spatial Clustering of Applications with Noise) を用いることができます。これは、密度に基づいたクラスタリングアルゴリズムで、形状の面を抽出するのに適しています。
    dbscan = DBSCAN(eps=1.0, min_samples=10)  # パラメータは必要に応じて調整してください
    labels = dbscan.fit_predict(points)
    unique_labels = np.unique(labels)
    clusters = [points[labels == label] for label in unique_labels if label != -1]
    return labels, clusters


# クラスタ分け（RANSACを使用）
def cluster_point_cloud_o3d(
    points, distance_threshold=0.1, ransac_n=3, num_iterations=1000, min_inliers=100
):
    point_cloud = o3d.geometry.PointCloud()
    point_cloud.points = o3d.utility.Vector3dVector(points)

    # 複数の平面を検出
    planes = []
    plane_labels = np.full(len(points), -1)  # 初期値は-1（クラスタ外の点）
    label_counter = 0
    remaining_points = point_cloud

    while len(remaining_points.points) >= ransac_n:  # RANSACの最小点数を考慮
        plane_model, inliers = remaining_points.segment_plane(
            distance_threshold=distance_threshold,
            ransac_n=ransac_n,
            num_iterations=num_iterations,
        )
        if len(inliers) < min_inliers:  # 面を構成する点が少ない場合は終了
            break
        planes.append((plane_model, inliers))
        plane_labels[inliers] = label_counter
        label_counter += 1
        remaining_points = remaining_points.select_by_index(inliers, invert=True)

    return (
        plane_labels,
        planes,
        [points[plane_labels == i] for i in range(label_counter)],
    )


# クラスタ分け
def cluster_point_cloud_by_surface(points, normals):
    clustering = DBSCAN(eps=0.1, min_samples=10, metric="cosine").fit(normals)
    return clustering.labels_, [
        points[clustering.labels_ == i]
        for i in np.unique(clustering.labels_)
        if i != -1
    ]


# 3D描画
def draw_point_cloud(display, points, labels):
    colors = ["RED", "GREEN", "BLUE1", "YELLOW", "MAGENTA", "BLACK"]  # クラスタごとの色
    for i, point in enumerate(points):
        vertex = BRepBuilderAPI_MakeVertex(gp_Pnt(*point)).Vertex()
        display.DisplayShape(vertex, update=True, color=colors[labels[i] % len(colors)])


# 平面を描画
def create_face_from_points(points):
    """点群から外形を作成し、平面を生成"""
    from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakePolygon

    polygon = BRepBuilderAPI_MakePolygon()
    for p in points:
        polygon.Add(gp_Pnt(*p))
    polygon.Close()
    return BRepBuilderAPI_MakeFace(polygon.Wire()).Face()


# メイン処理
display, start_display, add_menu, add_function_to_menu = init_display()

num_points = 1000
box = BRepPrimAPI_MakeBox(10, 10, 10).Shape()
points = generate_surface_point_cloud(box, num_points)
labels, planes, clusters = cluster_point_cloud_o3d(
    points, distance_threshold=0.1, ransac_n=3, num_iterations=1000, min_inliers=50
)

# 点群と平面をクラスタごとに描画
colors = ["RED", "GREEN", "BLUE1", "YELLOW", "CYAN1", "BLACK"]  # クラスタごとの色
for idx, (plane_model, inliers) in enumerate(planes):
    # 点を描画
    for i in inliers:
        vertex = BRepBuilderAPI_MakeVertex(gp_Pnt(*points[i])).Vertex()
        display.DisplayShape(vertex, update=True, color=colors[idx % len(colors)])

    # 外形を作成し平面を描画
    cluster_points = points[inliers[:3]]
    face = create_face_from_points(cluster_points)
    display.DisplayShape(face, update=True, color=colors[idx % len(colors)])

display.FitAll()
start_display()
