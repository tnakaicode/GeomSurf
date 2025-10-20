import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors
from OCC.Core.BRepPrimAPI import (
    BRepPrimAPI_MakeBox,
    BRepPrimAPI_MakeCylinder,
    BRepPrimAPI_MakeSphere,
)
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut
from OCC.Core.BRepIntCurveSurface import BRepIntCurveSurface_Inter
from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRepBndLib import brepbndlib
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE
from OCC.Core.gp import gp_Pnt, gp_Dir, gp_Lin, gp_Ax2
from OCC.Core.Quantity import Quantity_Color, Quantity_TOC_RGB
from OCC.Display.SimpleGui import init_display
import matplotlib.pyplot as plt


# バウンディングボックスの8点を取得し、全方向に大きくする
def get_bounding_box_corners(shape, margin=1.0):
    bbox = Bnd_Box()
    brepbndlib.Add(shape, bbox)
    xmin, ymin, zmin, xmax, ymax, zmax = bbox.Get()
    corners = [
        gp_Pnt(xmin - margin, ymin - margin, zmin - margin),
        gp_Pnt(xmin - margin, ymin - margin, zmax + margin),
        gp_Pnt(xmin - margin, ymax + margin, zmin - margin),
        gp_Pnt(xmin - margin, ymax + margin, zmax + margin),
        gp_Pnt(xmax + margin, ymin - margin, zmin - margin),
        gp_Pnt(xmax + margin, ymin - margin, zmax + margin),
        gp_Pnt(xmax + margin, ymax + margin, zmin - margin),
        gp_Pnt(xmax + margin, ymax + margin, zmax + margin),
    ]
    return corners


# 表面上に点群を生成
def generate_surface_point_cloud(shape, num_points, unit_sphere_centers):
    points = []
    explorer = TopExp_Explorer(shape, TopAbs_FACE)

    for center in unit_sphere_centers:
        explorer.ReInit()
        while explorer.More():
            face = explorer.Current()
            face_intersector = BRepIntCurveSurface_Inter()

            count = 0
            trials = 0
            while count < num_points and trials < num_points * 2:
                trials += 1
                # 単位球上に点を生成
                theta = np.random.uniform(0, np.pi)
                phi = np.random.uniform(0, 2 * np.pi)
                x = np.sin(theta) * np.cos(phi)
                y = np.sin(theta) * np.sin(phi)
                z = np.cos(theta)
                unit_sphere_point = gp_Pnt(
                    center.X() + x, center.Y() + y, center.Z() + z
                )
                direction = gp_Dir(x, y, z)
                line = gp_Lin(unit_sphere_point, direction)

                # 点と面の交差を判定
                face_intersector.Init(face, line, 1e-6)
                if face_intersector.More():
                    pnt = face_intersector.Pnt()
                    points.append([pnt.X(), pnt.Y(), pnt.Z()])
                    count += 1

            explorer.Next()

    return np.array(points)


# 法線ベクトルを計算
def calculate_normals(points, k=30):
    nn = NearestNeighbors(n_neighbors=k).fit(points)
    neighbors = nn.kneighbors(points, return_distance=False)
    normals = []
    for indices in neighbors:
        p = points[indices]
        p_mean = p.mean(axis=0)
        cov = np.cov((p - p_mean).T)
        eigenvalues, eigenvectors = np.linalg.eigh(cov)
        normal = eigenvectors[:, 0]
        normals.append(normal)
    return np.array(normals)


# 法線ベクトルの変化を計算
def calculate_normal_changes(normals, k=30):
    nn = NearestNeighbors(n_neighbors=k).fit(normals)
    neighbors = nn.kneighbors(normals, return_distance=False)
    normal_changes = []
    for i, indices in enumerate(neighbors):
        n = normals[indices]
        n_mean = n.mean(axis=0)
        change = np.linalg.norm(normals[i] - n_mean)
        normal_changes.append(change)
    return np.array(normal_changes)


# 法線ベクトルの変化に基づいてクラスタ分け
def cluster_points_by_normal_changes(points, normal_changes, eps=0.05, min_samples=10):
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(
        normal_changes.reshape(-1, 1))
    labels = db.labels_

    clusters = []
    unique_labels = set(labels)
    print(
        f"Number of clusters: {len(unique_labels) - (1 if -1 in unique_labels else 0)}"
    )

    for label in unique_labels:
        if label == -1:
            continue  # ノイズとして無視
        cluster_points = points[labels == label]
        clusters.append(cluster_points)

    return clusters


# クラスタごとの色を生成
def generate_cluster_colors(num_clusters):
    colors = plt.cm.get_cmap("hsv", num_clusters)
    return [colors(i) for i in range(num_clusters)]


if __name__ == "__main__":
    # メイン処理
    display, start_display, add_menu, add_function_to_menu = init_display()

    # 形状を作成
    box = BRepPrimAPI_MakeBox(10, 10, 10).Shape()
    cylinder = BRepPrimAPI_MakeCylinder(5, 8).Shape()
    sphere = BRepPrimAPI_MakeSphere(5).Shape()

    # ブール演算で形状を結合・減算
    combined_shape = BRepAlgoAPI_Fuse(box, cylinder).Shape()
    final_shape = BRepAlgoAPI_Cut(combined_shape, sphere).Shape()

    # 細い円柱で穴をあける
    thin_cylinder = BRepPrimAPI_MakeCylinder(
        gp_Ax2(gp_Pnt(2, 0, 0), gp_Dir(0, 0.1, 1)), 1, 15
    ).Shape()
    final_shape_with_hole = BRepAlgoAPI_Cut(final_shape, thin_cylinder).Shape()

    hole_sphere = BRepPrimAPI_MakeSphere(gp_Pnt(5, 5, 5), 2).Shape()
    final_shape_with_hole = BRepAlgoAPI_Cut(
        final_shape_with_hole, hole_sphere).Shape()

    display.DisplayShape(final_shape_with_hole, transparency=0.9)

    # バウンディングボックスの8点から単位球を作成
    unit_sphere_centers = get_bounding_box_corners(
        final_shape_with_hole, margin=1.5)
    for center in unit_sphere_centers:
        unit_sphere = BRepPrimAPI_MakeSphere(center, 1).Shape()
        # display.DisplayShape(unit_sphere, transparency=0.5)

    # 表面上の点群を生成
    num_points = 1000
    points = generate_surface_point_cloud(
        final_shape_with_hole, num_points, unit_sphere_centers
    )

    # 法線ベクトルを計算
    normals = calculate_normals(points)

    # 法線ベクトルの変化を計算
    normal_changes = calculate_normal_changes(normals)

    # 法線ベクトルの変化に基づいて点群をクラスタ分け
    clusters = cluster_points_by_normal_changes(points, normal_changes)

    # クラスタごとの色を生成
    colors = generate_cluster_colors(len(clusters))

    # 各クラスタを表示
    for cluster, color in zip(clusters, colors):
        quantity_color = Quantity_Color(*color[:3], Quantity_TOC_RGB)
        for point in cluster:
            display.DisplayShape(gp_Pnt(*point), color=quantity_color)

    display.FitAll()
    start_display()
