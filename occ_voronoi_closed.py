import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d

# 乱数生成のための種を設定
np.random.seed(42)

# 200~300のランダムな点を生成
num_points = np.random.randint(200, 301)
points = np.random.rand(num_points, 2)

# Voronoi図の計算
vor = Voronoi(points)


# 無限領域を閉じるための外角形の設定
def voronoi_finite_polygons_2d(vor, radius=None):
    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    if radius is None:
        radius = np.ptp(vor.points, axis=0).max()

    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]

        if all(v >= 0 for v in vertices):
            new_regions.append(vertices)
            continue

        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1

            if v1 >= 0:
                continue

            t = vor.points[p2] - vor.points[p1]
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            far_point = (
                vor.vertices[v2] + np.sign(np.dot(midpoint - center, n)) * n * radius
            )

            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        new_regions.append(new_region)

    return new_regions, np.asarray(new_vertices)


# 無限領域を閉じたVoronoi図の計算
regions, vertices = voronoi_finite_polygons_2d(vor)

# 結果のプロット
fig, ax = plt.subplots()
for region in regions:
    polygon = vertices[region]
    ax.fill(*zip(*polygon), edgecolor="black", alpha=0.5)

ax.plot(points[:, 0], points[:, 1], "o")
ax.plot(vertices[:, 0], vertices[:, 1], "o")
# ax.set_xlim(points[:, 0].min() - 0.1, points[:, 0].max() + 0.1)
# ax.set_ylim(points[:, 1].min() - 0.1, points[:, 1].max() + 0.1)
ax.set_title("Voronoi Diagram with Finite Polygons")
plt.show()
