import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
from sklearn.neighbors import KernelDensity


# 点群の生成関数
def generate_points(num_points, center, radius, x_range, y_range):
    points = []
    while len(points) < num_points:
        point = np.random.rand(2)
        point[0] = x_range[0] + point[0] * (x_range[1] - x_range[0])
        point[1] = y_range[0] + point[1] * (y_range[1] - y_range[0])
        if (
            np.linalg.norm(point - center) > radius
        ):  # 中心から半径radiusの円の外に点を配置
            points.append(point)
    return np.array(points)


# 点群の生成
num_points = 1000
center = np.array([0.5, 0.5])
radius = 0.3
x_range = (0, 1)
y_range = (0, 1)
points = generate_points(num_points, center, radius, x_range, y_range)

# 点群密度の計算
kde = KernelDensity(bandwidth=0.05).fit(points)
x, y = np.meshgrid(
    np.linspace(x_range[0], x_range[1], 100), np.linspace(y_range[0], y_range[1], 100)
)
xy = np.vstack([x.ravel(), y.ravel()]).T
density = np.exp(kde.score_samples(xy)).reshape(x.shape)

# Delaunay分割の実行
tri = Delaunay(points)


# 密度の薄い領域を交差するエッジを削除する関数
def remove_low_density_edges(tri, points, kde, density_threshold):
    mask = []
    for simplex in tri.simplices:
        coords = points[simplex]
        # エッジが通る中点の密度を計算
        edge_midpoints = [
            (coords[i] + coords[j]) / 2
            for i in range(len(coords))
            for j in range(i + 1, len(coords))
        ]
        densities = np.exp(kde.score_samples(edge_midpoints))
        if all(density >= density_threshold for density in densities):
            mask.append(True)
        else:
            mask.append(False)
    return mask


# 密度閾値の設定
density_threshold = 0.1

# 密度の薄い領域を交差するエッジを削除
mask = remove_low_density_edges(tri, points, kde, density_threshold)

# 分割結果と点群密度のプロット
plt.imshow(
    density,
    origin="lower",
    extent=(x_range[0], x_range[1], y_range[0], y_range[1]),
    cmap="Blues",
    alpha=0.5,
)
plt.colorbar(label="Density")
plt.triplot(points[:, 0], points[:, 1], tri.simplices[mask], "k-")
plt.plot(points[:, 0], points[:, 1], "o")
plt.gca().set_aspect("equal")
plt.show()
