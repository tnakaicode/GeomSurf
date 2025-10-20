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

import logging

logging.getLogger("matplotlib").setLevel(logging.ERROR)

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Compound
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCCUtils.Common import vertex2pnt
from OCCUtils.Construct import make_box, make_line, make_wire, make_edge
from OCCUtils.Construct import make_plane, make_polygon
from OCCUtils.Construct import point_to_vector, vector_to_point
from OCCUtils.Construct import dir_to_vec, vec_to_dir


def generate_point_cloud(num_points, mean_xy, cov_xy, mean_z, std_z):
    """
    点群を生成する関数。

    Parameters:
        num_points (int): 点の数
        mean_xy (list): XY平面上の2D正規分布の平均
        cov_xy (list): XY平面上の2D正規分布の共分散行列
        mean_z (float): Z方向の正規分布の平均
        std_z (float): Z方向の正規分布の標準偏差

    Returns:
        np.ndarray: 生成された点群 (num_points x 3)
    """
    # XY平面上の2D正規分布に従う点を生成
    xy_points = np.random.multivariate_normal(mean_xy, cov_xy, num_points)

    # Z方向に正規分布する値を生成
    z_points = np.random.normal(mean_z, std_z, num_points)

    # 点群を結合
    points = np.hstack((xy_points, z_points.reshape(-1, 1)))
    return points


if __name__ == "__main__":
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument(
        "--pxyz", dest="pxyz", default=[0.0, 0.0, 0.0], type=float, nargs=3
    )
    opt = parser.parse_args()
    print(opt)

    obj = dispocc(touch=True)
    axs = gp_Ax3()

    # パラメータ設定
    num_points = 1000  # 点の数
    mean_xy = [0.0, 0.0]  # XY平面上の2D正規分布の平均
    cov_xy = [[10.0, 0.2], [0.2, 5.0]]  # XY平面上の2D正規分布の共分散行列
    mean_z = 0.0  # Z方向の正規分布の平均
    std_z = 0.0  # Z方向の正規分布の標準偏差

    # 点群を生成
    points = generate_point_cloud(num_points, mean_xy, cov_xy, mean_z, std_z)

    for i in range(num_points):
        pnt = gp_Pnt(points[i, 0], points[i, 1], points[i, 2])
        obj.display.DisplayShape(pnt)

    obj.show_axs_pln(scale=10)
    obj.ShowOCC()
