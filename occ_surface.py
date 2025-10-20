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

from geomdl import fitting, BSpline, exchange
from geomdl.visualization import VisMPL

from OCC.Core.STEPControl import STEPControl_Writer, STEPControl_AsIs
from OCC.Core.IFSelect import IFSelect_RetDone
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Common
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.Core.gp import gp_Pnt
from OCC.Core.TopoDS import TopoDS_Shape


def generate_point_cloud(center, size, num_points):
    # 点群の生成（例としてランダムな点群を生成）
    return center + size * (np.random.rand(num_points, 3) - 0.5)


def create_nurbs_surface_from_point_cloud(point_cloud):
    # 点群をグリッド状に再サンプリング（余分な点を無視）
    grid_size = int(np.sqrt(len(point_cloud)))
    if grid_size < 2:
        raise ValueError("点群の数が少なすぎて曲面を生成できません。")
    points = point_cloud[:grid_size**2].reshape((grid_size, grid_size, 3))

    # NURBS曲面をフィッティング
    surf = fitting.approximate_surface(
        points, size_u=grid_size, size_v=grid_size, degree_u=3, degree_v=3)
    return surf


def visualize_surface(surface):
    # 曲面を可視化
    surface.vis = VisMPL.VisSurface()
    surface.render()


def extend_surface(surface, extension=2.0):
    # 曲面を拡張して交差する曲線を計算
    for ctrlpt in surface.ctrlpts:
        ctrlpt[0] += extension
        ctrlpt[1] += extension
        ctrlpt[2] += extension
    surface.evaluate()
    return surface


def compute_intersection_curve(surface1, surface2):
    # 交差する曲線を計算（OCCを使用）
    face1 = BRepBuilderAPI_MakeFace(surface1).Face()
    face2 = BRepBuilderAPI_MakeFace(surface2).Face()
    intersection = BRepAlgoAPI_Common(face1, face2)
    return intersection.Shape()


def export_to_step(shape, filename="output.step"):
    # STEPファイルに出力
    step_writer = STEPControl_Writer()
    step_writer.Transfer(shape, STEPControl_AsIs)
    status = step_writer.Write(filename)
    if status == IFSelect_RetDone:
        print(f"STEPファイルが正常に出力されました: {filename}")
    else:
        print("STEPファイルの出力に失敗しました。")


if __name__ == "__main__":

    # 点群1と点群2を生成
    point_cloud1 = generate_point_cloud(center=np.array([0, 0, 0]),
                                        size=10,
                                        num_points=100)

    point_cloud2 = generate_point_cloud(center=np.array([5, 5, 5]),
                                        size=10,
                                        num_points=100)

    # 点群1と点群2から曲面を生成
    nurbs_surface1 = create_nurbs_surface_from_point_cloud(point_cloud1)
    nurbs_surface2 = create_nurbs_surface_from_point_cloud(point_cloud2)

    # visualize_surface(nurbs_surface1)
    # visualize_surface(nurbs_surface2)

    extended_surface1 = extend_surface(nurbs_surface1)
    extended_surface2 = extend_surface(nurbs_surface2)

    intersection_curve = compute_intersection_curve(extended_surface1,
                                                    extended_surface2)

    export_to_step(intersection_curve)
