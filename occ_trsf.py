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
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3, gp_Trsf
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Compound
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCCUtils.Common import vertex2pnt
from OCCUtils.Construct import make_box, make_line, make_wire, make_edge
from OCCUtils.Construct import make_plane, make_polygon
from OCCUtils.Construct import point_to_vector, vector_to_point
from OCCUtils.Construct import dir_to_vec, vec_to_dir


def set_axs_trf(ax1=gp_Ax3(), ax2=gp_Ax3()):
    trf = gp_Trsf()
    if ax1 == None:
        trf.SetTransformation(ax2)
    else:
        trf.SetTransformation(ax2, ax1)
    return trf


def set_axs_loc(ax1=gp_Ax3(), ax2=gp_Ax3()):
    trf = set_axs_trf(ax1, ax2)
    loc = TopLoc_Location(trf)
    return loc


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
    axs = gp_Ax3(gp_Pnt(10, 10, 10), gp_Dir(0.1, 0.2, 1.0))

    pts = [
        gp_Pnt(-50, -50, 0),
        gp_Pnt(50, -50, 0),
        gp_Pnt(50, 50, 0),
        gp_Pnt(-50, 50, 0),
    ]
    poly = make_polygon(pts, closed=True)
    poly.Move(set_axs_loc(gp_Ax3(), axs))

    ax1 = gp_Ax3(gp_Pnt(100, 100, 0), gp_Dir(0, 0, 1))
    ax2 = gp_Ax3(gp_Pnt(175, 160, 0), gp_Dir(0, 0.1, 1), gp_Dir(0.5, 0.0, 0.0))
    poly.Move(set_axs_loc(axs, ax1))

    trf = gp_Trsf()
    trf.SetRotation(axs.Axis(), np.deg2rad(30))
    poly1 = poly.Moved(TopLoc_Location(trf))
    ax1_1 = ax1.Transformed(trf)

    trf = gp_Trsf()
    trf.SetRotation(ax1.Axis(), np.deg2rad(30))
    poly2 = poly.Moved(TopLoc_Location(trf))
    ax1_2 = ax1.Transformed(trf)

    poly3 = poly
    for i in range(5):
        trf = gp_Trsf()
        trf.SetRotation(ax2.Axis(), np.deg2rad(30))
        poly3 = poly3.Moved(TopLoc_Location(trf))
        obj.display.DisplayShape(poly3, color="YELLOW")
        # ax2.Transform(trf)

    obj.show_axs_pln()
    obj.show_axs_pln(ax1_1, scale=50)
    obj.show_axs_pln(ax1_2, scale=50)
    obj.display.DisplayShape(poly)
    obj.display.DisplayShape(poly1, color="BLUE1")
    obj.display.DisplayShape(poly2, color="GREEN")
    obj.ShowOCC()
