import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import argparse

sys.path.append(os.path.join("./"))
from src.base_occ import dispocc, spl_face

import logging

logging.getLogger("matplotlib").setLevel(logging.ERROR)

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Compound
from OCC.Core.TopLoc import TopLoc_Location
from OCCUtils.Topology import Topo
from OCCUtils.Construct import make_box, make_line, make_wire, make_edge
from OCCUtils.Construct import make_plane, make_polygon
from OCCUtils.Construct import point_to_vector, vector_to_point
from OCCUtils.Construct import dir_to_vec, vec_to_dir


def curvature(px, r, s):
    """( x + sx )**2 / 2*rx + ( y + sy )**2 / 2*ry"""
    if r == 0:
        py = np.zeros_like(px + s)
    else:
        py = (px + s) ** 2 / (2 * r)
    return py


if __name__ == "__main__":
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument("--rxy", dest="rxy", default=[0.0, 0.0], type=float, nargs=2)
    parser.add_argument("--lxy", dest="lxy", default=[200, 200], type=float, nargs=2)
    parser.add_argument("--nxy", dest="nxy", default=[500, 500], type=int, nargs=2)
    opt = parser.parse_args()
    print(opt, argvs)

    obj = dispocc(disp=False)
    if opt.dir != None:
        obj.tmpdir = "./stp_surf/"
    px = np.linspace(-1, 1, opt.nxy[0]) * opt.lxy[0] / 2
    py = np.linspace(-1, 1, opt.nxy[1]) * opt.lxy[1] / 2
    mesh = np.meshgrid(px, py)
    surf_x = curvature(mesh[0], opt.rxy[0], 0)
    surf_y = curvature(mesh[1], opt.rxy[1], 0)
    surf = surf_x + surf_y
    face = spl_face(*mesh, surf)
    obj.export_stp(face)
