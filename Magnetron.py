import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import argparse

sys.path.append(os.path.join("./"))
from src.base_occ import dispocc, set_trf

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Compound
from OCC.Core.TopLoc import TopLoc_Location
from OCCUtils.Topology import Topo
from OCCUtils.Construct import make_box, make_line, make_wire, make_edge
from OCCUtils.Construct import make_plane, make_polygon
from OCCUtils.Construct import point_to_vector, vector_to_point
from OCCUtils.Construct import dir_to_vec, vec_to_dir


def pnt_radi(radi=10.0, rad=0.0, axs=gp_Ax3()):
    x = radi * np.cos(rad)
    y = radi * np.sin(rad)
    z = 0
    pnt = gp_Pnt(x, y, z)
    pnt.Transform(set_trf(gp_Ax3(), axs))
    return pnt


if __name__ == '__main__':
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument("--pxyz", dest="pxyz",
                      default=[0.0, 0.0, 0.0], type="float", nargs=3)
    opt = parser.parse_args()
    print(opt, argvs)

    px = np.linspace(-1, 1, 100) * 100 + 50
    py = np.linspace(-1, 1, 200) * 100 - 50
    mesh = np.meshgrid(px, py)

    obj = dispocc(touch=True)
    r0 = 25.0
    r1 = 10.0
    r2 = 2.5
    wth = 0.1
    sft = 0.5
    num = 6
    pts = []
    for idx in range(num):
        d0 = 2 * np.pi * idx / num
        d1 = 2 * np.pi * (idx + 1) / num
        dg = (d1 - d0)
        t0 = sft - wth / 2
        t1 = sft + wth / 2
        pt1 = pnt_radi(r1, d0 + dg * sft)
        vc1 = gp_Vec(gp_Pnt(), pt1)
        print(pt1)
        ax1 = gp_Ax3()
        ax1.SetLocation(pt1)
        ax1.SetXDirection(vec_to_dir(vc1))
        pts += [pnt_radi(radi=r0, rad=d0 + dg * t)
                for t in np.linspace(0, t0, 10)]
        pts.append(pnt_radi(radi=r1 + r2, rad=d0 + dg * t0))
        pts += [pnt_radi(radi=r2, rad=-2 * np.pi * t, axs=ax1)
                for t in np.linspace(0.25, 0.75, 15)]
        pts.append(pnt_radi(radi=r1 + r2, rad=d0 + dg * t1))
        pts += [pnt_radi(radi=r0, rad=d0 + dg * t)
                for t in np.linspace(t1, 1, 10)]
    pts.append(pts[0])

    plate = obj.make_plate(pts, skin=1.0)
    print(plate)
    obj.display.DisplayShape(plate)
    # for pnt in pts:
    #    obj.display.DisplayShape(pnt)
    obj.ShowOCC()
