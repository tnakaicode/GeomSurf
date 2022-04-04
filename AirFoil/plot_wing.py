from OCC.Core.gp import gp_Pnt
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import argparse
from linecache import getline, clearcache

sys.path.append(os.path.join("../"))
from src.base_occ import dispocc

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCCUtils.Construct import make_polygon

if __name__ == '__main__':
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument("--pxyz", dest="pxyz",
                      default=[0.0, 0.0, 0.0], type=float, nargs=3)
    opt = parser.parse_args()
    print(opt, argvs)

    obj = dispocc()

    wingfile = "wing001.txt"
    wing_cfg = {
        'names': ('airfoil', 'dx', 'size', 'px', 'py', 'deg'),
        'formats': ('U15', float, float, float, float, float)
    }
    data = np.loadtxt(wingfile, dtype=wing_cfg, skiprows=2)
    for i, dat in enumerate(data):
        dat_file = "./uiuc_dat/" + dat[0] + ".dat"
        dx, s, px, py, d = dat[1], dat[2], dat[3], dat[4], dat[5]
        print(dat_file)
        xy = np.loadtxt(dat_file, skiprows=1)
        pts = []
        for xy_2d in xy:
            x, y = xy_2d[0] - px, xy_2d[1] - py
            pnt = gp_Pnt(x * s, y * s, dx)
            pts.append(pnt)
        air = make_polygon(pts, closed=True)
        obj.display.DisplayShape(air)
    obj.show_axs_pln()
    obj.ShowOCC()
