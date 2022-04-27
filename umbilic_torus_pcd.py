import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import argparse

from OCC.Core.gp import gp_Pnt, gp_PntMirror, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3

sys.path.append(os.path.join("../"))
from base_occ import dispocc

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

if __name__ == '__main__':
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument("--pxyz", dest="pxyz",
                      default=[0.0, 0.0, 0.0], type=float, nargs=3)
    opt = parser.parse_args()
    print(opt, argvs)

    u = np.linspace(-1, 1, 100) * np.pi
    v = np.linspace(-1, 1, 100) * np.pi
    u, v = np.meshgrid(u, v)

    x = np.sin(u) * (7 + np.cos(u / 3 - 2 * v) + 2 * np.cos(u / 3 + v))
    y = np.cos(u) * (7 + np.cos(u / 3 - 2 * v) + 2 * np.cos(u / 3 + v))
    z = np.sin(u / 3 - 2 * v) + 2 * np.sin(u / 3 + v)

    obj = dispocc(touch=True)
    for idx, val in np.ndenumerate(x):
        x1, y1, z1 = x[idx], y[idx], z[idx]
        pnt = gp_Pnt(x1, y1, z1)
        obj.display.DisplayShape(pnt)
    obj.ShowOCC()
