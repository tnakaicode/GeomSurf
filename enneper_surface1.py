import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import sys
import os
import time
import argparse

sys.path.append(os.path.join("../"))
from base import plot2d, plot3d

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

    obj = plot3d(aspect="equal")

    u = np.linspace(-1, 1, 100) * np.pi
    v = np.linspace(-1, 1, 100) * np.pi
    u, v = np.meshgrid(u, v)
    u, v = u.flatten(), v.flatten()

    x = u * (1 - u ** 2 / 3 + v ** 2) / 3
    y = -v * (1 - v ** 2 / 3 + u ** 2) / 3
    z = (u ** 2 - v ** 2) / 3

    tri = mtri.Triangulation(u, v)

    obj.axs.plot_trisurf(x, y, z, triangles=tri.triangles, cmap="jet")
    obj.set_axes_equal("xy")
    obj.axs.set_title('$Enneper Surface$')
    obj.SavePng()
