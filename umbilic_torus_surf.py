import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import argparse

from OCC.Core.gp import gp_Pnt, gp_PntMirror, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_ThruSections
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCCUtils.Construct import make_polygon

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

    u1 = np.linspace(-1, 1, 100) * np.pi
    v1 = np.linspace(-1, 1, 100) * np.pi
    u, v = np.meshgrid(u1, v1)

    x = np.sin(u) * (7 + np.cos(u / 3 - 2 * v) + 2 * np.cos(u / 3 + v))
    y = np.cos(u) * (7 + np.cos(u / 3 - 2 * v) + 2 * np.cos(u / 3 + v))
    z = np.sin(u / 3 - 2 * v) + 2 * np.sin(u / 3 + v)

    obj = dispocc(touch=True)
    api = BRepOffsetAPI_ThruSections()
    for i0, v0 in enumerate(v1):
        pts = []
        for j0, u0 in enumerate(u1):
            pts.append(gp_Pnt(x[i0, j0], y[i0, j0], z[i0, j0]))
        wire = make_polygon(pts)
        api.AddWire(wire)
    api.SetContinuity(1)
    api.Build()
    shp = api.Shape()
    print(list(TopologyExplorer(shp).faces()),
          TopologyExplorer(shp).number_of_faces())
    obj.display.DisplayShape(api.Shape(), transparency=0.9)
    obj.ShowOCC()
