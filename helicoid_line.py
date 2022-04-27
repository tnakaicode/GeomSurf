import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import argparse

sys.path.append(os.path.join("../"))
from base_occ import dispocc as plotocc

from OCC.Core.gp import gp_Pnt, gp_PntMirror, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_ThruSections, BRepOffsetAPI_MakeOffset, BRepOffsetAPI_MakeEvolved, BRepOffsetAPI_MakePipe, BRepOffsetAPI_MakePipeShell
from OCCUtils.Construct import make_polygon, make_circle, make_vertex, make_edge, make_wire
from OCCUtils.Construct import point_to_vector, vector_to_point
from OCCUtils.Construct import dir_to_vec, vec_to_dir

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

    obj = plotocc(touch=True)
    r0 = 25.0
    r1 = 100.0
    c = 50.0

    u0 = np.linspace(r0, r1, 200)
    v0 = np.linspace(-1, 1, 200) * 2 * np.pi

    api = BRepOffsetAPI_ThruSections()
    api.SetSmoothing(True)
    for vt in v0:
        p0_x, p0_y, p0_z = r0 * np.cos(vt), r0 * np.sin(vt), c * vt
        p1_x, p1_y, p1_z = r1 * np.cos(vt), r1 * np.sin(vt), c * vt
        p0 = gp_Pnt(p0_x, p0_y, p0_z)
        p1 = gp_Pnt(p1_x, p1_y, p1_z)
        edge = make_edge(p0, p1)
        api.AddWire(make_wire(edge))
        obj.display.DisplayShape(edge)
    api.Build()
    obj.display.DisplayShape(api.Shape())

    obj.show_axs_pln()
    obj.ShowOCC()
