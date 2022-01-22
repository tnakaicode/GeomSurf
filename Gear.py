import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
from optparse import OptionParser

sys.path.append(os.path.join("./"))
from src.base import plot2d, plotocc

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

if __name__ == '__main__':
    argvs = sys.argv
    parser = OptionParser()
    parser.add_option("--dir", dest="dir", default="./")
    parser.add_option("--pxyz", dest="pxyz",
                      default=[0.0, 0.0, 0.0], type="float", nargs=3)
    opt, argc = parser.parse_args(argvs)
    print(opt, argc)

    px = np.linspace(-1, 1, 100) * 100 + 50
    py = np.linspace(-1, 1, 200) * 100 - 50
    mesh = np.meshgrid(px, py)

    obj = dispocc()
    r0 = 10.0
    r1 = 11.0
    r2 = (r0 + r1) / 2
    num = 50
    dig = 51
    sft = 0.75
    pts = []
    for idx in range(num):
        d0 = 2 * np.pi * idx / num
        d1 = 2 * np.pi * (idx + 1) / num
        for t in [0.0, sft - 0.01, sft, 1.0]:
            rad = d0 + t * (d1 - d0)
            if t < sft:
                radi = r0
            else:
                radi = r1
            x = radi * np.cos(rad)
            y = radi * np.sin(rad)
            z = 0
            pts.append(gp_Pnt(x, y, z))
            print(t, d0, d1, rad, radi)
    pts.append(pts[0])

    plate = obj.make_plate(pts, skin=0.5)
    print(plate)
    obj.display.DisplayShape(plate)
    # for pnt in pts:
    #    obj.display.DisplayShape(pnt)
    obj.show_occ()
