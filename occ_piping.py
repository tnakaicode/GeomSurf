import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import argparse
from linecache import getline, clearcache, updatecache

basename = os.path.dirname(__file__)

sys.path.append(os.path.join("./"))
from src.base_occ import dispocc

import logging

logging.getLogger("matplotlib").setLevel(logging.ERROR)

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Circ
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Compound
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCCUtils.Common import vertex2pnt
from OCCUtils.Construct import make_box, make_line, make_wire, make_edge, make_face
from OCCUtils.Construct import make_plane, make_polygon
from OCCUtils.Construct import point_to_vector, vector_to_point
from OCCUtils.Construct import dir_to_vec, vec_to_dir

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
    axs = gp_Ax3()
    axs_pln = obj.make_plane_axs(axs, [-100, 100], [-100, 200])

    a0 = gp_Ax3(gp_Pnt(-50, 50, 0), gp_Dir(0, 0, -1), gp_Dir(-1, 0, 0))
    r0 = 5.0
    c0 = make_wire(make_edge(gp_Circ(a0.Ax2(), r0)))
    a1 = gp_Ax3(gp_Pnt(50, 150, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0))
    r1 = 10.0
    c1 = make_wire(make_edge(gp_Circ(a1.Ax2(), r1)))

    axs_pln = make_face(axs_pln, c0)
    # axs_pln = make_face(axs_pln, c1)

    obj.show_axs_pln()
    obj.show_axs_pln(a0, scale=r0)
    obj.show_axs_pln(a1, scale=r1)
    obj.display.DisplayShape(axs_pln)
    obj.ShowOCC()
