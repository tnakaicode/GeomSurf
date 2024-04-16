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
logging.getLogger('matplotlib').setLevel(logging.ERROR)

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Circ
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_MakePipe
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Compound
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCCUtils.Common import vertex2pnt
from OCCUtils.Construct import make_box, make_line, make_wire, make_edge
from OCCUtils.Construct import make_plane, make_polygon
from OCCUtils.Construct import point_to_vector, vector_to_point
from OCCUtils.Construct import dir_to_vec, vec_to_dir

if __name__ == '__main__':
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument("--pxyz", dest="pxyz",
                        default=[0.0, 0.0, 0.0], type=float, nargs=3)
    opt = parser.parse_args()
    print(opt)

    obj = dispocc(touch=True)
    axs = gp_Ax3()

    major = 2000
    minor = 1000
    radii = major + minor

    axs_major = gp_Ax2(gp_Pnt(0, 0, 0.0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0))
    axs_minor = gp_Ax2(gp_Pnt(0, major, 0.0), gp_Dir(1, 0, 0), gp_Dir(0, 1, 0))
    tok = gp_Ax3(gp_Pnt(0, 0, 0.0), gp_Dir(0, 1, 0), gp_Dir(-1, 0, 0))

    cir_major = make_edge(gp_Circ(axs_major, major), -np.pi / 2, np.pi / 2)
    cir_minor = make_edge(gp_Circ(axs_minor, minor), -np.pi / 2, np.pi / 2)
    api = BRepOffsetAPI_MakePipe(make_wire(cir_major), cir_minor)
    api.Build()
    shll = api.Shape()
    face = [f for f in TopologyExplorer(shll).faces()][0]
    surf = BRepAdaptor_Surface(face)

    obj.show_axs_pln()
    obj.display.DisplayShape(shll)
    obj.display.DisplayShape(cir_minor)
    obj.display.DisplayShape(cir_major)
    obj.ShowOCC()
