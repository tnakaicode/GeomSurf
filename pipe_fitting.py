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
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_MakePipe, BRepOffsetAPI_MakePipeShell
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Compound
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCCUtils.Common import vertex2pnt
from OCCUtils.Construct import make_box, make_line, make_wire, make_edge, make_face
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

    circ = make_edge(gp_Circ(gp_Ax2(gp_Pnt(0, 50, 0), gp_Dir(1, 0, 0)), 50),
                     np.pi / 2, np.pi / 2 + np.pi / 2)
    edg1 = make_edge(gp_Pnt(), gp_Pnt(0, 0, 100))

    wire = make_wire(edg1, circ)
    face = make_face(make_wire(make_edge(gp_Circ(gp_Ax2(), 50))))
    cylinder = BRepOffsetAPI_MakePipe(wire, face).Shape()

    obj.display.DisplayShape(cylinder, transparency=0.9)
    obj.display.DisplayShape(circ)
    obj.display.DisplayShape(edg1)
    # obj.show_axs_pln()
    obj.ShowOCC()
