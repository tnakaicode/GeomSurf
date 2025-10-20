import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import argparse
from linecache import getline, clearcache, updatecache

basename = os.path.dirname(__file__)

sys.path.append(os.path.join("./"))
from base_occ import dispocc

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepAdaptor import BRepAdaptor_Curve
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Compound
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Extend.TopologyUtils import TopologyExplorer, ordered_vertices_from_wire, ordered_edges_from_wire
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

    box1 = make_box(10, 10, 10)
    print(hash(box1))

    box2 = make_box(10, 10, 10)
    print(hash(box2))

    face = [f for f in TopologyExplorer(box2).faces()][0]
    wire = [w for w in TopologyExplorer(face).wires()][0]
    print([BRep_Tool.Pnt(v) for v in TopologyExplorer(wire).vertices()])
    print([BRep_Tool.Pnt(v) for v in ordered_vertices_from_wire(wire)])
    for e in ordered_edges_from_wire(wire):
        c, t0, t1 = BRep_Tool.Curve(e)
        print(c.Value(t0), c.Value(t1))
