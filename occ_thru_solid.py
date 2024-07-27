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
from OCC.Core.gp import gp_Circ, gp_Lin
from OCC.Core.BRepProj import BRepProj_Projection
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeSphere
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_ThruSections
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

    sph = BRepPrimAPI_MakeSphere(axs.Ax2(), 100).Shape()
    lin = gp_Lin(gp_Pnt(0, 0, 0), gp_Dir(0, 0.5, 1))
    c1 = gp_Circ(axs.Ax2(), 50.0)
    w1 = make_wire(make_edge(c1))

    proj = BRepProj_Projection(w1, sph, lin.Direction())
    proj_wire = []
    while proj.More():
        proj_wire.append(proj.Current())
        proj.Next()

    api = BRepOffsetAPI_ThruSections(True, False, 1.0E-6)
    api.AddWire(proj_wire[0])
    api.AddWire(proj_wire[1])
    api.Build()
    sol = api.Shape()

    obj.display.DisplayShape(sph, transparency=0.9)
    obj.display.DisplayShape(sol, transparency=0.9, color="BLUE1")
    obj.display.DisplayShape(w1)
    obj.display.DisplayShape(make_edge(lin, -150, 150))
    obj.show_axs_pln()
    obj.ShowOCC()
