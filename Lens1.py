import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import argparse

sys.path.append(os.path.join("./"))
from src.base_occ import dispocc, set_loc
from src.MakeLens import make_lens, make_lens1

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.BOPAlgo import BOPAlgo_Builder, BOPAlgo_Splitter
from OCC.Core.BRepFeat import BRepFeat_MakeDPrism
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakePrism
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Compound
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Core.TopAbs import TopAbs_EDGE, TopAbs_FACE, TopAbs_SHAPE, TopAbs_SHELL, TopAbs_SOLID
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Extend.DataExchange import read_step_file, write_step_file
from OCCUtils.Topology import Topo
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
    print(opt, argvs)

    obj = dispocc(disp=False)

    face1 = read_step_file("./stp_surf/gen_surf_030.stp")
    axis1 = gp_Ax3(gp_Pnt(0, 0, 0.0), gp_Dir(0, 0, 1))
    face1.Location(set_loc(gp_Ax3(), axis1))

    face2 = read_step_file("./stp_surf/gen_surf_030.stp")
    axis2 = gp_Ax3(gp_Pnt(0, 0, 3.5), gp_Dir(0, 0, 1))
    face2.Location(set_loc(gp_Ax3(), axis2))

    face3 = make_plane()
    axis3 = gp_Ax3(gp_Pnt(0, 0, 5.5), gp_Dir(0, 0, 1))
    face3.Location(set_loc(gp_Ax3(), axis3))

    axs = gp_Ax3()
    ax1 = gp_Ax3(gp_Pnt(0, 0, -10), axs.Direction())
    vec = gp_Vec(gp_Pnt(0, 0, -10), gp_Pnt(0, 0, 10))
    face = obj.make_EllipWire(rxy=[50.0, 50.0], axs=ax1, shft=0)
    body = BRepPrimAPI_MakePrism(face, vec).Shape()

    lens = make_lens1(body, face1, face2, face3)
    obj.create_tempdir(flag=1)
    for shp in lens:
        obj.export_stp(shp)
