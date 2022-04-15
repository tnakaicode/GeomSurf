import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import argparse
from linecache import getline, clearcache

basename = os.path.dirname(__file__)

sys.path.append(os.path.join("./"))
from src.base_occ import dispocc, spl_face
from src.geometry import curvature
from src.Coord import OCCSurfObj

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
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument("--pxyz", dest="pxyz",
                      default=[0.0, 0.0, 0.0], type=float, nargs=3)
    opt = parser.parse_args()
    print(opt, argvs)

    obj = dispocc(touch=True)

    beam0 = OCCSurfObj(name="beam0")
    beam0.beam = gp_Ax3(beam0.axs.Location(), gp_Dir(0.1, 0.1, 1.0))

    surf1 = OCCSurfObj(name="surf1")
    surf1.axs = gp_Ax3(gp_Pnt(0, 0, 50), gp_Dir(0, 1, 1))
    surf1.rot = gp_Ax3(surf1.axs.Ax2())
    surf1.SurfCurvature_Loc(lxy=[100, 100], rxy=[0, 0])

    surf1.reflect_beam(beam0.beam, tr=1)
    norm1 = obj.reflect_beam(surf1.surf, beam0.beam, tr=3)

    surf2 = OCCSurfObj(name="surf2")
    surf2.axs = gp_Ax3(gp_Pnt(0, 0, 50), gp_Dir(0, 1, 1))
    surf2.rot = gp_Ax3(surf2.axs.Ax2())
    surf2.SurfCurvature_Loc(lxy=[100, 100], rxy=[500, 0])
    norm2 = obj.reflect_beam(surf2.surf, beam0.beam, tr=4)

    obj.show_axs_pln(beam0.beam, scale=75)
    obj.show_axs_pln(surf1.beam, scale=25)
    obj.show_vec(norm1, scale=15)
    obj.show_vec(norm2, scale=15)
    obj.display.DisplayShape(surf1.surf, transparency=0.9)
    obj.display.DisplayShape(surf2.surf, transparency=0.9, color="RED")
    obj.display.View_Right()
    obj.ShowOCC()

    # Eye:  136.9416885375977 19.63805389404297 39.830164432525635
    # Pos:  12.313896179199276 19.63805389404297 39.830164432525635
    # Prj:  1.0 -0.0 -0.0
    # Up :  0.0 0.0 1.0
    #
    # Eye:  136.9416885375977 9.974650374071416 50.03858881008746
    # Pos:  12.313896179199276 9.974650374071416 50.03858881008746
    # Prj:  1.0 -0.0 -0.0
    # Up :  0.0 0.0 1.0
