import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import argparse

sys.path.append(os.path.join("./"))
from src.base_occ import dispocc, spl_face

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

    obj = dispocc(disp=False)
    if opt.dir != None:
        obj.tmpdir = "./stp_surf/"
    px = np.linspace(-1, 1, 200) * 150
    py = np.linspace(-1, 1, 200) * 150
    mesh = np.meshgrid(px, py)
    for rx in np.linspace(-1, 1, 6) * 1000:
        for ry in np.linspace(-1, 1, 6) * 1000:
            surf = mesh[0]**2 / rx + mesh[1]**2 / ry
            face = spl_face(*mesh, surf)
            obj.export_stp(face)
