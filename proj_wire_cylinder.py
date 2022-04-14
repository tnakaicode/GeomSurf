import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import argparse
from linecache import getline, clearcache

sys.path.append(os.path.join("./"))
from src.base_occ import dispocc, spl_face, set_loc

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Pnt, gp_Dir, gp_Vec
from OCC.Core.gp import gp_Cylinder
from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepProj import BRepProj_Projection
from OCC.Core.TopoDS import TopoDS_Vertex
from OCC.Extend.ShapeFactory import make_face, make_vertex
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCCUtils.Construct import vec_to_dir, dir_to_vec, vector_to_point
from OCCUtils.Construct import make_polygon

if __name__ == '__main__':
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument("--pxyz", dest="pxyz",
                        default=[0.0, 0.0, 0.0], type=float, nargs=3)
    opt = parser.parse_args()
    print(opt, argvs)

    face = make_face(gp_Cylinder(gp_Ax3(), 50), 0, 2 * np.pi, 0, 200)

    th = 15  # deg
    dx = np.cos(np.deg2rad(th))
    dy = np.sin(np.deg2rad(th))
    axs = gp_Ax3(gp_Pnt(0, 0, 100), gp_Dir(dx, dy, 0))
    pts = []
    pts.append(gp_Pnt(0, 0, 0))
    pts.append(gp_Pnt(10, 0, 0))
    pts.append(gp_Pnt(10, 10, 0))
    pts.append(gp_Pnt(0, 10, 0))
    poly = make_polygon(pts, closed=True)
    poly.Location(set_loc(gp_Ax3(), axs))

    obj = dispocc(touch=True)
    obj.display.DisplayShape(face, color="BLUE", transparency=0.9)
    obj.display.DisplayShape(poly)
    obj.show_axs_pln(axs, scale=10)

    i = 0
    idx = 0
    proj = BRepProj_Projection(poly, face, axs.Direction())
    # while proj.More() or i ==idx:
    #    proj_poly = proj.Current()
    #    i+=1
    #    proj.Next()

    #obj.display.DisplayShape(proj_poly, color="BLUE")

    while proj.More():
        brt = BRep_Tool()
        proj_poly = proj.Current()
        pnt = [brt.Pnt(v) for v in TopologyExplorer(proj_poly).vertices()]
        obj.display.DisplayShape(proj_poly)
        obj.display.DisplayMessage(pnt[0], f"{i:d}")
        i += 1
        proj.Next()

    obj.ShowOCC()
