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
from OCC.Core.gp import gp_Cylinder, gp_Pln
from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepProj import BRepProj_Projection
from OCC.Core.TopoDS import TopoDS_Vertex
from OCC.Extend.ShapeFactory import make_face, make_vertex
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCCUtils.Common import vertex2pnt
from OCCUtils.Construct import vec_to_dir, dir_to_vec
from OCCUtils.Construct import make_polygon

if __name__ == '__main__':
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument("--pxyz", dest="pxyz",
                        default=[0.0, 0.0, 0.0], type=float, nargs=3)
    opt = parser.parse_args()
    print(opt, argvs)

    obj = dispocc(touch=True)
    face1 = make_face(gp_Cylinder(
        gp_Ax3(gp_Pnt(-10, 5., 1), gp_Dir(0, 0, 1)), 50), 0, 2 * np.pi, 0, 200)
    face2 = make_face(gp_Cylinder(
        gp_Ax3(gp_Pnt(0, 10, 1), gp_Dir(0, 0.1, 1)), 25), -np.pi, np.pi, 0, 200)
    pln = make_face(
        gp_Pln(gp_Ax3(gp_Pnt(100, 0, 100), gp_Dir(1, 0, 0))), -100, 100, -100, 100)
    obj.selected_shape = [pln, face1, face2]
    face = obj.make_comp_selcted()

    th = -15  # deg
    dx = np.cos(np.deg2rad(th))
    dy = np.sin(np.deg2rad(th))
    vz = gp_Vec(dx, dy, 0)
    vy = gp_Vec(0, 0, 1)
    vx = vy.Crossed(vz)
    axs = gp_Ax3(gp_Pnt(-101, 50, 100), vec_to_dir(vz), vec_to_dir(vx))
    pts = []
    pts.append(gp_Pnt(0, 0, 0))
    pts.append(gp_Pnt(10, 0, 0))
    pts.append(gp_Pnt(10, 10, 0))
    pts.append(gp_Pnt(0, 10, 0))
    poly = make_polygon(pts, closed=True)
    poly.Location(set_loc(gp_Ax3(), axs))

    obj.display.DisplayShape(face, color="BLUE1", transparency=0.9)
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
        obj.display.DisplayMessage(pnt[0], f"{i:d}", height=25)
        i += 1
        proj.Next()

    # Projection order depends on the size of the UV parameters of the base shape being projected

    obj.display.View_Top()
    obj.ShowOCC()
