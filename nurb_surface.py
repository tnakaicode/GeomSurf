import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import argparse

sys.path.append(os.path.join("./"))
from src.base_occ import dispocc
from rnd_sample import cube01_sample, ball01_sample

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.BRepTools import BRepTools_Modification, BRepTools_NurbsConvertModification
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Compound, TopoDS_Vertex
from OCCUtils.Construct import make_vertex


if __name__ == '__main__':
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument("--pxyz", dest="pxyz",
                        default=[0.0, 0.0, 0.0], type=float, nargs=3)
    opt = parser.parse_args()
    print(opt, argvs)

    obj = dispocc(touch=True)
    brp = BRepTools_NurbsConvertModification()

    dat, seed = cube01_sample(2, 100, -1)
    dat = dat.T
    for xy in dat:
        pnt = gp_Pnt(*xy, xy[0]**2 + xy[1]**2)
        obj.display.DisplayShape(pnt)

    dat, seed = ball01_sample(1000, -1)
    dat = dat.T
    for xyz in dat:
        pnt = gp_Pnt(*xyz)
        obj.display.DisplayShape(pnt, color="RED")
    obj.show_axs_pln(scale=1.0)
    obj.ShowOCC()

    # for idx, x in enumerate (dat[0,:]):
    #    pnt = gp_Pnt(dat[0,idx], dat[1, idx], dat[2,idx])
    #    v = make_vertex(pnt)
    #    tol = 1.0 * 10e-06
    #    brp.NewPoint(v, pnt)
