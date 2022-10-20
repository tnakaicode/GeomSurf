import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os
import time
import argparse
from linecache import getline, clearcache

sys.path.append(os.path.join("../"))
from src.base_occ import dispocc
from AirFoil.selig import uiuc_database, coord_database

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_ThruSections
from OCCUtils.Construct import make_polygon, make_face

if __name__ == '__main__':
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument("--pxyz", dest="pxyz",
                      default=[0.0, 0.0, 0.0], type=float, nargs=3)
    opt = parser.parse_args()
    print(opt, argvs)

    
    obj = dispocc()

    df = pd.read_csv("./wing001.csv", skiprows=2)
    print(df)
    
    api = BRepOffsetAPI_ThruSections()
    #api.SetSmoothing(True)

    rims = []
    for i, data_type in enumerate(df["type"]):
        px, py, pz = df["px"].values[i], df["py"].values[i], df["pz"].values[i]
        sz, sx, sy = df["sz"].values[i], df["sx"].values[i], df["sy"].values[i]
        deg = df["deg"].values[i]
        print(i, data_type, df["airfoil"][i], px, py, pz, sz)
        if data_type == "UIUC":
            xy = uiuc_database(df["airfoil"][i])
        pts = []
        for xy_2d in xy:
            x2, y2 = xy_2d[0] - sx, xy_2d[1] - sy
            pnt = gp_Pnt(x2*sz, y2*sz, px)
            pts.append(pnt)
        rim = make_polygon(pts, True)
        rims.append (rim)
        api.AddWire(rim)
        obj.display.DisplayShape(rim)
    api.Build()
    print(api.Check())
    
    obj.selected_shape = [
        make_face(rims[0]),
        api.Shape(),
        make_face(rims[-1])
    ]
    surf = obj.make_comp_selcted()
    obj.export_stp(surf, "wing001.step")
    obj.display.DisplayShape(surf)
    
    obj.show_axs_pln()
    obj.ShowOCC()
