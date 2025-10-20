import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import argparse
from linecache import getline, clearcache, updatecache

basename = os.path.dirname(__file__)

sys.path.append(os.path.join("./"))
from src.base_occ import dispocc, spl_curv

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Elips
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Compound
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Core.TColgp import TColgp_Array1OfPnt, TColgp_HArray1OfPnt
from OCC.Core.TColStd import TColStd_HArray1OfReal
from OCC.Core.Geom import Geom_Ellipse
from OCC.Core.GeomAPI import GeomAPI_PointsToBSpline, GeomAPI_Interpolate
from OCC.Core.BRep import BRep_Tool
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

    obj = dispocc(disp=True, touch=True)
    axs = gp_Ax3()

    nt = 10
    t0 = 0
    rx = 11
    ry = 5

    pt = np.linspace(0, 2 * np.pi, nt)
    px = np.empty_like(pt)
    py = np.empty_like(pt)
    for i, t in enumerate(pt):
        x = rx * np.cos(t)
        y = ry * np.sin(t)
        px[i] = x * np.cos(np.deg2rad(t0)) + y * np.sin(np.deg2rad(t0))
        py[i] = x * -np.sin(np.deg2rad(t0)) + y * np.cos(np.deg2rad(t0))

    pts = TColgp_HArray1OfPnt(1, nt - 1)
    dat = TColStd_HArray1OfReal(1, nt - 1)
    for i in range(nt - 1):
        pts.SetValue(i + 1, gp_Pnt(px[i], py[i], 0))
    for i in range(nt - 1):
        dat.SetValue(i + 1, 1)

    api = GeomAPI_Interpolate(pts, True, 0.1E-3)
    api.Perform()
    crv = api.Curve()
    elp = make_edge(gp_Elips(axs.Ax2(), rx, ry))

    obj.display.DisplayShape(crv)
    obj.display.DisplayShape(elp, color="BLUE1")

    obj.new_2Dfig(aspect="equal")
    obj.axs.plot(px, py)
    obj.axs.set_xlim(-1500, 1500)
    obj.axs.set_ylim(-1500, 1500)
    obj.SavePng(obj.tempname + "_Plot.png")

    crv1 = obj.make_ellips(gp_Ax3(gp_Pnt(0, 0, 10),
                                  gp_Dir(0, 0, 1)),
                           [10, 11], 30, 200)
    crv2 = obj.make_ellips(gp_Ax3(gp_Pnt(0, 0, 20),
                                  gp_Dir(0, 0.1, 1)),
                           [5, 11], 30, 200)
    print(BRep_Tool.Curve(crv2))

    obj.display.DisplayShape(crv1)
    obj.display.DisplayShape(crv2)

    obj.show_axs_pln(scale=5)
    obj.ShowOCC()
