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
from OCC.Core.TColgp import TColgp_Array1OfPnt, TColgp_HArray1OfPnt
from OCC.Core.GeomPlate import GeomPlate_BuildAveragePlane
from OCC.Core.GeomLib import geomlib
from OCC.Core.Geom2dAdaptor import geom2dadaptor
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Compound
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCCUtils.Common import vertex2pnt
from OCCUtils.Construct import make_box, make_line, make_wire, make_edge, make_face
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
    axs = gp_Ax3(gp_Pnt(1, 1, 1), gp_Dir(0.1, 0.2, 0.3), gp_Dir(0.2, 0.5, 0.1))
    pln = obj.make_plane_axs(axs, [-10, 10], [-10, 10])
    pnt = TColgp_Array1OfPnt(1, 9)
    pnt.SetValue(1, gp_Pnt(0, 0, 0))
    pnt.SetValue(2, gp_Pnt(1, 1, 1.1))
    pnt.SetValue(3, gp_Pnt(2, 2, 2.1))
    pnt.SetValue(4, gp_Pnt(-1, -1, -1.1))
    pnt.SetValue(5, gp_Pnt(-1, 2, -2))
    pnt.SetValue(6, gp_Pnt(1, -2, 2))
    pnt.SetValue(7, gp_Pnt(2, -4, 4))
    pnt.SetValue(8, gp_Pnt(-2, 4, -4.1))
    pnt.SetValue(9, gp_Pnt(3, 3, 2.9))
    print(geomlib.Inertia(pnt, axs.Location(), axs.XDirection(), axs.Direction()))

    pnt_h = TColgp_HArray1OfPnt(pnt.Lower(), pnt.Upper())
    for i in range(pnt.Lower(), pnt.Upper()):
        pnt_h.SetValue(i, pnt.Value(i))
    api = GeomPlate_BuildAveragePlane(pnt_h, 5, 0.1e-3, 1, 1)
    pnt_pln = api.Plane()
    pnt_ax1 = pnt_pln.Axis()
    pnt_pln_face = make_face(pnt_pln.Pln(), *api.MinMaxBox())

    obj.show_axs_pln(axs, scale=10)
    obj.display.DisplayShape(pln, transparency=0.9)
    [obj.display.DisplayShape(pnt.Value(i))
     for i in range(pnt.Lower(), pnt.Upper())]
    obj.display.DisplayShape(pnt_pln_face, transparency=0.9)
    obj.show_axs_pln(gp_Ax3(pnt_ax1.Location(), pnt_ax1.Direction()), scale=10)
    obj.ShowOCC()
