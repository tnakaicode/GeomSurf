import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import argparse
from linecache import getline, clearcache, updatecache

basename = os.path.dirname(__file__)

sys.path.append(os.path.join("./"))
from base_occ import dispocc

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Lin
from OCC.Core.Geom import Geom_Line, Geom_Plane
from OCC.Core.GeomAPI import GeomAPI_IntCS, GeomAPI_IntSS, GeomAPI_ProjectPointOnSurf
from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepProj import BRepProj_Projection
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

    obj = dispocc(disp=False, touch=True)
    obj.proj_pnt_pln

    pnt = gp_Pnt(10, 11, 100)
    vec = gp_Dir(1, 2, 3)
    axs = gp_Ax3(gp_Pnt(0, 0, 1),
                 gp_Dir(0, 0, 1),
                 gp_Dir(1, 0, 0))
    print(vec.X(), vec.Y(), vec.Z())
    vec = gp_Vec(vec.XYZ()).Scaled(11)
    vec = gp_Vec(0, 0, 0)

    lin = gp_Lin(pnt, axs.Direction())
    api = GeomAPI_IntCS(Geom_Line(lin), Geom_Plane(axs))
    dst = np.inf
    sxy = pnt
    num = api.NbPoints()
    for i in range(num):
        p0 = api.Point(i + 1)
        dst0 = pnt.Distance(p0)
        if dst0 < dst:
            dst = dst0
            sxy = p0
    print(sxy)

    lin = gp_Lin(pnt.Translated(vec), axs.Direction())
    api = GeomAPI_IntCS(Geom_Line(lin), Geom_Plane(axs))
    dst = np.inf
    sxy = pnt
    num = api.NbPoints()
    for i in range(num):
        p0 = api.Point(i + 1)
        dst0 = pnt.Distance(p0)
        if dst0 < dst:
            dst = dst0
            sxy = p0
    print(sxy)
