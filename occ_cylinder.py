import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import argparse
from linecache import getline, clearcache

sys.path.append(os.path.join("./"))
from src.base_occ import dispocc, spl_curv_pts

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_XOY, gp_Pnt2d, gp_Dir2d, gp_Lin2d, gp_Vec2d
from OCC.Core.gp import gp_Circ
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_MakePipe
from OCC.Core.BRepAlgo import BRepAlgo_FaceRestrictor
from OCC.Core.BRepProj import BRepProj_Projection
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeEdge2d, BRepBuilderAPI_MakeFace
from OCC.Core.Geom import Geom_CylindricalSurface
from OCC.Core.GeomAPI import GeomAPI_PointsToBSplineSurface
from OCC.Core.GeomAPI import GeomAPI_PointsToBSpline
from OCC.Core.gce import gce_MakeTranslation
from OCC.Core.GCE2d import GCE2d_MakeSegment
from OCC.Core.Geom2d import Geom2d_Curve, Geom2d_Line
from OCC.Extend.DataExchange import write_step_file, read_step_file
from OCC.Extend.DataExchange import write_stl_file, read_stl_file
from OCCUtils.Construct import make_face, make_polygon, make_wire, make_edge, make_edge2d
from OCCUtils.Construct import dir_to_vec, vec_to_dir

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)


if __name__ == '__main__':
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument("--pxyz", dest="pxyz",
                        default=[0.0, 0.0, 0.0], type=float, nargs=3)
    opt = parser.parse_args()
    print(opt, argvs)

    obj = dispocc(touch=False)
    axs = gp_Ax3()
    axr = gp_Ax3(gp_Pnt(0, 0, 100),
                 gp_Dir(-1, 0, 0))
    cir = make_edge(gp_Circ(axr.Ax2(), 100), -np.pi / 2, np.pi / 2)
    lin = make_polygon([gp_Pnt(-200, 0, 0), gp_Pnt(200, 0, 0)])
    api = BRepOffsetAPI_MakePipe(lin, cir)
    api.Build()

    obj.display.DisplayShape(cir)
    obj.display.DisplayShape(lin)
    obj.display.DisplayShape(api.Shape())
    obj.show_axs_pln(axs)
    obj.ShowOCC()
