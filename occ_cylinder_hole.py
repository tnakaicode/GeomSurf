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
    surf = Geom_CylindricalSurface(gp_Ax3(), 100)
    face = make_face(surf.Cylinder(), 0, np.pi, -10, 50)

    dat = [
        [np.pi / 8, 5],
        [np.pi / 4, 8],
        [np.pi / 6, 10],
        [np.pi / 8, 10],
    ]
    edg = []
    for i, _ in enumerate(dat):
        i0, i1 = i, (i + 1) % (len(dat))
        p0, p1 = gp_Pnt2d(*dat[i0]), gp_Pnt2d(*dat[i1])
        li = gp_Lin2d(p0, gp_Dir2d(gp_Vec2d(p0, p1)))
        ci = GCE2d_MakeSegment(li, p0, p1).Value()
        ei = BRepBuilderAPI_MakeEdge(ci, surf, 0, p1.Distance(p0)).Edge()
        print(ei, ci, p1.Distance(p0))
        # obj.display.DisplayShape(ei)
        edg.append(ei)
    rim = make_wire(edg)

    api_face = BRepAlgo_FaceRestrictor()
    api_face.Init(face, False, True)
    api_face.Add(rim)
    api_face.Perform()
    print(api_face.IsDone())
    while api_face.More():
        print(api_face.Current())
        # obj.display.DisplayShape(api_face.Current(), transparency=0.9)
        api_face.Next()

    rim = make_polygon([gp_Pnt(0, 0, 0), gp_Pnt(20, 0, 0),
                       gp_Pnt(30, 0, 20), gp_Pnt(0, 0, 20)], closed=True)
    rim_proj = obj.proj_rim_pln(rim, face, gp_Ax3(gp_Pnt(0, 0, 0),
                                                  gp_Dir(0, 1, 0)), 1)
    api_face = BRepBuilderAPI_MakeFace(surf, rim_proj)
    # api_face.Init(face)
    # api_face.Add(rim)
    api_face.Build()
    # print(api_face.IsDone())
    # face = api_face.Shape()
    # face = make_face(face, rim)
    # obj.export_stp(face)
    obj.display.DisplayShape(rim)
    obj.display.DisplayShape(face)
    obj.display.DisplayShape(api_face.Shape(), color="RED")
    obj.ShowOCC()
