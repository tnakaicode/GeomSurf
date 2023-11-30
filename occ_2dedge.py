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
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3, gp_Trsf
from OCC.Core.gp import gp_Pnt2d, gp_Trsf2d
from OCC.Core.TColgp import TColgp_Array1OfPnt, TColgp_Array1OfPnt2d
from OCC.Core.TColStd import TColStd_Array1OfReal, TColStd_Array1OfInteger
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Compound
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Core.Geom import Geom_BezierCurve, Geom_BSplineCurve, Geom_Plane
from OCC.Core.Geom2d import Geom2d_BezierCurve, Geom2d_BSplineCurve
from OCC.Core.GeomAPI import GeomAPI_PointsToBSpline, GeomAPI_PointsToBSplineSurface
from OCC.Core.Geom2dAPI import Geom2dAPI_PointsToBSpline
from OCC.Core.GeomAdaptor import GeomAdaptor_Surface
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
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
    axs = gp_Ax3()
    ax1 = gp_Ax3(gp_Pnt(1, 0, 1), gp_Dir(0.1, 0.2, 1))
    ax2 = gp_Ax3(gp_Pnt(1, 0, 2), ax1.Direction())
    trf = gp_Trsf()
    trf.SetTransformation(ax1, gp_Ax3())

    dat = [[pt, np.sin(pt)/2+1, 0] for pt in np.linspace(0, 2 * np.pi, 10)]
    pts = [gp_Pnt(*d).Transformed(trf) for d in dat]
    pts_2d = [gp_Pnt2d(d[0], d[1]).Transformed(gp_Trsf2d(trf)) for d in dat]

    col_pts = TColgp_Array1OfPnt(1, len(dat))
    col_pts_2d = TColgp_Array1OfPnt2d(1, len(dat))
    for i, xyz in enumerate(dat):
        col_pts.SetValue(i + 1, pts[i])
        col_pts_2d.SetValue(i + 1, gp_Pnt2d(xyz[0], xyz[1]))

    curv = GeomAPI_PointsToBSpline(col_pts, 3, 8).Curve()
    curv_2d = Geom2dAPI_PointsToBSpline(col_pts_2d, 3, 8).Curve()

    pln = obj.make_plane_axs(ax1, rx=[-1, 10], ry=[1.5, -1.5])
    pln_surf = BRepAdaptor_Surface(pln, True).Surface()

    obj.display.DisplayShape(make_polygon(pts))
    # obj.display.DisplayShape(make_face(pln, make_polygon(pts)))
    obj.display.DisplayShape(make_face(pln, make_wire(make_edge(curv, curv.FirstParameter(), curv.LastParameter()))))
    # [obj.display.DisplayShape(p) for p in pts_2d]
    obj.display.DisplayShape(curv)
    obj.display.DisplayShape(curv_2d)
    obj.display.DisplayShape(make_edge(curv_2d, Geom_Plane(ax2)), color="BLUE1")

    obj.show_axs_pln(scale=1)
    obj.ShowOCC()
