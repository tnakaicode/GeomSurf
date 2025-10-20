import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import argparse
from linecache import getline, clearcache, updatecache

basename = os.path.dirname(__file__)

sys.path.append(os.path.join("./"))
from base_occ import dispocc, set_loc

import logging

logging.getLogger("matplotlib").setLevel(logging.ERROR)

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Circ
from OCC.Core.TColgp import TColgp_Array2OfPnt
from OCC.Core.GeomAbs import GeomAbs_G2
from OCC.Core.GeomAPI import GeomAPI_PointsToBSplineSurface
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Compound
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCCUtils.Common import vertex2pnt
from OCCUtils.Construct import make_box, make_line, make_wire, make_edge, make_face
from OCCUtils.Construct import make_plane, make_polygon
from OCCUtils.Construct import point_to_vector, vector_to_point
from OCCUtils.Construct import dir_to_vec, vec_to_dir


def spl_face(px, py, pz, axs=gp_Ax3()):
    nx, ny = px.shape
    pnt_2d = TColgp_Array2OfPnt(1, nx, 1, ny)
    for row in range(pnt_2d.LowerRow(), pnt_2d.UpperRow() + 1):
        for col in range(pnt_2d.LowerCol(), pnt_2d.UpperCol() + 1):
            i, j = row - 1, col - 1
            pnt = gp_Pnt(px[i, j], py[i, j], pz[i, j])
            pnt_2d.SetValue(row, col, pnt)
            # print (i, j, px[i, j], py[i, j], pz[i, j])

    api = GeomAPI_PointsToBSplineSurface(pnt_2d, 3, 8, GeomAbs_G2, 0.001)
    api.Interpolate(pnt_2d)
    # surface = BRepBuilderAPI_MakeFace(curve, 1e-6)
    # return surface.Face()
    face = BRepBuilderAPI_MakeFace(api.Surface(), 1e-6).Face()
    face.Location(set_loc(gp_Ax3(), axs))
    return face


if __name__ == "__main__":
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument(
        "--pxyz", dest="pxyz", default=[0.0, 0.0, 0.0], type=float, nargs=3
    )
    opt = parser.parse_args()
    print(opt)

    obj = dispocc(touch=True)
    axs = gp_Ax3()

    nx, ny = 501, 301
    lx, ly = 300, 400
    sx, sy = 0.0, 0.0
    rx, ry = 600.0, 500.0
    px = np.linspace(-1, 1, nx) * lx + sx
    py = np.linspace(-1, 1, ny) * ly + sy
    mesh = np.meshgrid(px, py)
    data = mesh[0] ** 2 / (2 * rx) + mesh[1] ** 2 / (2 * ry) + 10.0
    face = spl_face(*mesh, data, axs)

    pts = [
        gp_Pnt(-50, -50, 0),
        gp_Pnt(50, -50, 0),
        gp_Pnt(50, 50, 0),
        gp_Pnt(-50, 50, 0),
    ]
    wire = make_polygon(pts, closed=True)

    from OCC.Core.BRepProj import BRepProj_Projection
    from OCC.Core.GeomProjLib import geomprojlib
    from OCC.Core.ProjLib import projlib, ProjLib_Plane
    from OCC.Core.gp import gp_Pln

    pl = gp_Pln(gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0.1, 0.2, 0.8)))
    p2 = projlib.Project(pl, gp_Pnt(10, 11, 13))
    print(p2.X(), p2.Y())

    proj = BRepProj_Projection(wire, face, axs.Direction())
    obj.display.DisplayShape(wire)
    obj.display.DisplayShape(proj.Current())

    wire = make_wire(make_edge(gp_Circ(axs.Ax2(), 150)))
    proj = BRepProj_Projection(wire, face, gp_Dir(0.1, 0.2, 1.0))
    obj.display.DisplayShape(wire)
    obj.display.DisplayShape(proj.Current())
    obj.display.DisplayShape(make_face(face, proj.Current()), color="RED")

    obj.show_axs_pln()
    # obj.display.DisplayShape(face)
    obj.ShowOCC()
