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
from OCC.Core.gp import gp_Pln, gp_Trsf
from OCC.Core.ChFi2d import ChFi2d_AnaFilletAlgo
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Compound, TopoDS_Edge
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCCUtils.Common import vertex2pnt, curve_length, midpoint
from OCCUtils.Construct import make_box, make_line, make_wire, make_edge
from OCCUtils.Construct import make_plane, make_polygon, make_face
from OCCUtils.Construct import point_to_vector, vector_to_point
from OCCUtils.Construct import dir_to_vec, vec_to_dir, vertex2pnt


def edge_length(edg=TopoDS_Edge()):
    from OCC.Core.GProp import GProp_GProps
    from OCC.Core.BRepGProp import brepgprop
    prop = GProp_GProps()
    brepgprop.LinearProperties(edg, prop)
    return prop.Mass()


def edge_midpoint(edg=TopoDS_Edge()):
    from OCC.Core.BRepGProp import BRepGProp_EdgeTool
    from OCC.Core.BRepAdaptor import BRepAdaptor_Curve
    tool = BRepGProp_EdgeTool()
    u0 = tool.FirstParameter(BRepAdaptor_Curve(edg))
    u1 = tool.LastParameter(BRepAdaptor_Curve(edg))
    return tool.Value(BRepAdaptor_Curve(edg), (u0 + u1) / 2)


def edge_1stpoint(edg=TopoDS_Edge()):
    from OCC.Core.BRepGProp import BRepGProp_EdgeTool
    from OCC.Core.BRepAdaptor import BRepAdaptor_Curve
    tool = BRepGProp_EdgeTool()
    u0 = tool.FirstParameter(BRepAdaptor_Curve(edg))
    return tool.Value(BRepAdaptor_Curve(edg), u0)


def edge_endpoint(edg=TopoDS_Edge()):
    from OCC.Core.BRepGProp import BRepGProp_EdgeTool
    from OCC.Core.BRepAdaptor import BRepAdaptor_Curve
    tool = BRepGProp_EdgeTool()
    u1 = tool.LastParameter(BRepAdaptor_Curve(edg))
    return tool.Value(BRepAdaptor_Curve(edg), u1)


def edge_rndpont(edg=TopoDS_Edge(), n=3):
    from OCC.Core.BRepGProp import BRepGProp_EdgeTool
    from OCC.Core.BRepAdaptor import BRepAdaptor_Curve
    from OCC.Core.BRepGProp import brepgprop
    from OCC.Core.GProp import GProp_GProps
    from OCC.Core.GeomAPI import GeomAPI_PointsToBSpline
    from OCC.Core.TColgp import TColgp_Array1OfPnt
    tool = BRepGProp_EdgeTool()
    c = BRepAdaptor_Curve(edg)
    u0 = tool.FirstParameter(c)
    u1 = tool.LastParameter(c)
    if n == 2:
        return edg
    elif n >= 2:
        prop = GProp_GProps()
        brepgprop.LinearProperties(edg, prop)
        lth = prop.Mass()
        pts = TColgp_Array1OfPnt(1, n)
        for i, t in enumerate(np.linspace(u0, u1, n)):
            p = tool.Value(c, t)
            x = p.X() + np.random.uniform(0, lth / 10)
            y = p.Y() + np.random.uniform(0, lth / 10)
            pts.SetValue(i + 1, gp_Pnt(x, y, p.Z()))
        pts.SetValue(1, tool.Value(c, u0))
        pts.SetValue(n, tool.Value(c, u1))
        curv = GeomAPI_PointsToBSpline(pts).Curve()
        return make_edge(curv, 0, 1)
    elif n == 0:
        return tool.Value(c, u1)
    elif n == 1:
        return tool.Value(c, u1)


def make_fillet(e1=TopoDS_Edge(), e2=TopoDS_Edge(), radii=10, pln=gp_Pln()):
    f = ChFi2d_AnaFilletAlgo()
    f.Init(e1, e2, pln)
    if radii == 0 or radii == None or f.Perform(radii) != True:
        return [e1, e2]
    else:
        f2 = f.Result(e1, e2)
        print("f2", edge_length(f2))
        print("e1->f2", edge_endpoint(e1), edge_1stpoint(f2))
        print("f2->e2", edge_endpoint(f2), edge_1stpoint(e2))
        return [e1, f2, e2]


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

    data = [
        [[-50, -50, 0], 10],
        [[50, -50, 0], 10],
        [[50, 50, 0], 10],
        [[20, 60, 0], 500],
        [[40, -40, 0], 2],
        [[-50, 50, 0], 10],
    ]

    pts = [gp_Pnt(*dat) for (dat, r) in data]
    pln = gp_Pln(axs)
    edg = [make_edge(pts[i], pts[(i + 1) % len(pts)]) for i in range(len(pts))]
    edg[0] = edge_rndpont(edg[0], 5)

    poly = make_wire(edg)
    obj.display.DisplayShape(poly, color="BLUE1")

    obj.ShowOCC()
