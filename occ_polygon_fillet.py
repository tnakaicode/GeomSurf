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
from OCCUtils.Construct import make_plane, make_polygon
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


def make_fillet(e1=TopoDS_Edge(), e2=TopoDS_Edge(), radii=10):
    f = ChFi2d_AnaFilletAlgo()
    f.Init(e1, e2, pln)
    f.Perform(radii)
    f2 = f.Result(e1, e2)
    print("f2", edge_length(f2))
    print("e1->f2", edge_endpoint(e1), edge_1stpoint(f2))
    print("f2->e2", edge_endpoint(f2), edge_1stpoint(e2))
    return e1, f2, e2


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

    pts = [
        gp_Pnt(-50, -50, 0),
        gp_Pnt(50, -50, 0),
        gp_Pnt(50, 50, 0),
        gp_Pnt(20, 60, 0),
        gp_Pnt(40, -40, 0),
        gp_Pnt(-50, 50, 0),
    ]
    edg = [make_edge(pts[i], pts[(i + 1) % len(pts)]) for i in range(len(pts))]
    fil = []
    pln = gp_Pln()
    radii = 5.0
    poly = make_wire(edg)
    obj.display.DisplayShape(edg, color="BLUE1")

    edg = [make_edge(pts[0], pts[1])]
    for i, p in enumerate(pts + [pts[0]]):
        i1, i2 = (i + 1) % (len(pts)), (i + 2) % (len(pts))
        e1, e2 = edg[-1], make_edge(pts[i1], pts[i2])
        e1, f2, e2 = make_fillet(e1, e2, radii)
        edg[-1] = e1
        edg += [f2, e2]
    edg.pop(0)
    edg.pop(-1)
    edg.pop(-1)

    # f = ChFi2d_AnaFilletAlgo()
    # edg = []
    # e0 = make_edge(pts[0], pts[1])
    # e1 = make_edge(pts[1], pts[2])
    # f.Init(e0, e1, pln)
    # f.Perform(radii)
    # f01 = f.Result(e0, e1)
    # edg += [e0, f01, e1]
    # print("0->1->2", pts[0], pts[1], pts[2], edge_length(f01))
    # print("e0->f01", edge_endpoint(e0), edge_1stpoint(f01))
    # print("f01->e1", edge_1stpoint(e1), edge_endpoint(f01))

    # for i, p in enumerate(pts):
    #    i1, i2, i3 = i, (i + 1) % (len(pts)), (i + 2) % (len(pts))
    #    p1, p2, p3 = pts[i1], pts[i2], pts[i3]
#
    #    # Making the edges
    #    ed1 = make_edge(p3, p2)
    #    ed2 = make_edge(p2, p1)
#
    #    # Making the 2dFillet
    #    f = ChFi2d_AnaFilletAlgo(ed1, ed2, pln)
    #    f.Perform(radii)
    #    fillet2d = f.Result(ed1, ed2)
#
    #    #edg += [ed1, fillet2d]
    #    obj.display.DisplayShape(make_wire([ed1, fillet2d, ed2]), color=obj.colors[0])

    # for i, e in enumerate(edg[:2]):
    #    print(i, e)
    #    print(edge_length(edg[i]), edge_length(edg[i+1]))
    #    f.Init(edg[i].Reversed(), edg[i+1].Reversed(), pln)
    #    #try:
    #    #    f.Init(edg[i], edg[i+1], pln)
    #    #except:
    #    #    f.Init(edg[i].Reversed(), edg[i+1], pln)
    #    print(f.Perform(radii))
    #    fi = f.Result(edg[i], edg[i+1])
    #    print(edge_length(fi))
    #    print(edge_length(edg[i]), edge_length(edg[i+1]))

    # e2 = make_edge(pts[2], pts[3])
    # print("e1", edge_1stpoint(e1), edge_endpoint(e1))
    # print("e2", edge_1stpoint(e2), edge_endpoint(e2))
    # f = ChFi2d_AnaFilletAlgo()
    # f.Init(e1, e2, pln)
    # f.Perform(radii)
    # f12 = f.Result(e1, e2)
    # edg += [f12, e2]
    # print("1->2->3", pts[1], pts[2], pts[3], edge_length(f12))
    # print("e1->f12", edge_endpoint(e1), edge_1stpoint(f12))
    # print("f12->e3", edge_1stpoint(e2), edge_endpoint(f12))

    # e2, e3 = edg[-1], make_edge(pts[4], pts[3])
    # f.Init(e2, e3, pln)
    # f.Perform(radii)
    # f23 = f.Result(e2, e3)
    # edg += [f23, e3]
    #
    # edg.insert(1, fil[0])
    # edg.insert(3, fil[1])

    poly = make_wire(edg)
    obj.display.DisplayShape(poly)
    [obj.display.DisplayShape(edge_midpoint(e)) for e in edg + fil]
    obj.ShowOCC()
