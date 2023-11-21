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
from OCCUtils.Construct import dir_to_vec, vec_to_dir


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
        gp_Pnt(0, 60, 0),
        gp_Pnt(-50, 50, 0),
    ]
    edg = [make_edge(pts[(i + 1) % len(pts)], pts[i]) for i in range(len(pts))]
    fil = []
    pln = gp_Pln()
    radii = 20.0
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

    f = ChFi2d_AnaFilletAlgo()
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

    e1, e2 = edg[1], edg[2]

    f.Init(e1, e2, pln)
    f.Perform(radii)
    fil.append(f.Result(e1, e2))
    #
    # edg.insert(1, fil[0])
    # edg.insert(3, fil[1])

    # poly = make_wire([edg[0], fil[0], edg[1], fil[1], edg[2], edg[3], edg[4]])
    # obj.display.DisplayShape(poly)
    obj.display.DisplayShape(edg, color="BLUE1")
    obj.display.DisplayShape(fil)
    [obj.display.DisplayShape(edge_midpoint(e)) for e in edg + fil]
    obj.ShowOCC()
