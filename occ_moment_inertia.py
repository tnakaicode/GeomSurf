import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import argparse
from linecache import getline, clearcache, updatecache

basename = os.path.dirname(__file__)

sys.path.append(os.path.join("./"))
from base_occ import dispocc, set_trf

import logging

logging.getLogger("matplotlib").setLevel(logging.ERROR)

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3, gp_XYZ
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Compound
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Core.Bnd import Bnd_Box, Bnd_OBB
from OCC.Core.BRep import BRep_Builder
from OCC.Core.BRepBndLib import brepbndlib
from OCC.Core.BRepGProp import brepgprop
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCC.Core.GProp import GProp_GProps
from OCC.Extend.DataExchange import read_step_file
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCCUtils.Construct import make_box, make_line, make_wire, make_edge
from OCCUtils.Construct import make_plane, make_polygon, make_vertex, make_box
from OCCUtils.Construct import point_to_vector, vector_to_point
from OCCUtils.Construct import dir_to_vec, vec_to_dir

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

    datfile = "./PlotGallery/occ_gallery/assets/models/bunny.pcd"
    dat = np.loadtxt(datfile, skiprows=10)
    pts = [gp_Pnt(*xyz) for xyz in dat]
    ver = [make_vertex(p) for p in pts]

    build = BRep_Builder()
    shape = TopoDS_Compound()
    build.MakeCompound(shape)
    for v in ver:
        build.Add(shape, v)

    # shape = read_step_file("./PlotGallery/occ_gallery/assets/models/as1_pe_203.stp")

    obj.display.DisplayShape(shape)
    # obj.display.DisplayShape(ver)

    # ver = [
    #    make_vertex(gp_Pnt(0, 0, 0)),
    #    make_vertex(gp_Pnt(1, 1, 1)),
    # ]
    ax0 = gp_Ax3(gp_Pnt(-1, -1, -1), gp_Dir(0.1, 0.1, 0.5))

    bbox = Bnd_Box()
    bbox.Add(ax0.Location(), ax0.Direction())
    brepbndlib.AddOptimal(shape, bbox, True, True)
    # bbox.Set(axs.Location(), axs.Direction())
    # for v in ver:
    #    brepbndlib.Add(v, bbox)
    print(bbox.IsOpen())
    print("Xmax", bbox.IsOpenXmax())
    print("Xmin", bbox.IsOpenXmin())
    print("Ymax", bbox.IsOpenYmax())
    print("Ymin", bbox.IsOpenYmin())
    print("Zmax", bbox.IsOpenZmax())
    print("Zmin", bbox.IsOpenZmin())
    # bbox = bbox.FinitePart()
    xmin, ymin, zmin, xmax, ymax, zmax = bbox.Get()
    p0 = bbox.CornerMin()
    p1 = bbox.CornerMax()
    # obj.display.DisplayShape(make_box(p0, p1), color="BLUE1", transparency=0.9)
    # obj.display.DisplayShape(p0, color="RED")
    # obj.display.DisplayShape(p1, color="RED")
    print(p0)
    print(p1)

    bobb = Bnd_OBB()
    bbox.Set(ax0.Location(), ax0.Direction())
    brepbndlib.AddOBB(shape, bobb, True, True, True)
    p2 = gp_Pnt(bobb.Center() - gp_XYZ(bobb.XHSize(), bobb.YHSize(), bobb.ZHSize()))
    print(bobb.Center())
    print(bobb.Position())
    ax_obb = bobb.Position()
    # ax_obb.YReverse()
    p0_obb = gp_Pnt(-bobb.XHSize(), -bobb.YHSize(), -bobb.ZHSize())
    p1_obb = gp_Pnt(bobb.XHSize(), bobb.YHSize(), bobb.ZHSize())
    trf = set_trf(axs, ax_obb)
    p0_obb.Transform(trf)
    p1_obb.Transform(trf)
    box_obb = make_box(
        gp_Pnt(-bobb.XHSize(), -bobb.YHSize(), -bobb.ZHSize()),
        bobb.XHSize() * 2,
        bobb.YHSize() * 2,
        bobb.ZHSize() * 2,
    )
    box_obb.Move(TopLoc_Location(trf))
    obj.display.DisplayShape(box_obb, transparency=0.9)
    print(p0_obb)
    print(p1_obb)
    print(dir_to_vec(ax0.Direction()))
    print(dir_to_vec(ax_obb.Direction()))
    print(dir_to_vec(ax_obb.XDirection()))
    print(dir_to_vec(ax_obb.YDirection()))

    # Sets this bounding box so that it bounds the half-line defined by point p and direction d,
    # i.e. all points m defined by m=p+u*d, where u is greater than or equal to 0, are inside the bounding volume.
    # this involves first setting this box to be void and then adding the half-line.

    # このバウンディングボックスを、点pと方向dで定義される半直線を境界とするように、
    # すなわち、m=p+u*dで定義されるすべての点m（uは0以上）がバウンディングボリューム内に入るように設定します。

    # < gp_Vec: -0.22113772383568767, -0.3976834059076125, -0.8904751629116291, magnitude: 1.0 >
    # < gp_Vec: 0.8695499303663634, 0.3330186365376008, -0.36466629446452803, magnitude: 1.0 >
    # < gp_Vec: 0.4415665586257905, -0.8549540902202636, 0.27216259463603343, magnitude: 1.0 >

    # < gp_Vec: 0.0, 0.7071067811865475, 0.7071067811865475, magnitude: 0.9999999999999999 >
    # < gp_Vec: 0.9848287924639173, -0.08043250350807839, -0.15376235531977758, magnitude: 1.0 >
    # < gp_Vec: -0.17319370016637425, -0.40055534077678684, -0.8997551673638071, magnitude: 1.0 >
    # < gp_Vec: 0.010779228021646171, 0.9127354662521968, -0.4084088354683847, magnitude: 0.9999999999999999 >

    prop = GProp_GProps()
    brepgprop.VolumePropertiesGK(shape, prop)
    print(prop.Mass())
    print(prop.MomentOfInertia(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1))))
    print(prop.MomentOfInertia(gp_Ax1(gp_Pnt(0.1, 0.1, 0.1), gp_Dir(0, 1, 1))))
    print(prop.RadiusOfGyration(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1))))
    print(prop.StaticMoments())

    obj.show_axs_pln(ax_obb, scale=0.1)
    obj.ShowOCC()
