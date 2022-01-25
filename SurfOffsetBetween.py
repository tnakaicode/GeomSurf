import numpy as np
import matplotlib.pyplot as plt
import sys
import time
import os
from optparse import OptionParser

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax3
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakePrism
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_ThruSections, BRepOffsetAPI_MakeOffset, BRepOffsetAPI_MakeEvolved, BRepOffsetAPI_MakePipe, BRepOffsetAPI_MakePipeShell
from OCC.Core.BRepOffset import BRepOffset_MakeOffset, BRepOffset_Skin, BRepOffset_Interval
from OCC.Core.GeomAbs import GeomAbs_C0, GeomAbs_C1, GeomAbs_C2
from OCC.Core.GeomAbs import GeomAbs_G1, GeomAbs_G2
from OCC.Core.GeomAbs import GeomAbs_Intersection, GeomAbs_Arc
from OCCUtils.Construct import make_polygon

sys.path.append(os.path.join('./'))
from src.base_occ import dispocc, occ_to_grasp_cor, set_loc


def make_offsetface(dat, axs=gp_Ax3()):
    pts = []
    for xyz in dat + [pts[0]]:
        pts.append(gp_Pnt(*xyz))
    face = dispocc.make_FaceByOrder(pts)


if __name__ == '__main__':
    obj = dispocc()

    ax1 = gp_Ax3()
    pt1 = np.loadtxt(obj.rootname + "_pln1.txt")
    print(pt1)
    pts = []
    for xyz in pt1 + [pt1[0]]:
        pts.append(gp_Pnt(*xyz))
    br1 = make_polygon(pts, closed=True)
    br1.Location(set_loc(gp_Ax3(), ax1))
    fc1 = obj.make_FaceByOrder(pts)
    fc1.Location(set_loc(gp_Ax3(), ax1))
    print(fc1)
    obj.display.DisplayShape(br1)
    obj.display.DisplayShape(fc1)

    ax2 = gp_Ax3()
    ax2.SetLocation(gp_Pnt(0, 0, 25))
    pt2 = np.loadtxt(obj.rootname + "_pln2.txt")
    print(pt2)
    pts = []
    for xyz in pt2 + [pt2[0]]:
        pts.append(gp_Pnt(*xyz))
    br2 = make_polygon(pts, closed=True)
    br2.Location(set_loc(gp_Ax3(), ax2))
    fc2 = obj.make_FaceByOrder(pts)
    fc2.Location(set_loc(gp_Ax3(), ax2))
    print(fc2)
    obj.display.DisplayShape(br2)
    obj.display.DisplayShape(fc2)

    api = BRepOffsetAPI_ThruSections()
    api.SetSmoothing(True)
    api.AddWire(br1)
    api.AddWire(br2)
    api.Build()
    shp = api.Shape()
    obj.display.DisplayShape(shp)
    #obj.export_stp(shp)
    obj.ShowOCC()
