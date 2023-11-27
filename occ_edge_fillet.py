import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import argparse
from itertools import accumulate, repeat, takewhile
from linecache import getline, clearcache, updatecache

basename = os.path.dirname(__file__)

sys.path.append(os.path.join("./"))
from src.base_occ import dispocc

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3, gp_Trsf, gp_Translation
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Compound, TopoDS_Solid, TopoDS_Shell, TopoDS_Face
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Core.ChFi3d import ChFi3d_Rational
from OCC.Core.BRep import BRep_Builder
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_Sewing
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Fuse
from OCC.Core.BRepFill import BRepFill_Filling
from OCC.Core.BRepFilletAPI import BRepFilletAPI_MakeFillet, BRepFilletAPI_MakeChamfer
from OCC.Core.LocOpe import LocOpe_FindEdges
from OCC.Core.Message import Message_ProgressRange
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCCUtils.Common import vertex2pnt
from OCCUtils.Construct import make_box, make_line, make_wire, make_edge
from OCCUtils.Construct import make_plane, make_polygon, make_face
from OCCUtils.Construct import point_to_vector, vector_to_point
from OCCUtils.Construct import dir_to_vec, vec_to_dir


def find_fillet_edge(sewed=TopoDS_Shape(), face1=TopoDS_Face(), face2=TopoDS_Face()):
    # Find Edge that the two faces share
    find_edge = LocOpe_FindEdges(face1, face2)
    find_edge.InitIterator()
    face_edge = find_edge.EdgeTo()

    # Find Edge of sewed shape that is shared with common_edge
    find_edge = LocOpe_FindEdges(face_edge, sewed)
    find_edge.InitIterator()
    face_edge = find_edge.EdgeTo()
    print(list(TopologyExplorer(sewed).faces_from_edge(face_edge)))
    return face_edge


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

    pts1 = [
        gp_Pnt(-20, 0, 0),
        gp_Pnt(30, 0, 0),
        gp_Pnt(30, 20, 0),
        gp_Pnt(-20, 30, 0),
    ]
    face1 = make_face(make_polygon(pts1, closed=True))
    vert1 = list(TopologyExplorer(face1).vertices())

    face2 = make_face(make_polygon(pts1, closed=True))
    trsf2 = gp_Trsf()
    trsf2.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0)), np.deg2rad(60))
    face2.Move(TopLoc_Location(trsf2))
    vert2 = list(TopologyExplorer(face2).vertices())

    face3 = make_face(make_polygon(
        [vert1[1], gp_Pnt(30, 20, 0), gp_Pnt(31, 20, 10)], closed=True))
    face4 = make_face(make_polygon(
        [vert1[0], vert1[3], vert2[3]], closed=True))

    # Make Shell by only two faces that are sewed
    sew = BRepBuilderAPI_Sewing()
    for face in [face1, face2, face3, face4]:
        sew.Add(face)
    sew.Perform()
    sewed = sew.SewedShape()

    fillet_edge12 = find_fillet_edge(sewed, face1, face2)
    fillet_edge14 = find_fillet_edge(sewed, face1, face4)

    pr = Message_ProgressRange()
    fillet = BRepFilletAPI_MakeFillet(sewed, 0)
    fillet.Add(5, fillet_edge12)
    fillet.Add(3, fillet_edge14)
    fillet.Build()
    if fillet.IsDone():
        obj.display.DisplayShape(fillet.Shape())
        # obj.export_stp(fillet.Shape())
    else:
        print(fillet.IsDone())

    print(list(takewhile(lambda x: x <= 20, accumulate(repeat(5), lambda x, _: x + 3))))

    # obj.display.DisplayShape([face1, face2])
    obj.display.DisplayShape(fillet_edge12, color="BLUE1")
    obj.display.DisplayShape(fillet_edge14, color="BLUE1")
    obj.display.DisplayShape(sewed, transparency=0.7)
    obj.ShowOCC()
