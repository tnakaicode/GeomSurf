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
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Compound, TopoDS_Solid, TopoDS_Shell
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Core.ChFi3d import ChFi3d_Rational
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

    box_trsf = gp_Trsf()
    box_trsf.SetRotation(
        gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), np.deg2rad(0))
    box = make_box(gp_Pnt(0, 0, 0), 20, 20, 20)
    box_faces = list(TopologyExplorer(box).faces())
    print(box)
    print(box_faces[0].HashCode(1000000))
    print([vertex2pnt(v) for v in TopologyExplorer(box_faces[0]).vertices()])
    box_faces[0].Move(TopLoc_Location(box_trsf), True)
    # box_faces[0].Reverse()
    print(box_faces[0].HashCode(1000000))
    print([vertex2pnt(v) for v in TopologyExplorer(box_faces[0]).vertices()])
    obj.selected_shape = [box_faces[0],
                          box_faces[2]]
    # obj.selected_shape = box_faces
    box = obj.make_shell_selcted()
    print(box)

    pts1 = [
        gp_Pnt(-20, 0, 0),
        gp_Pnt(30, 0, 0),
        gp_Pnt(30, 20, 0),
        gp_Pnt(-20, 30, 0),
    ]
    face1 = make_face(make_polygon(pts1, closed=True))
    edge1 = list(TopologyExplorer(face1).vertices())

    face2 = make_face(make_polygon(pts1, closed=True))
    trsf2 = gp_Trsf()
    trsf2.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0)), np.deg2rad(60))
    face2.Move(TopLoc_Location(trsf2))
    edge2 = list(TopologyExplorer(face2).vertices())

    face3 = make_face(make_polygon(
        [edge1[1], edge1[2], edge2[2]], closed=True))
    face4 = make_face(make_polygon(
        [edge1[0], edge1[3], edge2[3]], closed=True))

    obj.selected_shape = [face1, face2, face3, face4]
    face = obj.make_shell_selcted()
    faces = list(TopologyExplorer(face).faces())
    find_edge = LocOpe_FindEdges(box_faces[0], box_faces[2])
    find_edge.InitIterator()
    face_edge = find_edge.EdgeTo()
    print([vertex2pnt(v)
          for v in TopologyExplorer(find_edge.EdgeFrom()).vertices()])
    print([vertex2pnt(v)
          for v in TopologyExplorer(find_edge.EdgeTo()).vertices()])
    # face_edge = list(TopologyExplorer(box).edges())[0]

    pr = Message_ProgressRange()
    fillet = BRepFilletAPI_MakeFillet(box, 0)
    fillet.Add(face_edge)
    fillet.SetRadius(5, 1, 1)
    print(list(TopologyExplorer(box).faces_from_edge(face_edge)))
    fillet.Build()
    if fillet.IsDone():
        obj.display.DisplayShape(fillet.Shape())
    else:
        print(fillet.IsDone())

    print(list(takewhile(lambda x: x <= 20, accumulate(repeat(5), lambda x, _: x + 3))))

    # obj.display.DisplayShape([face1, face2])
    obj.display.DisplayShape(face_edge, color="BLUE1")
    obj.display.DisplayShape(box, transparency=0.7)
    obj.ShowOCC()
