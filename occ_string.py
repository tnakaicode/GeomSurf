import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import argparse
from linecache import getline, clearcache

basename = os.path.dirname(__file__)

sys.path.append(os.path.join("../"))
from src.base_occ import dispocc

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Compound
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Core.BRepOffset import BRepOffset_MakeOffset, BRepOffset_Skin, BRepOffset_Interval
from OCC.Core.GeomAbs import GeomAbs_C0, GeomAbs_Intersection, GeomAbs_Arc
from OCC.Core.Addons import text_to_brep, register_font, Font_FA_Regular
from OCCUtils.Topology import Topo
from OCCUtils.Construct import make_box, make_line, make_wire, make_edge
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

    # register the font
    register_font("../font/arialbi.ttf")
    register_font("../font/msgothic.ttc")
    register_font("../font/arialuni.ttf")
    text = r"あ" + "\n" + r"い"
    text = r"あ"
    brep_string = text_to_brep(text, "MS Gothic",
                               Font_FA_Regular, 12.0, True)
    print(brep_string)
    api = BRepOffset_MakeOffset(
        brep_string, 1.0, 1.0E-5, BRepOffset_Skin,
        False, False, GeomAbs_Arc,
        True, True
    )
    # const TopoDS_Shape & 	S,
    # const Standard_Real 	Offset,
    # const Standard_Real 	Tol,
    # const BRepOffset_Mode 	Mode = BRepOffset_Skin,
    # const Standard_Boolean 	Intersection = Standard_False,
    # const Standard_Boolean 	SelfInter = Standard_False,
    # const GeomAbs_JoinType 	Join = GeomAbs_Arc,
    # const Standard_Boolean 	Thickening = Standard_False,
    # const Standard_Boolean 	RemoveIntEdges = Standard_False
    print(api.Error())
    # BRepOffset_NoError = 0
    # BRepOffset_UnknownError = 1
    # BRepOffset_BadNormalsOnGeometry = 2
    # BRepOffset_C0Geometry = 3
    # BRepOffset_NullOffset = 4
    # BRepOffset_NotConnectedShell = 5
    sold_string = api.Shape()

    obj.display.DisplayShape(brep_string)
    obj.display.DisplayShape(sold_string)

    # obj.show_axs_pln(scale=15)
    obj.ShowOCC()
