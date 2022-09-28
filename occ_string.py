import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import argparse
from linecache import getline, clearcache

basename = os.path.dirname(__file__)

sys.path.append(os.path.join("./"))
from src.base_occ import dispocc
from src.base_occ import set_loc

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


def boolean_cut(shapeToCutFrom, cuttingShape):
    from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Cut

    cut = BRepAlgoAPI_Cut(shapeToCutFrom, cuttingShape)
    print(cut.Check())

    shp = cut.Shape()
    return shp

def boolean_common(shape1, shape2):
    from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Common
    from OCC.Core.BOPAlgo import BOPAlgo_PaveFiller

    com = BRepAlgoAPI_Common(shape1, shape2)
    print(com.Check())

    shp = com.Shape()
    return shp

def boolean_cut_bycommon(shapeToCutFrom, cuttingShape):
    from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Cut

    com = boolean_common(shapeToCutFrom, cuttingShape)
    cut = BRepAlgoAPI_Cut(shapeToCutFrom, com)
    print(cut.Check())

    shp = cut.Shape()
    return shp

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
        gp_Pnt(0, -10, 0),
        gp_Pnt(100, -10, 0),
        gp_Pnt(100, 20, 0),
        gp_Pnt(0, 20, 0)
    ]
    pts_axs = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1))
    plate = obj.make_plate(pts, skin=-5.0, axs=pts_axs)
    # obj.display.DisplayShape(plate)

    surf = obj.make_trimmedcylinder(axs, 100, 50, [-np.pi / 6, np.pi / 6])
    api = BRepOffset_MakeOffset(
        surf, 5.0, 1.0E-5, BRepOffset_Skin,
        False, False, GeomAbs_Arc,
        True, True
    )
    surf_soild = api.Shape()
    # obj.display.DisplayShape(surf_soild)

    # register the font
    register_font("../font/arialbi.ttf")
    register_font("../font/msgothic.ttc")
    register_font("../font/arialuni.ttf")
    text_axs = gp_Ax3(gp_Pnt(0, 0, 1), gp_Dir(0, 0, 1))
    text = r"あ" + "\n" + r"い"
    text = r"あい"
    text = "Aib"
    brep_string = text_to_brep(text, "MS Gothic",
                               Font_FA_Regular, 12.0, True)
    brep_string.Location(set_loc(axs, text_axs))
    print(brep_string)
    api = BRepOffset_MakeOffset(
        brep_string, -2.0, 1.0E-5, BRepOffset_Skin,
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
    if api.Error() == 5:
        obj.selected_shape = []
        for text_face in Topo(brep_string).faces():
            api = BRepOffset_MakeOffset(
                text_face, -2.0, 1.0E-5, BRepOffset_Skin,
                False, False, GeomAbs_Arc,
                True, True
            )
            obj.selected_shape.append(api.Shape())
        sold_string = obj.make_comp_selcted()
        obj.selected_shape = []
    else:
        sold_string = api.Shape()

    # obj.display.DisplayShape(brep_string)
    #obj.display.DisplayShape(sold_string)

    plate_counterbore = boolean_cut(plate, sold_string)
    # obj.display.DisplayShape(plate_counterbore)

    surf_counterbore = boolean_common(surf_soild, sold_string)
    obj.display.DisplayShape(surf_counterbore)

    # obj.show_axs_pln(scale=15)
    obj.ShowOCC()
