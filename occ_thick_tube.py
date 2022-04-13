import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import argparse

sys.path.append(os.path.join("../"))
from base_occ import dispocc, set_loc

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax3
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_ThruSections, BRepOffsetAPI_MakeOffset, BRepOffsetAPI_MakeEvolved, BRepOffsetAPI_MakePipe, BRepOffsetAPI_MakePipeShell
from OCC.Core.BRepOffset import BRepOffset_MakeOffset, BRepOffset_Skin, BRepOffset_Interval
from OCC.Core.GeomAbs import GeomAbs_C0
from OCC.Core.GeomAbs import GeomAbs_Intersection, GeomAbs_Arc
from OCCUtils.Construct import make_polygon
from OCCUtils.Construct import vec_to_dir, dir_to_vec

from OCC.Core.Quantity import (Quantity_Color, Quantity_TOC_RGB, Quantity_TOC_HLS,
                               Quantity_NOC_WHITE,
                               Quantity_NOC_BLACK, Quantity_NOC_BLUE1,
                               Quantity_NOC_CYAN1, Quantity_NOC_RED,
                               Quantity_NOC_GREEN, Quantity_NOC_ORANGE, Quantity_NOC_YELLOW)


def poly_rotate_shape(poly, axis=gp_Ax3(), num=50):
    api = BRepOffsetAPI_ThruSections()
    print(poly.Location().Transformation())
    print(dir_to_vec(axis.Direction()))
    for idx, phi in enumerate(np.linspace(0, 2 * np.pi, num)):
        ax = axis.Rotated(axis.Axis(), phi)
        poly_i = poly.Located(set_loc(axis, ax))
        api.AddWire(poly_i)
    api.Build()
    return api.Shape()


def get_color_from_name(color_name):
    ''' from the string 'WHITE', returns Quantity_Color
    WHITE.
    color_name is the color name, case insensitive.
    '''
    enum_name = 'Quantity_NOC_%s' % color_name.upper()
    if enum_name in globals():
        color_num = globals()[enum_name]
    elif enum_name + '1' in globals():
        color_num = globals()[enum_name + '1']
        print('Many colors for color name %s, using first.' % color_name)
    else:
        color_num = Quantity_NOC_WHITE
        print('Color name not defined. Use White by default')
    return Quantity_Color(color_num)


if __name__ == '__main__':
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument("--pxyz", dest="pxyz",
                      default=[0.0, 0.0, 0.0], type=float, nargs=3)
    opt = parser.parse_args()
    print(opt, argvs)

    px = np.linspace(-1, 1, 100) * 100 + 50
    py = np.linspace(-1, 1, 200) * 100 - 50
    mesh = np.meshgrid(px, py)

    p0 = gp_Pnt(0, 50.0, 0)
    p1 = gp_Pnt(0, 60.0, 100.0)
    p2 = gp_Pnt(0, 50.0, 200.0)
    axis = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1))
    poly = make_polygon([p0, p1, p2])
    face = poly_rotate_shape(poly, axis)
    sold = BRepOffset_MakeOffset(
        face, 1.0, 1.0E-5, BRepOffset_Skin,
        False, True, GeomAbs_Arc, True, True
    ).Shape()

    obj = dispocc(touch=True)
    obj.display.DisplayShape(poly)
    obj.display.DisplayShape(face)
    obj.display.DisplayShape(sold, transparency=0.9, color="BLUE1")
    #obj.display.DisplayShape(sold, transparency=0.9, color=255)
    obj.show_axs_pln(scale=50)
    obj.ShowOCC()

    # if color:
    # if isinstance(color, str):
    #     color = get_color_from_name(color)
    # elif isinstance(color, int):
    #     color = Quantity_Color(color)

    #
    # Creates a color according to the definition system theType. Quantity_TOC_RGB:
    #    theR1 the value of Red within range [0.0; 1.0]
    #    theR2 the value of Green within range [0.0; 1.0]
    #    theR3 the value of Blue within range [0.0; 1.0]
    #
    # Quantity_TOC_HLS:
    #    theR1 is the Hue (H) angle in degrees within range [0.0; 360.0], 0.0 being Red. Value -1.0 is a special value reserved for grayscale color (S should be 0.0).
    #    theR2 is the Lightness (L) within range [0.0; 1.0]
    #    theR3 is the Saturation (S) within range [0.0; 1.0]
    #
    # color should either be a string ( "BLUE" ) or a Quantity_Color(0.1, 0.8, 0.1) got
    #
    # yellow(255,255,0)
    # green(0,255,0)
    # cyan(0,255,255)
    # red(255,0,0)
    # blue(0,0,255)
    # red(255,0,0)
    # magenta(255,0,255)
