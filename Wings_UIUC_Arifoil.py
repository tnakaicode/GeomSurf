import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import urllib.request as urllib2  # Python3
from optparse import OptionParser

sys.path.append(os.path.join("./"))
from base import plot2d, plotocc, spl_face, set_loc

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.BOPAlgo import BOPAlgo_Builder, BOPAlgo_Splitter
from OCC.Core.BRepFeat import BRepFeat_MakeDPrism
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakePrism
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Compound
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Core.TopAbs import TopAbs_EDGE, TopAbs_FACE, TopAbs_SHAPE, TopAbs_SHELL, TopAbs_SOLID
from OCC.Core.TopExp import TopExp_Explorer
from OCCUtils.Topology import Topo
from OCCUtils.Construct import make_box, make_line, make_wire, make_edge
from OCCUtils.Construct import make_plane, make_polygon
from OCCUtils.Construct import point_to_vector, vector_to_point
from OCCUtils.Construct import dir_to_vec, vec_to_dir


def uiuc_database(name="dae51"):
    uiuc_url = 'http://m-selig.ae.illinois.edu/ads/coord_seligFmt/'
    foil_dat_url = uiuc_url + '{}.dat'.format(name)

    data = []
    fp = urllib2.urlopen(foil_dat_url)
    fp_lines = fp.readlines()
    for idx, line in enumerate(fp_lines[1:]):
        data.append([float(v) for v in line.split()])
    return np.array(data)


if __name__ == '__main__':
    argvs = sys.argv
    parser = OptionParser()
    parser.add_option("--dir", dest="dir", default="./")
    parser.add_option("--pxyz", dest="pxyz",
                      default=[0.0, 0.0, 0.0], type="float", nargs=3)
    opt, argc = parser.parse_args(argvs)
    print(opt, argc)

    obj = plotocc(touch=True)
    dae_data = uiuc_database()

    axs = gp_Ax3()
    uic = "dae51"
    pnt = gp_Pnt(0, 10, 0)
    ax1 = gp_Ax3(pnt, axs.YDirection(),
                 axs.XDirection().Reversed())
    ui1 = make_polygon([gp_Pnt(*100 * xy, 0)
                        for xy in uiuc_database(uic)], closed=True)
    ui1.Location(set_loc(gp_Ax3(), ax1))
    obj.display.DisplayShape(ui1)
    obj.display.DisplayMessage(pnt, uic)

    uic = "geminism"
    pnt = gp_Pnt(0, 20, 0)
    ax1.SetLocation(pnt)
    ui1 = make_polygon([gp_Pnt(*100 * xy, 0)
                        for xy in uiuc_database(uic)], closed=True)
    ui1.Location(set_loc(gp_Ax3(), ax1))
    obj.display.DisplayShape(ui1)
    obj.display.DisplayMessage(pnt, uic)

    uic = "naca0006"
    pnt = gp_Pnt(0, 30, 0)
    ax1.SetLocation(pnt)
    ui1 = make_polygon([gp_Pnt(*100 * xy, 0)
                        for xy in uiuc_database(uic)], closed=True)
    ui1.Location(set_loc(gp_Ax3(), ax1))
    obj.display.DisplayShape(ui1)
    obj.display.DisplayMessage(pnt, uic)

    uic = "naca001034a08cli0.2"
    pnt = gp_Pnt(0, 40, 0)
    ax1.SetLocation(pnt)
    ui1 = make_polygon([gp_Pnt(*100 * xy, 0)
                        for xy in uiuc_database(uic)], closed=True)
    ui1.Location(set_loc(gp_Ax3(), ax1))
    obj.display.DisplayShape(ui1)
    obj.display.DisplayMessage(pnt, uic)

    uic = "ncambre"
    pnt = gp_Pnt(0, 0, 0)
    ax1.SetLocation(pnt)
    ui1 = make_polygon([gp_Pnt(*100 * xy, 0)
                        for xy in uiuc_database(uic)], closed=True)
    ui1.Location(set_loc(gp_Ax3(), ax1))
    obj.display.DisplayShape(ui1)
    obj.display.DisplayMessage(pnt, uic)

    uic = "ah81k144wfKlappe"
    pnt = gp_Pnt(0, -10, 0)
    ax1.SetLocation(pnt)
    ui1 = make_polygon([gp_Pnt(*100 * xy, 0)
                        for xy in uiuc_database(uic)], closed=True)
    ui1.Location(set_loc(gp_Ax3(), ax1))
    obj.display.DisplayShape(ui1)
    obj.display.DisplayMessage(pnt, uic)

    uic = "eiffel430"
    pnt = gp_Pnt(0, -20, 0)
    ax1.SetLocation(pnt)
    ui1 = make_polygon([gp_Pnt(*100 * xy, 0)
                        for xy in uiuc_database(uic)], closed=True)
    ui1.Location(set_loc(gp_Ax3(), ax1))
    obj.display.DisplayShape(ui1)
    obj.display.DisplayMessage(pnt, uic)

    obj.show_axs_pln(axs, scale=10.0)
    obj.show()
