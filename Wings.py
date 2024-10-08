import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import urllib.request as urllib2  # Python3
import argparse

sys.path.append(os.path.join("./"))
from src.base_occ import dispocc, spl_face, set_loc

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

def uiuc_dae_data (name="dae51"):
    foil_dat_url = 'http://m-selig.ae.illinois.edu/ads/coord_seligFmt/{}.dat'.format(name)
    
    data = []
    fp = urllib2.urlopen(foil_dat_url)
    fp_lines = fp.readlines()
    for idx, line in enumerate (fp_lines[1:]):
        print(line.split())
        data.append([float(v) for v in line.split()])
    return np.array(data)


if __name__ == '__main__':
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument("--pxyz", dest="pxyz",
                      default=[0.0, 0.0, 0.0], type=float, nargs=3)
    opt = parser.parse_args()
    print(opt, argvs)

    px = np.linspace(-1, 1, 100) * 110 / 2
    py = np.linspace(-1, 1, 200) * 110 / 2
    mesh = np.meshgrid(px, py)
    dae_data = uiuc_dae_data()

    obj = dispocc()
    axs = gp_Ax3()
    ax1 = gp_Ax3(gp_Pnt(0, 0, -10), axs.Direction())
    vec = gp_Vec(gp_Pnt(0, 0, -10), gp_Pnt(0, 0, 10))
    face = obj.make_EllipWire(rxy=[50.0, 50.0], axs=ax1, skin=0)
    body = BRepPrimAPI_MakePrism(face, vec).Shape()
    data1 = (mesh[0]**2 / 750 + mesh[1]**2 / 750) + 6.0
    data2 = (mesh[0]**2 / 1000 + mesh[1]**2 / 1000) - 6.0
    face1 = spl_face(*mesh, data1, axs=axs)
    face2 = spl_face(*mesh, data2, axs=axs)

    splitter = BOPAlgo_Splitter()
    splitter.AddArgument(body)
    splitter.AddTool(face1)
    splitter.AddTool(face2)
    splitter.Perform()
    print(splitter.Arguments())
    print(splitter.ShapesSD())
    exp = TopExp_Explorer(splitter.Shape(), TopAbs_SOLID)
    shp = []
    while exp.More():
        print(exp.Current())
        shp.append(exp.Current())
        exp.Next()

    obj.show_axs_pln(axs, scale=20)
    obj.display.DisplayShape(face1, transparency=0.9, color="RED")
    obj.display.DisplayShape(face2, transparency=0.9, color="BLUE1")
    obj.display.DisplayShape(shp[1])
    # obj.display.DisplayShape(splitter.Shape())
    # obj.display.DisplayShape(sol)
    obj.ShowOCC()
