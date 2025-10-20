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
from OCC.Core.gp import gp_Lin
from OCC.Core.Geom import Geom_Line, Geom_Plane
from OCC.Core.GeomAPI import GeomAPI_IntCS, GeomAPI_IntSS, GeomAPI_ProjectPointOnSurf
from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepProj import BRepProj_Projection
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Compound
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCCUtils.Common import vertex2pnt
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

    obj = dispocc(disp=True, touch=True)

    axis = gp_Ax3()
    axis.Rotate(axis.Axis(), np.deg2rad(0))

    beam = gp_Ax3(gp_Pnt(1000, 2000, 3000),
                  gp_Dir(0, 0, 1))
    beam.SetDirection(vec_to_dir(gp_Vec(beam.Location(), axis.Location())))
    beam.SetYDirection(axis.Direction())
    beam_lin = make_edge(gp_Lin(beam.Axis()), -1, 2)
    obj.display.DisplayShape(beam_lin, color="BLACK")

    print("Beam")
    tor, pol = obj.calc_tor_angle1(axis, beam, False)
    tor, pol = obj.calc_tor_angle2(axis, beam, False)

    print("Beam-1")
    beam1 = beam.Rotated(gp_Ax1(beam.Location(),
                                beam.XDirection()),
                         pol + np.deg2rad(-5))
    obj.calc_tor_angle1(axis, beam1, False)
    obj.calc_tor_angle2(axis, beam1, True)

    print("Beam-2")
    beam2 = beam.Rotated(gp_Ax1(beam.Location(),
                                beam.XDirection()),
                         pol)
    beam2.Rotate(gp_Ax1(beam2.Location(),
                        beam2.YDirection()),
                 np.deg2rad(2))
    obj.calc_tor_angle1(axis, beam2, False)
    obj.calc_tor_angle2(axis, beam2, False)

    print("Beam-3")
    beam3 = beam.Rotated(gp_Ax1(beam.Location(),
                                beam.XDirection()),
                         pol)
    beam3.Rotate(gp_Ax1(beam2.Location(),
                        beam2.YDirection()),
                 np.deg2rad(5))
    beam3.Rotate(gp_Ax1(beam2.Location(),
                        beam2.XDirection()),
                 np.deg2rad(5))
    obj.calc_tor_angle1(axis, beam3, False)
    obj.calc_tor_angle2(axis, beam3, False)

    print("Beam-4")
    beam4 = obj.prop_axs(beam3, 100, "z")
    obj.calc_tor_angle1(axis, beam4, False)
    obj.calc_tor_angle2(axis, beam4, False)

    obj.ShowOCC()
