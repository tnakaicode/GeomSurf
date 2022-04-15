import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import argparse
from linecache import getline, clearcache

sys.path.append(os.path.join("./"))
from src.base_occ import dispocc, spl_curv_pts

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_XOY, gp_Pnt2d, gp_Dir2d, gp_Lin2d
from OCC.Core.BRepProj import BRepProj_Projection
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge
from OCC.Core.Geom import Geom_CylindricalSurface
from OCC.Core.GeomAPI import GeomAPI_PointsToBSplineSurface
from OCC.Core.GeomAPI import GeomAPI_PointsToBSpline
from OCC.Core.gce import gce_MakeTranslation
from OCC.Core.GCE2d import GCE2d_MakeSegment
from OCC.Extend.DataExchange import write_step_file, read_step_file
from OCC.Extend.DataExchange import write_stl_file, read_stl_file
from OCCUtils.Construct import make_face, make_polygon, make_wire, make_edge
from OCCUtils.Construct import dir_to_vec, vec_to_dir

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)


def make_spiral(r=1, z=1.0):
    pn = int(100 * r)
    pt = np.linspace(0, np.pi, pn)
    pz = np.linspace(0, z, pn)
    pts = []
    for i, t in enumerate(pt):
        x = np.cos(r * t)
        y = np.sin(r * t)
        z = pz[i]
        pts.append(gp_Pnt(x, y, z))
    curv = make_edge(spl_curv_pts(pts))
    edge = make_edge(pts[0], pts[-1])
    return make_wire(curv, edge)


def make_helix(r=1, z=1.0):
    # Build an helix
    aCylinder = Geom_CylindricalSurface(gp_Ax3(), r)
    aLine2d = gp_Lin2d(gp_Pnt2d(0.0, 0.0), gp_Dir2d(z, z))
    aSegment = GCE2d_MakeSegment(aLine2d, 0.0, np.pi * 2.0)

    return BRepBuilderAPI_MakeEdge(aSegment.Value(), aCylinder, 0.0, 2 * np.pi).Edge()


if __name__ == '__main__':
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument("--pxyz", dest="pxyz",
                        default=[0.0, 0.0, 0.0], type=float, nargs=3)
    opt = parser.parse_args()
    print(opt, argvs)

    obj = dispocc(touch=True)
    axs = gp_Ax3()
    surf = obj.make_cylinder_surf(
        axs, radii=1, hight=2.1, rng=[0, 2 * np.pi], xyz="z")
    sprl = make_spiral(r=2)
    helx = make_helix(r=2)
    proj = BRepProj_Projection(sprl, surf, axs.Location())
    sprl_proj = proj.Current()
    trim = make_face(surf, sprl_proj)

    # https://dev.opencascade.org/doc/occt-7.5.0/refman/html/class_b_rep_builder_a_p_i___make_edge.html#details
    #
    #proj = BRepProj_Projection(sprl, surf, axs.Location())
    #i = 0
    # while proj.More():
    #    shpe = proj.Current()
    #    obj.display.DisplayShape(shpe, color=obj.colors[i % 5])
    #    proj.Next()
    #    i += 1
    #    print(i)

    # obj.display.DisplayShape(surf)
    obj.display.DisplayShape(sprl, color="BLUE1")
    obj.display.DisplayShape(sprl_proj, color="GREEN")
    obj.display.DisplayShape(helx, color="RED")
    obj.display.DisplayShape(trim, transparency=0.9)

    obj.ShowOCC()
