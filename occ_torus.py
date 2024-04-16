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
from OCC.Core.gp import gp_Trsf
from OCC.Core.gp import gp_Circ
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_MakePipe
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

    obj = dispocc(touch=True)
    axs = gp_Ax3()
    tok = gp_Ax3(gp_Pnt(0, 0, 0.0), gp_Dir(0, 1, 0), gp_Dir(-1, 0, 0))

    z = 0
    major = 2000
    minor = 5000
    radii = 3000

    pnt_major = gp_Pnt(0, 0, z)
    axs_major = gp_Ax2(pnt_major, gp_Dir(0, 0, 1), gp_Dir(1, 0, 0))
    rng_major = [0 + np.pi / 3, np.pi - np.pi / 3]

    pnt_minor = gp_Pnt(0, radii - minor, z)
    axs_minor = gp_Ax2(pnt_minor, gp_Dir(1, 0, 0), gp_Dir(0, 1, 0))
    rng_minor = [-np.pi / 2 + np.pi / 3, np.pi / 2 - np.pi / 3]

    cir_major = make_edge(gp_Circ(axs_major, major), *rng_major)
    cir_minor = make_edge(gp_Circ(axs_minor, minor), *rng_minor)
    api = BRepOffsetAPI_MakePipe(make_wire(cir_major), cir_minor)
    api.Build()
    shll = api.Shape()
    face = [f for f in TopologyExplorer(shll).faces()][0]
    surf = BRepAdaptor_Surface(face)

    trsf = gp_Trsf()
    trsf.SetDisplacement(gp_Ax3(), tok)
    nu, nv = 3, 4
    us, ue = surf.FirstUParameter(), surf.LastUParameter()
    vs, ve = surf.FirstVParameter(), surf.LastVParameter()
    pu = np.linspace(us, ue, nu)
    pv = np.linspace(vs, ve, nv)
    mesh = np.meshgrid(pu, pv)
    data = np.zeros_like(mesh[0])
    for i, u in enumerate(pu):
        for j, v in enumerate(pv):
            p0 = surf.Value(u, v)
            obj.display.DisplayShape(p0)
            p0.Transform(trsf)
            print(p0)
            # mesh[0][j, i] = p0.X()
            # mesh[1][j, i] = p0.Y()
            data[j, i] = p0.Z()

    obj.show_axs_pln(tok)
    obj.display.DisplayShape(shll)
    obj.display.DisplayShape(cir_minor)
    obj.display.DisplayShape(cir_major)
    obj.ShowOCC()
