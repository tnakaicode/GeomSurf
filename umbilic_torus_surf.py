import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import argparse

from OCC.Core.gp import gp_Pnt, gp_PntMirror, gp_Vec, gp_Dir, gp_Pln
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.GeomAbs import GeomAbs_Shape
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Section
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_ThruSections
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCCUtils.Construct import make_polygon

sys.path.append(os.path.join("../"))
from base_occ import dispocc

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

if __name__ == '__main__':
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument("--pxyz", dest="pxyz",
                        default=[0.0, 0.0, 0.0], type=float, nargs=3)
    opt = parser.parse_args()
    print(opt, argvs)

    u1 = np.linspace(-1, 1, 100) * np.pi
    v1 = np.linspace(-1, 1, 100) * np.pi
    u, v = np.meshgrid(u1, v1)

    x = np.sin(u) * (7 + np.cos(u / 3 - 2 * v) + 2 * np.cos(u / 3 + v))
    y = np.cos(u) * (7 + np.cos(u / 3 - 2 * v) + 2 * np.cos(u / 3 + v))
    z = np.sin(u / 3 - 2 * v) + 2 * np.sin(u / 3 + v)

    obj = dispocc(touch=True)
    api = BRepOffsetAPI_ThruSections()
    for i0, v0 in enumerate(v1):
        pts = []
        for j0, u0 in enumerate(u1):
            pts.append(gp_Pnt(x[i0, j0], y[i0, j0], z[i0, j0]))
        wire = make_polygon(pts)
        api.AddWire(wire)
    api.SetContinuity(3)
    api.Build()
    shp = api.Shape()
    print(api.Continuity())
    print(list(TopologyExplorer(shp).faces()),
          TopologyExplorer(shp).number_of_faces())
    obj.display.DisplayShape(api.Shape(), transparency=0.9)

    axs = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0.1, 0, 1))
    pln = gp_Pln(axs)
    section = BRepAlgoAPI_Section(shp, pln)
    section.Build()
    print(section.SectionEdges().Size())
    obj.display.DisplayShape(section.Shape())

    axs = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0.5))
    pln = gp_Pln(axs)
    section = BRepAlgoAPI_Section(shp, pln)
    section.Build()
    print(section.SectionEdges().Size())
    obj.display.DisplayShape(section.Shape(), color="BLUE1")

    obj.ShowOCC()
