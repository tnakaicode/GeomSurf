import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import argparse

sys.path.append(os.path.join("./"))
from src.base_occ import dispocc

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Mat, gp_XYZ
from OCC.Core.gp import gp_Lin, gp_Sphere
from OCC.Core.gp import gp_GTrsf, gp_Trsf
from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeShape, BRepBuilderAPI_MakeFace
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_GTransform
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeSphere
from OCC.Core.Geom import Geom_ToroidalSurface, Geom_SphericalSurface
from OCC.Core.GeomLib import GeomLib_Tool
from OCC.Core.ShapeAnalysis import ShapeAnalysis_Surface
from OCC.Core.GeomLProp import GeomLProp_SLProps, GeomLProp_SurfaceTool
from OCCUtils.Construct import make_n_sided, make_n_sections
from OCCUtils.Construct import make_edge


class Torus (dispocc):

    def __init__(self):
        super().__init__()
        self.axs = gp_Ax3()
        self.radi = [500, 200]
        self.rxyz = [1.5, 1.2, 1.5]
        mat = gp_Mat(
            self.rxyz[0], 0, 0,
            0, self.rxyz[1], 0,
            0, 0, self.rxyz[2]
        )
        gtrf = gp_GTrsf(mat, gp_XYZ(0, 0, 0))
        #self.t = Geom_ToroidalSurface(self.axs, *self.radi)
        self.t = Geom_SphericalSurface(self.axs, 100.0)
        self.face = BRepBuilderAPI_MakeFace(self.t, 1e-6).Face()
        self.face = BRepBuilderAPI_GTransform(self.face, gtrf).Shape()
        self.surf = BRep_Tool.Surface(self.face)
        self.prop = GeomLProp_SLProps(self.surf, 0.0, 0.0, 1, 1.0)
        self.export_stp(self.face)
        print(self.t.UPeriod())

    def get_prof(self, uv=[0, 0]):
        u, v = uv
        u1, v1 = 2 * np.pi * u, 2 * np.pi * v
        p, vu, vv = gp_Pnt(), gp_Vec(), gp_Vec()
        self.surf.D1(u1, v1, p, vu, vv)
        self.prop.SetParameters(u1, v1)
        vx = vu.Normalized()
        vy = vv.Normalized()
        vz = vx.Crossed(vy)
        print(u, v)
        print(p)
        print(vx)
        print(vy)
        print(vz)
        #print(self.prop.GaussianCurvature())
        #print(self.prop.MaxCurvature())
        #print(self.prop.MinCurvature())
        #print(self.prop.MeanCurvature())
        self.display.DisplayShape(p)
        self.display.DisplayVector(vz.Scaled(20), p)
        return p, vu, vv, vz

    def ShowTorus(self):
        self.display.DisplayShape(self.face, transparency=0.7)
        self.show_axs_pln(scale=100)
        self.ShowOCC()


if __name__ == '__main__':
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    opt = parser.parse_args()
    print(opt, argvs)

    obj = Torus()
    obj.get_prof([0, 0])
    obj.get_prof([0.1, 0.3])
    obj.ShowTorus()
