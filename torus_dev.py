import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import argparse

basename = os.path.dirname(__file__) + "/"

sys.path.append(os.path.join("./"))
from src.base_occ import dispocc
from src.base_widget import UVWidget

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

from PyQt5 import QtWidgets, QtGui, QtCore

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Mat, gp_XYZ
from OCC.Core.gp import gp_Lin, gp_Sphere
from OCC.Core.gp import gp_GTrsf, gp_Trsf
from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeShape, BRepBuilderAPI_MakeFace
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_GTransform
from OCC.Core.Geom import Geom_ToroidalSurface
from OCC.Core.GeomLib import GeomLib_Tool
from OCC.Core.ShapeAnalysis import ShapeAnalysis_Surface
from OCC.Core.GeomLProp import GeomLProp_SLProps, GeomLProp_SurfaceTool
from OCCUtils.Construct import make_n_sided, make_n_sections
from OCCUtils.Construct import make_edge, dir_to_vec, vec_to_dir




class Torus (dispocc, Geom_ToroidalSurface):

    def __init__(self, temp=True, disp=True, touch=False):
        super().__init__(temp, disp, touch)
        self.axs = gp_Ax3()
        self.radi = [500, 200]
        self.rxyz = [1.0, 1.0, 1.0]
        mat = gp_Mat(
            self.rxyz[0], 0, 0,
            0, self.rxyz[1], 0,
            0, 0, self.rxyz[2]
        )
        gtrf = gp_GTrsf(mat, gp_XYZ(0, 0, 0))
        self.t = Geom_ToroidalSurface(self.axs, *self.radi)
        self.face = BRepBuilderAPI_MakeFace(self.t, 1e-6).Face()
        self.face = BRepBuilderAPI_GTransform(self.face, gtrf).Shape()
        self.surf = BRep_Tool.Surface(self.face)
        self.prop = GeomLProp_SLProps(self.surf, 0.0, 0.0, 2, 0.1E-3)
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
        axs = gp_Ax3(p,vec_to_dir(vz), vec_to_dir(vx))
        print(u, v)
        print(p)
        print(vx)
        print(vy)
        print(vz)
        if self.prop.IsCurvatureDefined():
            d0, d1 = gp_Dir(0,0,1), gp_Dir(0,0,1)
            print("Gaussian :", 1 / self.prop.GaussianCurvature(), self.prop.GaussianCurvature())
            print("Max      :", 1 / self.prop.MaxCurvature(), self.prop.MaxCurvature())
            print("Min      :", 1 / self.prop.MinCurvature(), self.prop.MinCurvature())
            print("Mean     :", 1 / self.prop.MeanCurvature(), self.prop.MeanCurvature())
            self.prop.CurvatureDirections(d0, d1)
            print("Max:", dir_to_vec(d0))
            print("Min:", dir_to_vec(d1))
            print(d0.Dot(d1))
        else:
            print("NO Curvature defined")
        self.show_axs_pln(axs, scale=25)
        self.display.DisplayShape(p)
        self.display.DisplayVector(vz.Scaled(20), p, update=True)
        return p, vu, vv, vz

    def set_uv_widget(self):
        self.uvWidget = UVWidget(None)
        self.uvWidget.setWindowTitle("Torus UV Parameter")
        self.uvWidget.uval_labl.setText(f"U: 0~")
        self.uvWidget.vval_labl.setText(f"V: 0~")

        self.uvWidget.closeButton.clicked.connect(self.get_uv_val)
        self.uvWidget.closeButton.clicked.connect(self.uvWidget.close_self)
        self.uvWidget.show()

    def get_uv_val(self):
        u = self.text2float(self.uvWidget.uval_text.text(), 0.0)
        v = self.text2float(self.uvWidget.vval_text.text(), 0.0)
        self.get_prof([u, v])

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
    obj.add_menu("SPlop")
    obj.add_function("SPlop", obj.set_uv_widget)
    obj.get_prof([0, 0])
    obj.get_prof([0, 0.3])
    obj.get_prof([0.1, 0])
    obj.get_prof([0.1, 0.3])
    obj.ShowTorus()
