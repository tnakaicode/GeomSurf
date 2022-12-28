import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import argparse

basename = os.path.dirname(__file__) + "/"

sys.path.append(os.path.join("./"))
from src.base_occ import dispocc

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
from OCCUtils.Construct import make_edge


class NewWidget(QtWidgets.QWidget):

    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent)
        self.left = 200
        self.top = 200
        self.width = 125
        self.height = 125
        self.setGeometry(self.left, self.top, self.width, self.height)

        self.setting = QtCore.QSettings(basename + "temp_setting/Widget.ini",
                                        QtCore.QSettings.IniFormat)
        self.setting.setFallbacksEnabled(False)
        self.move(self.setting.value("pos", self.pos()))
        self.resize(self.setting.value("size", self.size()))
        font = self.font()
        font.setPointSize(self.setting.value("font", 9, int))
        self.setFont(font)

        self.closeButton = QtWidgets.QPushButton("Set & Close", self)
        self.uval_text = QtWidgets.QLineEdit("0.0", self)
        self.vval_text = QtWidgets.QLineEdit("0.0", self)

        self.layout = QtWidgets.QVBoxLayout()
        self.layout.addWidget(self.closeButton)
        self.layout.addWidget(self.uval_text)
        self.layout.addWidget(self.vval_text)
        self.setLayout(self.layout)

    def close_self(self):
        if self.setting != None:
            # Write window size and position to config file
            self.setting.setValue("size", self.size())
            self.setting.setValue("pos", self.pos())
            self.setting.setValue("font", self.font().pointSize())
        self.close()


class Torus (dispocc, Geom_ToroidalSurface):

    def __init__(self, temp=True, disp=True, touch=False):
        super().__init__(temp, disp, touch)
        self.axs = gp_Ax3()
        self.radi = [500, 200]
        self.rxyz = [1.0, 1.1, 1.1]
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
        print(u, v)
        print(p)
        print(vx)
        print(vy)
        print(vz)
        if self.prop.IsCurvatureDefined():
            print("Gaussian :", 1 / self.prop.GaussianCurvature())
            print("Max      :", 1 / self.prop.MaxCurvature())
            print("Min      :", 1 / self.prop.MinCurvature())
            print("Mean     :", 1 / self.prop.MeanCurvature())
        else:
            print("NO Curvature defined")
        self.display.DisplayShape(p)
        self.display.DisplayVector(vz.Scaled(20), p)
        return p, vu, vv, vz

    def set_uv_widget(self):
        self.uvWidget = NewWidget(None)

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
