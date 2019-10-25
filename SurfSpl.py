import numpy as np
import matplotlib.pyplot as plt
import sys

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Lin, gp_Sphere
from OCC.Core.gp import gp_Mat, gp_GTrsf, gp_XYZ
from OCC.Core.gp import gp_Trsf, gp_GTrsf
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeSphere
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_GTransform

from base import plotocc


def trsf_scale(axs=gp_Ax3(), scale=1):
    trf = gp_Trsf()
    trf.SetDisplacement(gp_Ax3(), axs)
    return trf


def gen_ellipsoid(axs=gp_Ax3(), rxyz=[10, 20, 30]):
    mat = gp_Mat(
        rxyz[0], 0, 0,
        0, rxyz[1], 0,
        0, 0, rxyz[2]
    )
    gtrf = gp_GTrsf(mat, gp_XYZ(0, 0, 0))
    sphere = BRepPrimAPI_MakeSphere(axs.Ax2(), 1).Solid()
    ellips = BRepBuilderAPI_GTransform(sphere, gtrf).Shape()
    return ellips


if __name__ == '__main__':
    obj = plotocc()
    obj.show_axs_pln(scale=20)
    # obj.show_ball()

    obj.display.DisplayShape(gen_ellipsoid(), transparency=0.0)

    axs = gp_Ax3(
        gp_Pnt(3, 0, 0), gp_Dir(1, 1, 1)
    )
    obj.display.DisplayShape(gen_ellipsoid(axs), transparency=0.1)
    obj.show_axs_pln(axs, scale=20)
    obj.show()
