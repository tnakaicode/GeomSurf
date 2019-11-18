import numpy as np
import matplotlib.pyplot as plt
import sys

from OCC.Display.SimpleGui import init_display
from OCC.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.gp import gp_Lin, gp_Sphere
from OCC.gp import gp_Mat, gp_XYZ
from OCC.gp import gp_Trsf, gp_GTrsf
from OCC.BRepFill import BRepFill_Filling
from OCC.GeomAbs import GeomAbs_G1, GeomAbs_G2
from OCC.GeomAbs import GeomAbs_C0, GeomAbs_C1, GeomAbs_C2, GeomAbs_C3
from OCC.Extend.DataExchange import write_step_file
from OCCUtils.Construct import make_edge

from base import plotocc


def fill_surface():
    n_sided = BRepFill_Filling()
    p0 = gp_Pnt(-20, -20, 0)
    p1 = gp_Pnt(+20, -20, 10)
    p2 = gp_Pnt(+20, +20, 0)
    p3 = gp_Pnt(-20, +20, 0)
    p4 = gp_Pnt(-10, -10, +5)
    p5 = gp_Pnt(-10, +10, -5)
    p6 = gp_Pnt(+10, -10, -10)
    p7 = gp_Pnt(+10, +10, +10)
    p8 = gp_Pnt(-15, -15, +2)
    p9 = gp_Pnt(-15, +15, -15)
    p10 = gp_Pnt(+15, -15, -2)
    p11 = gp_Pnt(+15, +15, +50)

    n_sided.Add(make_edge(p0, p1), GeomAbs_C0)
    n_sided.Add(make_edge(p1, p2), GeomAbs_C0)
    n_sided.Add(make_edge(p2, p3), GeomAbs_C0)
    n_sided.Add(make_edge(p3, p0), GeomAbs_C0)
    n_sided.Add(p4)
    n_sided.Add(p5)
    n_sided.Add(p6)
    n_sided.Add(p7)
    n_sided.Add(p8)
    n_sided.Add(p9)
    n_sided.Add(p10)
    n_sided.Add(p11)
    n_sided.Build()
    write_step_file(n_sided.Face(), "./FillSurf.stp")
    return n_sided.Face()


if __name__ == '__main__':
    obj = plotocc()
    obj.show_axs_pln(scale=20)
    obj.display.DisplayShape(fill_surface())
    obj.show()
