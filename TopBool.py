import numpy as np
import matplotlib.pyplot as plt
import sys

from OCC.Display.SimpleGui import init_display
from OCC.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.gp import gp_Lin, gp_Pln, gp_Sphere
from OCC.gp import gp_Mat, gp_XYZ
from OCC.gp import gp_Trsf, gp_GTrsf
from OCC.BRepPrimAPI import BRepPrimAPI_MakeWedge
from OCCUtils.Construct import make_n_sided, make_n_sections
from OCCUtils.Construct import make_edge

from base import plotocc




if __name__ == '__main__':
    obj = plotocc()
    obj.show_axs_pln(scale=20)

    Wedge = BRepPrimAPI_MakeWedge(60., 100., 80., 20.).Shape()
    obj.display.DisplayShape(Wedge)

    
    obj.show()
