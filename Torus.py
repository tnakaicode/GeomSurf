import numpy as np
import matplotlib.pyplot as plt
import sys

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Lin, gp_Sphere
from OCC.Core.Geom import Geom_ToroidalSurface
from OCCUtils.Construct import make_n_sided, make_n_sections
from OCCUtils.Construct import make_edge

from base import plotocc


if __name__ == '__main__':
    obj = plotocc()

    t = Geom_ToroidalSurface(gp_Ax3(), 100, 50)

    obj.display.DisplayShape(t, transparency=0.7)
    obj.display.DisplayShape(t.Value(0, 0))
    obj.display.DisplayShape(t.Value(0, 2*np.pi*(1/3)))
    obj.display.DisplayShape(t.Value(0, 2*np.pi*(2/3)))
    obj.display.DisplayShape(t.Value(2*np.pi*(1/4), 0))
    obj.display.DisplayShape(t.Value(2*np.pi*(2/4), 0))
    obj.display.DisplayShape(t.Value(2*np.pi*(3/4), 0))
    obj.show_axs_pln(scale=100)
    obj.show()
