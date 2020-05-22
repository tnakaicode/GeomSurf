import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
from optparse import OptionParser

sys.path.append(os.path.join("./"))
from base import plot2d, plotocc

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Lin, gp_Sphere
from OCC.Core.Geom import Geom_ToroidalSurface
from OCCUtils.Construct import make_n_sided, make_n_sections
from OCCUtils.Construct import make_edge


class Torus (plotocc):

    def __init__(self):
        super().__init__()
        self.t = Geom_ToroidalSurface(gp_Ax3(), 500, 100)
        print(self.t.UPeriod())

    def get_prof(self, uv=[0, 0]):
        u, v = uv
        p, vu, vv = gp_Pnt(), gp_Vec(), gp_Vec()
        self.t.D1(2 * np.pi * u, 2 * np.pi * v, p, vu, vv)
        vx = vu.Normalized()
        vy = vv.Normalized()
        vz = vx.Crossed(vy)
        print(u, v)
        print(p)
        print(vx)
        print(vy)
        print(vz)
        self.display.DisplayShape(p)
        self.display.DisplayVector(vz.Scaled(20), p)
        return p, vu, vv, vz

    def ShowTorus(self):
        self.display.DisplayShape(self.t, transparency=0.7)
        self.show()


if __name__ == '__main__':
    argvs = sys.argv
    parser = OptionParser()
    parser.add_option("--dir", dest="dir", default="./")
    opt, argc = parser.parse_args(argvs)
    print(opt, argc)

    obj = Torus()
    obj.get_prof([0, 0])
    obj.get_prof([0.1, 0.3])
    obj.ShowTorus()
