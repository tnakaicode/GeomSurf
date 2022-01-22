import numpy as np
import matplotlib.pyplot as plt
import sys

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Lin, gp_Pln, gp_Sphere
from OCC.Core.gp import gp_Mat, gp_XYZ
from OCC.Core.gp import gp_Trsf, gp_GTrsf
from OCC.Core.BRepFill import BRepFill_Filling
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeWedge
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakePrism
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Common
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Fuse
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Section
from OCC.Core.BRepSweep import BRepSweep_Prism
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_MakePipe
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_ThruSections
from OCCUtils.Construct import make_n_sided, make_n_sections
from OCCUtils.Construct import make_edge

from src.base import dispocc, gen_ellipsoid


if __name__ == '__main__':
    obj = dispocc()
    obj.show_axs_pln(scale=20)

    elli1 = gen_ellipsoid(axs=gp_Ax3(gp_Pnt(-90, 0, 0),
                                     gp_Dir(0, 0, 1)), rxyz=[100, 100, 100])
    Wedge = BRepPrimAPI_MakeWedge(60., 100., 80., 20.).Shape()

    obj.display.DisplayShape(elli1)
    obj.display.DisplayShape(Wedge)

    obj.show_occ()
