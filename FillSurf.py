import numpy as np
import matplotlib.pyplot as plt
import sys

from OCC.Display.SimpleGui import init_display
from OCC.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.gp import gp_Lin, gp_Sphere
from OCC.gp import gp_Mat, gp_XYZ
from OCC.gp import gp_Trsf, gp_GTrsf
from OCC.FEmTool import FEmTool_Curve, FEmTool_ProfileMatrix
from OCC.BRepFill import BRepFill_Filling, BRepFill_Draft, BRepFill_Pipe
from OCC.GeomAbs import GeomAbs_G1, GeomAbs_G2
from OCC.GeomAbs import GeomAbs_C0, GeomAbs_C1, GeomAbs_C2, GeomAbs_C3
from OCC.Extend.DataExchange import write_step_file
from OCCUtils.Construct import make_edge

from base import plotocc


def fill_surface():
    n_sided = BRepFill_Filling()

    n_sided.SetResolParam(4, 30, 5, True)
    # Sets the parameters used for resolution.
    # The default values of these parameters have been chosen for a good ratio quality/performance.
    #
    # Degree: it is the order of energy criterion to minimize for computing the deformation of the surface.
    #   The default value is 3.
    #   The recommanded value is i+2
    #   where i is the maximum order of the constraints.
    #
    # NbPtsOnCur: it is the average number of points for discretisation of the edges.
    #
    # NbIter: it is the maximum number of iterations of the process.
    #   For each iteration the number of discretisation points is increased.
    #
    # Anisotropie:

    n_sided.SetConstrParam()
    # Sets the values of Tolerances used to control the constraint.
    # Tol2d:
    # Tol3d:
    #   it is the maximum distance allowed between the support surface and the constraints
    # TolAng: it is the maximum angle allowed between the normal of the surface and the constraints
    # TolCurv: it is the maximum difference of curvature allowed between the surface and the constraint

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
