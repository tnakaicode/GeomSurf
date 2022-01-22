import numpy as np
import matplotlib.pyplot as plt
import sys

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Lin, gp_Pln, gp_Sphere
from OCC.Core.gp import gp_Mat, gp_XYZ
from OCC.Core.gp import gp_Trsf, gp_GTrsf
from OCC.Core.TColgp import TColgp_Array1OfPnt, TColgp_Array2OfPnt
from OCC.Core.TColgp import TColgp_HArray1OfPnt, TColgp_HArray2OfPnt
from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepFill import BRepFill_Filling
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeSphere
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeEdge
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_GTransform
from OCC.Core.Geom import Geom_BSplineSurface
from OCC.Core.GeomAPI import geomapi
from OCC.Core.GeomAPI import GeomAPI_PointsToBSplineSurface
from OCC.Core.GeomAPI import GeomAPI_PointsToBSpline
from OCC.Core.GeomAPI import GeomAPI_Interpolate
from OCC.Core.GeomAPI import GeomAPI_IntCS, GeomAPI_IntSS
from OCC.Core.GeomAbs import GeomAbs_C0, GeomAbs_C1, GeomAbs_C2
from OCC.Core.GeomAbs import GeomAbs_G1, GeomAbs_G2
from OCC.Core.GeomLProp import GeomLProp_CLProps, GeomLProp_CurveTool
#from OCC.Core.GeomFill import GeomFill_BSplineCurves
#from OCC.Core.GeomFill import GeomFill_StretchStyle, GeomFill_CoonsStyle, GeomFill_CurvedStyle
from OCCUtils.Construct import make_n_sided, make_n_sections
from OCCUtils.Construct import make_edge

from src.base_occ import dispocc, spl_curv_pts


if __name__ == '__main__':
    obj = dispocc()
    obj.show_axs_pln(scale=1)

    pt = np.linspace(0, 1, 100)
    num = 2
    omg = 0.5
    idx = 2
    pts = []
    for t in pt[:-1]:
        ft = 1 - (idx - 1) * np.sin(2 * np.pi * t)
        if ft >= 1:
            ut = omg**ft
        else:
            ut = 1 - (1 - omg)**(2 - ft)
        x = np.cos(num * ut * 2 * np.pi)
        y = np.sin(num * ut * 2 * np.pi)
        #z = 2 * (t % 0.5)
        z = 2 * t
        pnt = gp_Pnt(x, y, z)
        pts.append(pnt)

    crv_p, crv = spl_curv_pts(pts)

    print(GeomLProp_CurveTool.FirstParameter(crv))
    print(GeomLProp_CurveTool.LastParameter(crv))

    obj.display.DisplayShape(crv)
    for pnt in pts[::10]:
        obj.display.DisplayShape(pnt)

    obj.show_occ()
