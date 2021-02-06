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
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
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

from src.base import plotocc, spl_face


def trsf_scale(axs=gp_Ax3(), scale=1):
    trf = gp_Trsf()
    trf.SetDisplacement(gp_Ax3(), axs)
    return trf


if __name__ == '__main__':
    obj = plotocc(touch=True)
    obj.show_axs_pln(scale=5)

    px = np.linspace(-1, 1, 100) * 2 * np.pi * 2.5
    py = np.linspace(-1, 1, 150) * 2 * np.pi * 2.5
    mesh = np.meshgrid(px, py)
    surf1 = np.sin(mesh[0] / 4) * 0.5
    surf2 = np.sin(mesh[1] / 2)
    face1 = spl_face(*mesh, surf1)
    face2 = spl_face(*mesh, surf2)

    api = GeomAPI_IntSS(BRep_Tool.Surface(
        face1), BRep_Tool.Surface(face2), 1.0e-7)
    print(api.NbLines())

    obj.display.DisplayShape(face1)
    obj.display.DisplayShape(face2)
    for i in range(api.NbLines()):
        crv = api.Line(i + 1)
        p0 = GeomLProp_CurveTool.FirstParameter(crv)
        p1 = GeomLProp_CurveTool.LastParameter(crv)
        print(p0, p1)
        obj.display.DisplayShape(crv)

    obj.show()
