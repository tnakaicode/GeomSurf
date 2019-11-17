import numpy as np
import matplotlib.pyplot as plt
import sys

from OCC.Display.SimpleGui import init_display
from OCC.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.gp import gp_Lin, gp_Pln, gp_Sphere
from OCC.gp import gp_Mat, gp_XYZ
from OCC.gp import gp_Trsf, gp_GTrsf
from OCC.TColgp import TColgp_Array1OfPnt, TColgp_Array2OfPnt
from OCC.TColgp import TColgp_HArray1OfPnt, TColgp_HArray2OfPnt
from OCC.BRepFill import BRepFill_Filling
from OCC.BRepPrimAPI import BRepPrimAPI_MakeSphere
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.BRepBuilderAPI import BRepBuilderAPI_GTransform
from OCC.GeomAPI import geomapi
from OCC.GeomAPI import GeomAPI_PointsToBSplineSurface
from OCC.GeomAPI import GeomAPI_PointsToBSpline
from OCC.GeomAPI import GeomAPI_Interpolate
from OCC.GeomAbs import GeomAbs_C0, GeomAbs_C1, GeomAbs_C2
from OCC.GeomAbs import GeomAbs_G1, GeomAbs_G2
from OCC.GeomFill import GeomFill_BSplineCurves
from OCC.GeomFill import GeomFill_StretchStyle, GeomFill_CoonsStyle, GeomFill_CurvedStyle
from OCCUtils.Construct import make_n_sided, make_n_sections
from OCCUtils.Construct import make_edge

from base import plotocc, spl_face


def trsf_scale(axs=gp_Ax3(), scale=1):
    trf = gp_Trsf()
    trf.SetDisplacement(gp_Ax3(), axs)
    return trf





if __name__ == '__main__':
    obj = plotocc()
    obj.show_axs_pln(scale=5)

    px = np.linspace(-1, 1, 10) * 2*np.pi
    py = np.linspace(-1, 1, 15) * 2*np.pi
    mesh = np.meshgrid(px, py)
    surf = np.sin(mesh[0])
    
    obj.display.DisplayShape(spl_face(*mesh, np.sin(mesh[0])))
    obj.display.DisplayShape(spl_face(*mesh, np.sin(mesh[1])))
    
    obj.show()
