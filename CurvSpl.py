import numpy as np
import matplotlib.pyplot as plt
import sys

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Lin, gp_Sphere
from OCC.Core.gp import gp_Mat, gp_GTrsf, gp_XYZ
from OCC.Core.gp import gp_Trsf, gp_GTrsf
from OCC.Core.TColgp import TColgp_Array1OfPnt, TColgp_Array2OfPnt
from OCC.Core.TColgp import TColgp_HArray1OfPnt, TColgp_HArray2OfPnt
from OCC.Core.BRepFill import BRepFill_Filling
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeSphere
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_GTransform
from OCC.Core.GeomAPI import geomapi
from OCC.Core.GeomAPI import GeomAPI_PointsToBSplineSurface
from OCC.Core.GeomAPI import GeomAPI_PointsToBSpline
from OCC.Core.GeomAPI import GeomAPI_Interpolate
from OCC.Core.GeomAbs import GeomAbs_C0, GeomAbs_C1, GeomAbs_C2
from OCC.Core.GeomAbs import GeomAbs_G1, GeomAbs_G2
from OCCUtils.Construct import make_n_sided, make_n_sections
from OCCUtils.Construct import make_edge

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


def spl_face(px, py, pz):
    nx, ny = px.shape
    pnt_2d = TColgp_Array2OfPnt(1, nx, 1, ny)
    for row in range(pnt_2d.LowerRow(), pnt_2d.UpperRow() + 1):
        for col in range(pnt_2d.LowerCol(), pnt_2d.UpperCol() + 1):
            i, j = row - 1, col - 1
            pnt = gp_Pnt(px[i, j], py[i, j], pz[i, j])
            pnt_2d.SetValue(row, col, pnt)
            #print (i, j, px[i, j], py[i, j], pz[i, j])

    api = GeomAPI_PointsToBSplineSurface(pnt_2d, 3, 8, GeomAbs_G2, 0.001)
    api.Interpolate(pnt_2d)
    #surface = BRepBuilderAPI_MakeFace(curve, 1e-6)
    # return surface.Face()
    return BRepBuilderAPI_MakeFace(api.Surface(), 1e-6).Face()


if __name__ == '__main__':
    obj = plotocc()
    obj.show_axs_pln(scale=20)

    px = np.linspace(-1, 1, 10) * 100 / 2
    py = np.linspace(-1, 1, 15) * 120 / 2
    mesh = np.meshgrid(px, py)
    surf = mesh[0]**2 / 100 + mesh[1]**2 / 1000
    surf[7, 5] = 50
    
    p0 = gp_Pnt(mesh[0][0, 0], mesh[1][0, 0], surf[0, 0])
    p1 = gp_Pnt(mesh[0][0, -1], mesh[1][0, -1], surf[0, -1])
    p2 = gp_Pnt(mesh[0][-1, 0], mesh[1][-1, 0], surf[-1, 0])
    p3 = gp_Pnt(mesh[0][-1, 1], mesh[1][-1, -1], surf[-1, -1])

    pt = np.linspace(0, 2 * np.pi, 100)
    pts = []
    p_array = TColgp_Array1OfPnt(1, 100)
    for idx, t in enumerate(pt):
        x = 30.0 * np.cos(t)
        y = 30.0 * np.sin(t)
        z = 15 * np.cos(3 * t)
        pnt = gp_Pnt(x, y, z)
        pts.append(pnt)
        p_array.SetValue(idx + 1, pnt)
    api = GeomAPI_PointsToBSpline(p_array)
    #curve = geomapi.To3d(api.Curve())

    for idx, t in enumerate(pts):
        obj.display.DisplayShape(t)
    obj.display.DisplayShape(api.Curve())

    obj.show()
