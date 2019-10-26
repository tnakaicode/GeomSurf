import numpy as np
import matplotlib.pyplot as plt
import sys

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Lin, gp_Sphere
from OCC.Core.gp import gp_Mat, gp_GTrsf, gp_XYZ
from OCC.Core.gp import gp_Trsf, gp_GTrsf
from OCC.Core.TColgp import TColgp_Array2OfPnt
from OCC.Core.BRepFill import BRepFill_Filling
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeSphere
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_GTransform
from OCC.Core.GeomAPI import GeomAPI_PointsToBSplineSurface
from OCC.Core.GeomAPI import GeomAPI_Interpolate
from OCC.Core.GeomAbs import GeomAbs_C2, GeomAbs_C0, GeomAbs_G1, GeomAbs_G2
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
    # obj.show_ball()

    obj.display.DisplayShape(gen_ellipsoid(), transparency=0.0)

    axs = gp_Ax3(
        gp_Pnt(3, 0, 0), gp_Dir(1, 1, 1)
    )
    #obj.display.DisplayShape(gen_ellipsoid(axs), transparency=0.1)

    px = np.linspace(-1, 1, 10) * 100 / 2
    py = np.linspace(-1, 1, 15) * 120 / 2
    mesh = np.meshgrid(px, py)
    surf = mesh[0]**2 / 100 + mesh[1]**2 / 1000
    surf[7, 5] = 50
    obj.display.DisplayShape(spl_face(*mesh, surf))

    p0 = gp_Pnt(mesh[0][0, 0], mesh[1][0, 0], surf[0, 0])
    p1 = gp_Pnt(mesh[0][0, -1], mesh[1][0, -1], surf[0, -1])
    p2 = gp_Pnt(mesh[0][-1, 0], mesh[1][-1, 0], surf[-1, 0])
    p3 = gp_Pnt(mesh[0][-1, 1], mesh[1][-1, -1], surf[-1, -1])

    n_sided = BRepFill_Filling()
    n_sided.Add(make_edge(p0, p1), GeomAbs_C2)
    n_sided.Add(make_edge(p0, p2), GeomAbs_C2)
    n_sided.Add(make_edge(p1, p3), GeomAbs_C2)
    n_sided.Add(make_edge(p2, p3), GeomAbs_C2)
    n_sided.Add(gp_Pnt(0, 0, 50))

    """n_sided.Add(p0)
    n_sided.Add(p1)
    n_sided.Add(p2)
    n_sided.Add(p3)"""
    n_sided.Build()
    face = n_sided.Face()
    obj.display.DisplayShape(face)

    obj.show_axs_pln(axs, scale=20)
    obj.show()
