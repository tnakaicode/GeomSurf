import numpy as np
import matplotlib.pyplot as plt
import sys

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.TColgp import TColgp_Array2OfPnt
from OCC.Core.GeomAPI import GeomAPI_PointsToBSplineSurface
from OCC.Core.GeomAbs import GeomAbs_C0, GeomAbs_C1, GeomAbs_C2
from OCC.Core.GeomAbs import GeomAbs_G1, GeomAbs_G2
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakePolygon
from OCC.Core.BRepProj import BRepProj_Projection


def spl_face(px, py, pz):
    nx, ny = px.shape
    pnt_2d = TColgp_Array2OfPnt(1, nx, 1, ny)
    for row in range(pnt_2d.LowerRow(), pnt_2d.UpperRow() + 1):
        for col in range(pnt_2d.LowerCol(), pnt_2d.UpperCol() + 1):
            i, j = row - 1, col - 1
            pnt = gp_Pnt(px[i, j], py[i, j], pz[i, j])
            pnt_2d.SetValue(row, col, pnt)
            #print (i, j, px[i, j], py[i, j], pz[i, j])

    curv = GeomAPI_PointsToBSplineSurface(
        pnt_2d, 3, 8, GeomAbs_G2, 0.001).Surface()
    surf = BRepBuilderAPI_MakeFace(curv, 1e-6).Face()
    return surf


def make_polygon(args, closed=False):
    poly = BRepBuilderAPI_MakePolygon()
    for pt in args:
        # support nested lists
        if isinstance(pt, list) or isinstance(pt, tuple):
            for i in pt:
                poly.Add(i)
        else:
            poly.Add(pt)
    if closed:
        poly.Close()
    poly.Build()

    result = poly.Wire()
    return result


if __name__ == '__main__':
    display, start_display, add_menu, add_function = init_display()

    axs = gp_Ax3()

    px = np.linspace(-1, 1, 100) * 600
    py = np.linspace(-1, 1, 100) * 600
    mesh = np.meshgrid(px, py)
    surf = mesh[0]**2 / 1000 + mesh[1]**2 / 2000
    face = spl_face(*mesh, surf)

    pts = []
    pts.append(gp_Pnt(-100, -200, 0))
    pts.append(gp_Pnt(+200, -300, 0))
    pts.append(gp_Pnt(+300, +400, 0))
    pts.append(gp_Pnt(-400, +500, 0))
    poly = make_polygon(pts, closed=True)

    proj = BRepProj_Projection(poly, face, axs.Direction())
    proj_poly = proj.Current()

    display.DisplayShape(face, color="RED", transparency=0.9)
    display.DisplayShape(poly)
    display.DisplayShape(proj_poly, color="BLUE1")
    display.FitAll()
    display.View.Dump("Proj_Wire_to_Surf.png")
    start_display()
