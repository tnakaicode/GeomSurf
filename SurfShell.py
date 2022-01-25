import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
from linecache import getline, clearcache
from optparse import OptionParser

sys.path.append(os.path.join("./"))
from src.base_occ import dispocc, set_trf

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Lin, gp_Sphere
from OCC.Core.gp import gp_Mat, gp_XYZ
from OCC.Core.gp import gp_Trsf, gp_GTrsf
from OCC.Core.TColgp import TColgp_Array1OfPnt, TColgp_Array2OfPnt
from OCC.Core.Geom import Geom_BezierSurface, Geom_BSplineSurface, Geom_RectangularTrimmedSurface
from OCC.Core.GeomAPI import GeomAPI_IntCS, GeomAPI_IntSS
from OCC.Core.GeomAPI import GeomAPI_PointsToBSplineSurface, GeomAPI_ProjectPointOnSurf
from OCC.Core.GeomAbs import GeomAbs_C0, GeomAbs_C1, GeomAbs_C2
from OCC.Core.GeomAbs import GeomAbs_G1, GeomAbs_G2
from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeShell
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeSolid


def spl_face(px, py, pz, axs=gp_Ax3()):
    nx, ny = px.shape
    pnt_2d = TColgp_Array2OfPnt(1, nx, 1, ny)
    for row in range(pnt_2d.LowerRow(), pnt_2d.UpperRow() + 1):
        for col in range(pnt_2d.LowerCol(), pnt_2d.UpperCol() + 1):
            i, j = row - 1, col - 1
            pnt = gp_Pnt(px[i, j], py[i, j], pz[i, j])
            pnt_2d.SetValue(row, col, pnt)
            #print (i, j, px[i, j], py[i, j], pz[i, j])

    api = GeomAPI_PointsToBSplineSurface(pnt_2d, 0, 0, GeomAbs_C0, 0.001)
    api.Interpolate(pnt_2d)
    #surface = BRepBuilderAPI_MakeFace(curve, 1e-6)
    # return surface.Face()
    surf = api.Surface()
    surf.Transform(set_trf(gp_Ax3(), axs))
    face = BRepBuilderAPI_MakeFace(surf, 1e-6).Face()
    shel = BRepBuilderAPI_MakeShell(surf).Shell()
    return surf, face, shel


if __name__ == '__main__':
    argvs = sys.argv
    parser = OptionParser()
    parser.add_option("--dir", dest="dir", default="./")
    parser.add_option("--pxyz", dest="pxyz",
                      default=[0.0, 0.0, 0.0], type="float", nargs=3)
    opt, argc = parser.parse_args(argvs)
    print(opt, argc)

    px = np.linspace(-1, 1, 100) * 150
    py = np.linspace(-1, 1, 200) * 150
    axis1 = gp_Ax3(gp_Pnt(0, 0, -5), gp_Dir(0.0, 0.05, 0.9))
    mesh1 = np.meshgrid(px, py)
    data1 = mesh1[0]**2 / 1000 + mesh1[1]**2 / 2000
    surf1, face1, shel1 = spl_face(*mesh1, data1, axis1)

    px = np.linspace(-1, 1, 100) * 150
    py = np.linspace(-1, 1, 200) * 150
    axis2 = gp_Ax3(gp_Pnt(0, 0, 5), gp_Dir(0.1, 0.05, 0.9))
    mesh2 = np.meshgrid(px, py)
    data2 = mesh2[0]**2 / 2000 - mesh2[1]**2 / 1000
    surf2, face2, shel2 = spl_face(*mesh2, data2, axis2)

    solid = BRepBuilderAPI_MakeSolid(shel1, shel2).Solid()
    surf_face1 = BRep_Tool.Surface(face1)
    surf_face2 = BRep_Tool.Surface(face2)
    intss = GeomAPI_IntSS(surf_face1, surf_face2, 0.1)
    print(intss.NbLines())
    intss_curv1 = intss.Line(1)

    obj = dispocc(touch=True)
    obj.show_axs_pln(axis1, scale=50, name="Axis1")
    obj.show_axs_pln(axis2, scale=50, name="Axis2")
    obj.display.DisplayShape(surf1, color="RED", transparency=0.9)
    obj.display.DisplayShape(surf2, color="BLUE1", transparency=0.9)
    obj.display.DisplayShape(solid, color="GREEN", transparency=0.9)
    obj.display.DisplayShape(intss_curv1)
    obj.export_stp(solid)
    #obj.export_stl(solid, linear_deflection=0.1, angular_deflection=0.1)
    obj.ShowOCC()
