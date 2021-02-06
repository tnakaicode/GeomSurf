import numpy as np
import matplotlib.pyplot as plt
import sys

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Lin, gp_Sphere
from OCC.Core.gp import gp_Mat, gp_XYZ
from OCC.Core.gp import gp_Trsf, gp_GTrsf
from OCC.Core.TopoDS import TopoDS_Shape
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Core.TColgp import TColgp_Array1OfPnt, TColgp_Array2OfPnt
from OCC.Core.TColgp import TColgp_HArray1OfPnt, TColgp_HArray2OfPnt
from OCC.Core.BRepFill import BRepFill_Filling
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeSphere
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_GTransform
from OCC.Core.Geom import Geom_SphericalSurface, Geom_OsculatingSurface
from OCC.Core.Geom import Geom_Line
from OCC.Core.GeomAPI import geomapi
from OCC.Core.GeomAPI import GeomAPI_IntCS
from OCC.Core.GeomAPI import GeomAPI_ExtremaCurveSurface
from OCC.Core.GeomAPI import GeomAPI_PointsToBSplineSurface
from OCC.Core.GeomAPI import GeomAPI_PointsToBSpline
from OCC.Core.GeomAPI import GeomAPI_Interpolate
from OCC.Core.GeomAbs import GeomAbs_C0, GeomAbs_C1, GeomAbs_C2
from OCC.Core.GeomAbs import GeomAbs_G1, GeomAbs_G2
from OCCUtils.Topology import TopExp_Explorer, Topo
from OCCUtils.Topology import TopAbs_FACE
from OCCUtils.Construct import make_n_sided, make_n_sections
from OCCUtils.Construct import make_edge

from OCC.Core.GeomAPI import GeomAPI_IntCS
from OCC.Core.ShapeAnalysis import ShapeAnalysis_ShapeTolerance
from OCC.Core.BRepTools import breptools_UVBounds
from OCC.Core.BRep import BRep_Tool_Surface

from src.base import plotocc


def set_trf(ax1=gp_Ax3(), ax2=gp_Ax3()):
    trf = gp_Trsf()
    trf.SetTransformation(ax2, ax1)
    return trf


def set_loc(ax1=gp_Ax3(), ax2=gp_Ax3()):
    trf = set_trf(ax1, ax2)
    loc = TopLoc_Location(trf)
    return loc


def trsf_scale(axs=gp_Ax3(), scale=1):
    trf = gp_Trsf()
    trf.SetDisplacement(gp_Ax3(), axs)
    return trf


def gen_ellipsoid(axs=gp_Ax3(), rxyz=[10, 20, 30]):
    sphere = BRepPrimAPI_MakeSphere(gp_Ax2(), 1).Solid()
    loc = set_loc(gp_Ax3(), axs)
    mat = gp_Mat(
        rxyz[0], 0, 0,
        0, rxyz[1], 0,
        0, 0, rxyz[2]
    )
    gtrf = gp_GTrsf(mat, gp_XYZ(0, 0, 0))
    ellips = BRepBuilderAPI_GTransform(sphere, gtrf).Shape()
    ellips.Location(loc)
    return ellips


def gen_ellipsoid_geom(axs=gp_Ax3(), rxyz=[10, 20, 30]):
    sphere = gp_Sphere(axs, 1)
    api = Geom_SphericalSurface(sphere)
    api.Continuity()
    api.Transform()
    print(api.Area())

    loc = set_loc(gp_Ax3(), axs)
    mat = gp_Mat(
        rxyz[0], 0, 0,
        0, rxyz[1], 0,
        0, 0, rxyz[2]
    )
    gtrf = gp_GTrsf(mat, gp_XYZ(0, 0, 0))
    ellips = BRepBuilderAPI_GTransform(sphere, gtrf).Shape()
    ellips.Location(loc)
    return ellips


if __name__ == '__main__':
    obj = plotocc(touch=True)
    # obj.show_axs_pln(scale=20)
    # obj.show_ball()

    axs = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(1, 1, 1))
    elp = gen_ellipsoid(axs, [10, 20, 30])
    obj.display.DisplayShape(elp, transparency=0.7, color="BLUE")
    obj.show_axs_pln(axs, scale=20)

    axs = gp_Ax3(gp_Pnt(30, 0, 0), gp_Dir(1, 1, 0))
    elp = gen_ellipsoid(axs, [10, 20, 30])
    obj.display.DisplayShape(elp, transparency=0.7, color="BLUE")
    obj.show_axs_pln(axs, scale=20)

    #elp = gen_ellipsoid_geom(axs, [10, 20, 30])

    top_api = Topo(elp)
    print(top_api.number_of_faces())
    for face in top_api.faces():
        elp_face = face

    print(elp_face)
    elp_surf = BRep_Tool_Surface(elp_face)
    print(elp_surf)

    lin = Geom_Line(gp_Ax1(axs.Location(), gp_Dir(0, 1, 1)))
    api = GeomAPI_IntCS(lin, elp_surf)
    obj.display.DisplayShape(lin)
    print(api.Point(1))
    print(api.Point(2))

    obj.show()
