import numpy as np
import matplotlib.pyplot as plt
import sys
import os

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Lin, gp_Sphere
from OCC.Core.gp import gp_Mat, gp_XYZ
from OCC.Core.gp import gp_Trsf, gp_GTrsf
from OCC.Core.TColgp import TColgp_Array1OfPnt, TColgp_Array2OfPnt
from OCC.Core.TColGeom import TColGeom_HArray2OfSurface
from OCC.Core.TopTools import TopTools_ListOfShape
from OCC.Core.Geom import Geom_BezierSurface, Geom_BSplineSurface, Geom_RectangularTrimmedSurface
from OCC.Core.GeomAPI import GeomAPI_PointsToBSplineSurface, GeomAPI_ProjectPointOnSurf
from OCC.Core.GeomAbs import GeomAbs_C0, GeomAbs_C1, GeomAbs_C2
from OCC.Core.GeomAbs import GeomAbs_G1, GeomAbs_G2
from OCC.Core.GeomProjLib import geomprojlib_ProjectOnPlane
from OCC.Core.ShapeExtend import ShapeExtend_CompositeSurface
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.Core.BRepProj import BRepProj_Projection
from OCC.Core.BRepOffset import BRepOffset_Analyse
from OCC.Core.BRep import BRep_Tool
from OCC.Core.GeomPlate import GeomPlate_MakeApprox, GeomPlate_BuildPlateSurface
from OCC.Core.GeomToStep import GeomToStep_MakeBoundedSurface
from OCC.Core.BRepMesh import BRepMesh_IncrementalMesh
from OCC.Core.BRepFill import BRepFill_CurveConstraint
from OCC.Core.BRepFeat import BRepFeat_SplitShape
from OCC.Core.BRepAlgo import BRepAlgo_Cut, BRepAlgo_Section, BRepAlgo_Common
from OCC.Core.BOPAlgo import BOPAlgo_Tools, BOPAlgo_Builder
from OCC.Core.BOPTools import BOPTools_AlgoTools
from OCC.Core.STEPControl import STEPControl_Reader, STEPControl_Writer, STEPControl_AsIs
from OCC.Core.Interface import Interface_Static_SetCVal
from OCC.Core.IFSelect import IFSelect_RetDone, IFSelect_ItemsByEntity
from OCC.Extend.DataExchange import write_step_file, read_step_file, write_stl_file
from OCCUtils.Construct import make_edge, make_polygon

from src.base import dispocc, set_loc, set_trf


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
    return surf, face


def bez_face(px, py, pz, axs=gp_Ax3()):
    nx, ny = px.shape
    pnt_2d = TColgp_Array2OfPnt(1, nx, 1, ny)
    for row in range(pnt_2d.LowerRow(), pnt_2d.UpperRow() + 1):
        for col in range(pnt_2d.LowerCol(), pnt_2d.UpperCol() + 1):
            i, j = row - 1, col - 1
            pnt = gp_Pnt(px[i, j], py[i, j], pz[i, j])
            pnt_2d.SetValue(row, col, pnt)
            #print (i, j, px[i, j], py[i, j], pz[i, j])

    surf = Geom_BezierSurface(pnt_2d)
    surf.Transform(set_trf(gp_Ax3(), axs))
    # api.Interpolate(pnt_2d)
    #surface = BRepBuilderAPI_MakeFace(curve, 1e-6)
    # return surface.Face()
    face = BRepBuilderAPI_MakeFace(surf, 1e-6).Face()
    return surf, face


if __name__ == '__main__':
    obj = dispocc(touch=True)
    obj.show_axs_pln(obj.base_axs, scale=1)

    # OCCT dev
    # https://dev.opencascade.org/doc/overview/html/occt_user_guides__modeling_data.html
    #
    # BRepBuilderAPI_MakeFace
    # https://www.opencascade.com/doc/occt-7.5.0/refman/html/class_b_rep_builder_a_p_i___make_face.html
    #
    # Bezier
    # https://old.opencascade.com/doc/occt-7.5.0/refman/html/class_geom___bezier_surface.html#details
    # the degree for a Geom_BezierSurface is limited to a value of (25) which is defined and controlled by the system

    px = np.linspace(-1, 1, 100) * 300
    py = np.linspace(-1, 1, 100) * 300
    mesh = np.meshgrid(px, py)
    data1 = mesh[0]**2 / 2000 + mesh[1]**2 / 1000

    px = np.linspace(-1, 1, 25) * 150
    py = np.linspace(-1, 1, 25) * 150
    mesh = np.meshgrid(px, py)
    data2 = mesh[0]**2 / 1000 + mesh[1]**2 / 2000

    pts = []
    pts.append(gp_Pnt(-50, -60, 0))
    pts.append(gp_Pnt(-60, +70, 0))
    pts.append(gp_Pnt(+70, +80, 0))
    pts.append(gp_Pnt(+80, -90, 0))
    poly = make_polygon(pts, closed=True)

    axis1 = gp_Ax3(gp_Pnt(-150, +50.0, 0.0), gp_Dir(0, 0, 1))
    surf1, face1 = spl_face(*mesh, data1, axs=axis1)
    poly1 = poly.Located(set_loc(gp_Ax3(), axis1))
    poly1_proj = obj.proj_rim_pln(poly1, face1, axis1)
    print(face1)

    #split = BRepFeat_SplitShape(face1)
    #split.Add(poly1_proj, face1)
    # split.Build()
    #surf1_trim1 = split.Left()
    #surf1_trim2 = split.DirectLeft()
    # print(split.Check())
    #print(surf1_trim1.Size(), surf1_trim2.Size())

    face1_holl = BRepBuilderAPI_MakeFace(face1, poly1_proj).Face()
    print(face1_holl)
    api = BRepAlgo_Cut(face1, face1_holl)
    api.Build()
    face1_trim = api.Shape()
    print(face1_trim)

    axis2 = gp_Ax3(gp_Pnt(+100, -50.0, 0.0), gp_Dir(0, 0, 1))
    surf2, face2 = bez_face(*mesh, data2, axs=axis2)
    u0, u1, v0, v1 = surf2.Bounds()
    surf2_uiso = surf2.UIso(0.5)
    surf2_viso = surf2.VIso(0.5)
    poly2 = poly.Located(set_loc(gp_Ax3(), axis2))
    poly2_proj = obj.proj_rim_pln(poly2, face2, axis1)
    print(face2, u0, u0, v0, v1)
    # write_stl_file(face2_stl, obj.tempname + "_face2_001.stp",
    #               linear_deflection=0.1, angular_deflection=0.1)
    # write_stl_file(face2_stl, obj.tempname + "_face2_002.stp",
    #               linear_deflection=0.001, angular_deflection=0.001)

    surf2_trim = Geom_RectangularTrimmedSurface(
        surf2, 0.1, 0.8, 0.5, 0.7, True, True)

    face3 = read_step_file("surf1.step")
    print(face3)
    #surf3 = BRep_Tool.Surface(face3)
    # print(surf3)

    mesh = BRepMesh_IncrementalMesh(face2, 0.01, True, 0.01, True)
    mesh.Perform()
    face3_mesh = mesh.Shape()

    obj.display.DisplayShape(surf1, color="BLUE1", transparency=0.5)
    obj.display.DisplayShape(poly1)
    obj.display.DisplayShape(poly1_proj)
    obj.display.DisplayShape(face1_holl, color="GREEN", transparency=0.5)
    #obj.display.DisplayShape(face1_trim, color="BLACK", transparency=0.5)
    obj.display.DisplayShape(surf2, color="RED", transparency=0.9)
    obj.display.DisplayShape(poly2)
    obj.display.DisplayShape(poly2_proj)
    obj.display.DisplayShape(face3, color="YELLOW", transparency=0.9)
    obj.display.DisplayShape(surf2_trim, color="RED", transparency=0.5)
    obj.display.DisplayShape(surf2_uiso)
    obj.display.DisplayShape(surf2_viso)
    obj.show_occ()
