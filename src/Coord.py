import numpy as np
import matplotlib.pyplot as plt
import json
import glob
import sys
import time
import os
import scipy.constants as cnt
from scipy.integrate import simps
import argparse
from linecache import getline, clearcache

from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Pnt, gp_Dir, gp_Vec
from OCC.Core.gp import gp_Pln, gp_Circ, gp_Lin, gp_Elips
from OCC.Core.gp import gp_Trsf, gp_Quaternion
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Core.Geom import Geom_Curve, Geom_Plane, Geom_Line
from OCC.Core.GeomAPI import GeomAPI_PointsToBSplineSurface, GeomAPI_ProjectPointOnSurf, GeomAPI_ProjectPointOnCurve
from OCC.Core.GeomAPI import GeomAPI_IntCS, GeomAPI_IntSS
from OCC.Core.GeomLProp import GeomLProp_SurfaceTool, GeomLProp_CurveTool
from OCC.Core.GeomAbs import GeomAbs_C3, GeomAbs_G2
from OCC.Core.TColgp import TColgp_Array2OfPnt
from OCC.Core.TopOpeBRepTool import TopOpeBRepTool_TOOL
from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepTools import breptools
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
from OCC.Core.BRepIntCurveSurface import BRepIntCurveSurface_Inter
from OCC.Core.BRepLProp import BRepLProp_SLProps
from OCC.Core.BRepProj import BRepProj_Projection
from OCCUtils.Topology import Topo
from OCCUtils.Construct import vec_to_dir, dir_to_vec
from OCCUtils.Construct import point_to_vector, vector_to_point
from OCCUtils.Construct import make_polygon, make_plane, make_edge, make_line

sys.path.append(os.path.join('../'))
from src.base_occ import dispocc
from src.base_occ import float_to_string, pnt_to_xyz, rotate_xyz, set_loc, set_trf
from src.base_occ import get_axs, get_deg, grasp_sfc, occ_to_grasp_cor, occ_to_grasp_cor_ref, surf_spl_pcd
from src.geometry import curvature


class CoordSys (dispocc):

    def __init__(self):
        # dispocc.__init__(self, touch=False)
        self.axs = gp_Ax3()
        self.org = gp_Ax3()

    def SetAxs(self, pxyz=[0, 0, 0], rxyz=[0, 0, 0]):
        self.axs = gp_Ax3(gp_Pnt(*pxyz), gp_Dir(0, 0, 1))
        rotate_xyz(self.axs, rxyz[0], "x")
        rotate_xyz(self.axs, rxyz[1], "y")
        rotate_xyz(self.axs, rxyz[2], "z")

    def TrfAxs(self, pxyz=[0, 0, 0], rxyz=[0, 0, 0]):
        self.axs.Translate(gp_Pnt(), gp_Pnt(*pxyz))
        rotate_xyz(self.axs, rxyz[0], "x")
        rotate_xyz(self.axs, rxyz[1], "y")
        rotate_xyz(self.axs, rxyz[2], "z")

    def TrfAxs_Ref(self, pxyz=[0, 0, 0], rxyz=[0, 0, 0], ax=gp_Ax3()):
        trf = gp_Trsf()
        axs = gp_Ax3()
        axs.SetLocation(gp_Pnt(*pxyz))
        rotate_xyz(axs, rxyz[0], "x")
        rotate_xyz(axs, rxyz[1], "y")
        rotate_xyz(axs, rxyz[2], "z")
        trf.SetTransformation(axs, gp_Ax3())
        self.axs = ax.Transformed(trf)

    def export_axs(self, dirname="./", name="surf"):
        filename = dirname + name + ".cor"
        occ_to_grasp_cor(self.axs, name, filename)

    def import_axs(self, dirname="./", name="surf"):
        filename = dirname + name + ".cor"
        return get_axs(filename)

    def get_diff(self):
        vx = dir_to_vec(self.org.XDirection())
        vy = dir_to_vec(self.org.YDirection())
        vz = dir_to_vec(self.org.Direction())
        pln_x = Geom_Plane(self.org.Location(), self.org.YDirection())
        pln_y = Geom_Plane(self.org.Location(), self.org.XDirection())

        vec = dir_to_vec(self.axs.Direction())
        vec_p = gp_Pnt((gp_Vec(self.org.Location().XYZ()) + vec).XYZ())
        pnt_x = GeomAPI_ProjectPointOnSurf(vec_p, pln_x).Point(1)
        pnt_y = GeomAPI_ProjectPointOnSurf(vec_p, pln_y).Point(1)
        vec_x = gp_Vec(self.org.Location(), pnt_x)
        vec_y = gp_Vec(self.org.Location(), pnt_y)
        deg_x = vec_x.AngleWithRef(vz, vy)
        deg_y = vec_y.AngleWithRef(vz, vx)

        p0 = self.org.Location()
        p1 = self.axs.Location()
        print(self.org.Location().Coord())
        print(self.axs.Location().Coord())

        txt = ""
        txt += "deg_x = {:+.2f}\n".format(np.rad2deg(deg_x))
        txt += "deg_y = {:+.2f}\n".format(np.rad2deg(deg_y))
        txt += "dist = {:.2f}\n".format(p0.Distance(p1))
        print(txt)


class OCCSurfObj(object):

    def __init__(self, name="surf"):
        self.name = name
        self.rot = gp_Ax3()
        self.axs = gp_Ax3()
        self.rim = make_edge(gp_Circ(self.axs.Ax2(), 100))
        self.pln = dispocc.make_plane_axs(self.axs)
        self.surf = make_plane(self.axs.Location(),
                               dir_to_vec(self.axs.Direction()),
                               -500, 500, -500, 500)

        # Beam
        num = 4
        self.pt = np.linspace(0, 2 * np.pi, num + 1)[:-1]
        self.wxy = [10, 10]
        self.beam = gp_Ax3(self.axs.Ax2())
        self.beam_wxy = self.make_beam_wxy()
        self.beam_pts = self.make_beam_pts()

    def MovTrfSurf(self, trf=gp_Trsf(), rot=False):
        if rot == True:
            self.axs.Transform(trf)
        self.axs.Transform(trf)
        self.rim.Move(TopLoc_Location(trf))
        self.pln.Move(TopLoc_Location(trf))
        self.surf.Move(TopLoc_Location(trf))

    def get_trsf(self):
        self.trf = gp_Trsf()
        self.trf.SetTransformation(self.axs, gp_Ax3())
        return self.trf

    def get_beam_local(self):
        axs = self.beam.Transformed(set_trf(self.axs, gp_Ax3()))
        print(self.name, "-", "beam")
        get_deg(self.axs, dir_to_vec(self.beam.Direction()))
        print(axs.Location())
        print(gp_Vec(axs.Direction()))
        return axs

    def get_vxyz(self):
        trf = gp_Trsf()
        trf.SetTransformation(gp_Ax3(), self.rot)
        # trf.Invert()
        ax1 = self.axs.Transformed(trf)
        d_z = ax1.Direction()
        print(dir_to_vec(d_z),
              get_deg(self.axs, dir_to_vec(self.rot.Direction())))
        return [d_z.X(), d_z.Y(), d_z.Z()]

    def RotateSurf(self, deg=0.0, axs="z", rot=False):
        if axs == "x":
            ax1 = gp_Ax1(self.rot.Location(), self.rot.XDirection())
        elif axs == "y":
            ax1 = gp_Ax1(self.rot.Location(), self.rot.YDirection())
        elif axs == "z":
            ax1 = gp_Ax1(self.rot.Location(), self.rot.Direction())
        else:
            ax1 = gp_Ax1(self.rot.Location(), self.rot.Direction())
        trf = gp_Trsf()
        trf.SetRotation(ax1, np.deg2rad(deg))
        self.MovTrfSurf(trf, rot)

    def MovXYZSurf(self, dst=0.0, axs="z", rot=False):
        if axs == "x":
            vec = dir_to_vec(self.rot.XDirection())
        elif axs == "y":
            vec = dir_to_vec(self.rot.YDirection())
        elif axs == "z":
            vec = dir_to_vec(self.rot.Direction())
        else:
            vec = dir_to_vec(self.rot.Direction())
        trf = gp_Trsf()
        trf.SetTranslation(vec.Scaled(dst))
        self.MovTrfSurf(trf, rot)

    def MovVecSurf(self, vec=gp_Vec(0, 0, 1), rot=False):
        trf = gp_Trsf()
        trf.SetTranslation(vec)
        self.MovTrfSurf(trf, rot)

    def RotateAxis(self, deg=0.0, axs="z"):
        if axs == "x":
            ax1 = gp_Ax1(self.rot.Location(), self.rot.XDirection())
        elif axs == "y":
            ax1 = gp_Ax1(self.rot.Location(), self.rot.YDirection())
        elif axs == "z":
            ax1 = gp_Ax1(self.rot.Location(), self.rot.Direction())
        else:
            ax1 = gp_Ax1(self.rot.Location(), self.rot.Direction())
        rot = self.rot.Rotated(ax1, np.deg2rad(deg))
        trf = gp_Trsf()
        trf.SetDisplacement(self.rot, rot)
        self.axs.Transform(trf)

    def RotateAxis_Ax(self, ax=gp_Ax3(), deg=0.0, axs="z"):
        if axs == "x":
            ax1 = gp_Ax1(ax.Location(), ax.XDirection())
        elif axs == "y":
            ax1 = gp_Ax1(ax.Location(), ax.YDirection())
        elif axs == "z":
            ax1 = gp_Ax1(ax.Location(), ax.Direction())
        else:
            ax1 = gp_Ax1(ax.Location(), ax.Direction())
        rot = ax.Rotated(ax1, np.deg2rad(deg))
        trf = gp_Trsf()
        trf.SetDisplacement(ax, rot)
        ax.Transform(trf)

    def make_mesh(self, lxy=[100, 100], sxy=[0, 0], nxy=[200, 200]):
        px = np.linspace(-1, 1, nxy[0]) * lxy[0] / 2 + sxy[0]
        py = np.linspace(-1, 1, nxy[1]) * lxy[1] / 2 + sxy[1]
        self.mesh = np.meshgrid(px, py)

    def make_surf_curvature(self, rxy=[1000, 1000], sxy=[0, 0]):
        curv_x = curvature(self.mesh[0], rxy[0], sxy[0])
        curv_y = curvature(self.mesh[1], rxy[1], sxy[1])
        self.surf_data = curv_x + curv_y
        self.surf, self.surf_pts = surf_spl_pcd(*self.mesh, self.surf_data)
        self.surf.Move(TopLoc_Location(self.get_trsf()))

    def make_beam_wxy(self):
        rim = dispocc.make_EllipWire(None, self.wxy, 0, self.beam, None)
        beam_wxy = dispocc.proj_rim_pln(None, rim, self.surf, self.axs)
        return beam_wxy

    def make_beam_pts(self):
        pts = []
        for i, t in enumerate(self.pt):
            x, y = self.wxy[0] * np.cos(t), self.wxy[1] * np.sin(t)
            pnt = gp_Pnt(x, y, 0)
            pnt.Transform(set_trf(gp_Ax3(), self.beam))
            pnt = dispocc.proj_pnt_pln(None, pnt, self.surf, self.beam)
            pts.append(pnt)
        return pts

    def get_surf_mesh(self):
        nx = self.surf_pts.UpperRow()
        ny = self.surf_pts.UpperCol()
        x, y, z = np.zeros((nx, ny)), np.zeros((nx, ny)), np.zeros((nx, ny))
        for row in range(self.surf_pts.LowerRow(), self.surf_pts.UpperRow() + 1):
            for col in range(self.surf_pts.LowerCol(), self.surf_pts.UpperCol() + 1):
                i, j = row - 1, col - 1
                pnt = self.surf_pts.Value(row, col)
                x[i, j], y[i, j], z[i, j] = pnt.Coord()
        return [x, y], z

    def get_surf_uvpnt(self, uv=[0, 0]):
        surf = BRep_Tool.Surface(self.surf)
        pnt = gp_Pnt()
        surf.D0(uv[0], uv[1], pnt)
        return pnt

    def reflect_beam(self, beam0=gp_Ax3(), tr=0):
        v0 = dir_to_vec(beam0.Direction())
        v1 = dir_to_vec(beam0.XDirection())
        p0 = beam0.Location()
        lin = gp_Lin(beam0.Axis())
        api = BRepIntCurveSurface_Inter()
        api.Init(self.surf, lin, 1.0E-9)
        dst = np.inf
        num = 0
        sxy = p0
        uvw = [0, 0, 0]
        fce = None
        while api.More():
            p1 = api.Pnt()
            dst1 = p0.Distance(p1)
            if dst1 < dst and api.W() > 1.0E-6:
                dst = dst1
                uvw = [api.U(), api.V(), api.W()]
                sxy = api.Pnt()
                fce = api.Face()
                api.Next()
            else:
                api.Next()
        if fce == None:
            fce = make_plane(self.axs.Location(), dir_to_vec(self.axs.Direction()),
                             -10000, 10000, -10000, 10000)
            api.Init(fce, lin, 1.0E-9)
            uvw = [api.U(), api.V(), api.W()]

        print(self.name, *uvw)
        u, v, w = uvw
        surf = BRepAdaptor_Surface(fce)
        prop = BRepLProp_SLProps(surf, u, v, 2, 1.0E-9)
        p1, vx, vy = prop.Value(), prop.D1U(), prop.D1V()
        vz = vx.Crossed(vy)
        if vz.Dot(v0) > 0:
            vz.Reverse()
        vx.Normalize()
        vy.Normalize()
        vz.Normalize()

        self.beam = gp_Ax3(
            p1,
            vec_to_dir(v0.Reversed()),
            vec_to_dir(v1.Reversed())
        )
        self.norm = gp_Ax3(
            p1,
            vec_to_dir(vz),
            vec_to_dir(vx)
        )
        if tr == 0:
            # Reflect
            self.beam.Mirror(self.norm.Ax2())
            if self.beam.Direction().Dot(self.norm.Direction()) < 0:
                self.beam.ZReverse()
        elif tr == 1:
            # Transporse
            self.beam.ZReverse()
            self.beam.XReverse()
            # if self.beam.Direction().Dot(self.norm.Direction()) < 0:
            #    self.beam.ZReverse()
        # print(self.beam.Direction().Dot(self.norm.Direction()))


if __name__ == "__main__":
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="../temp/")
    parser.add_argument("--surf", dest="surf", default="surf")
    parser.add_argument("--pxyz", dest="pxyz",
                        default=(0, 0, 0), type=float, nargs=3)
    parser.add_argument("--rxyz", dest="rxyz",
                        default=(0, 0, 0), type=float, nargs=3)
    opt = parser.parse_args()
    print(opt, argvs)

    obj = CoordSys()
    obj.SetAxs(opt.pxyz, opt.rxyz)

    obj.axs = obj.import_axs(obj.tmpdir, name="surf")
    obj.TrfAxs([10, 20, 30], [1, 2, 3])

    obj.axs = obj.import_axs(obj.tmpdir, name="surf")
    obj.TrfAxs([10, 10, 0], [45.0, 0.0, 0.0])

    obj.axs = gp_Ax3()
    obj.org = obj.import_axs(obj.tmpdir, name="surf2")
    obj.TrfAxs_Ref(pxyz=[10, 20, 0], rxyz=[1.0, 0.0, 0.0])
    obj.get_diff()

    obj.axs = gp_Ax3()
    obj.org = obj.import_axs(obj.tmpdir, name="surf2")
    obj.TrfAxs(pxyz=[10, 20, 0], rxyz=[1.0, 0.0, 0.0])
    obj.get_diff()
