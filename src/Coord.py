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
from OCC.Core.gp import gp_Pln, gp_Circ, gp_Lin
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
from base_occ import dispocc
from base_occ import float_to_string, pnt_to_xyz, rotate_xyz, set_loc, set_trf
from base_occ import get_axs, grasp_sfc, occ_to_grasp_cor, occ_to_grasp_cor_ref, surf_spl_pcd


def curvature(px, r, s):
    """( x + sx )**2 / 2*rx + ( y + sy )**2 / 2*ry"""
    if (r == 0):
        py = np.zeros_like(px + s)
    else:
        py = (px + s)**2 / (2 * r)
    return py


def get_deg(axs, vec):
    vx = dir_to_vec(axs.XDirection())
    vy = dir_to_vec(axs.YDirection())
    vz = dir_to_vec(axs.Direction())
    pln_x = Geom_Plane(axs.Location(), axs.YDirection())
    pln_y = Geom_Plane(axs.Location(), axs.XDirection())
    vec_p = gp_Pnt((gp_Vec(axs.Location().XYZ()) + vec).XYZ())
    pnt_x = GeomAPI_ProjectPointOnSurf(vec_p, pln_x).Point(1)
    pnt_y = GeomAPI_ProjectPointOnSurf(vec_p, pln_y).Point(1)
    vec_x = gp_Vec(axs.Location(), pnt_x)
    vec_y = gp_Vec(axs.Location(), pnt_y)
    deg_x = vec_x.AngleWithRef(vz, vy)
    deg_y = vec_y.AngleWithRef(vz, vx)
    print(np.rad2deg(deg_x), np.rad2deg(deg_y))
    return deg_x, deg_y


class CoordSys (dispocc):

    def __init__(self):
        #dispocc.__init__(self, touch=False)
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
        self.axs = gp_Ax3()
        self.rot = self.axs
        self.rim = make_edge(gp_Circ(self.axs.Ax2(), 100))
        self.pln = dispocc.make_plane_axs(self.axs)
        self.surf = make_plane(
            self.axs.Location(), dir_to_vec(self.axs.Direction()),
            -500, 500, -500, 500)
        self.face = self.BuilFace()
        self.name = name

    def get_trsf(self):
        self.trf = gp_Trsf()
        self.trf.SetTransformation(self.axs, gp_Ax3())
        return self.trf

    def get_vxyz(self):
        trf = gp_Trsf()
        trf.SetTransformation(gp_Ax3(), self.rot)
        # trf.Invert()
        ax1 = self.axs.Transformed(trf)
        d_z = ax1.Direction()
        print(dir_to_vec(d_z), get_deg(
            self.axs, dir_to_vec(self.rot.Direction())))
        return [d_z.X(), d_z.Y(), d_z.Z()]

    def BuilFace(self):
        proj = BRepProj_Projection(self.rim, self.surf, self.axs.Direction())
        bild = BRepBuilderAPI_MakeFace(self.surf, proj.Current())
        return bild.Face()

    def UpdateFace(self):
        self.face = self.BuilFace()

    def RotateFace(self, deg=0.0, axs="z"):
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
        self.face.Move(TopLoc_Location(trf))

    def RotateSurf(self, deg=0.0, axs="z"):
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
        self.rot.Transform(trf)
        self.axs.Transform(trf)
        self.surf.Move(TopLoc_Location(trf))

    def RotateSurf2(self, deg=0.0, axs="z"):
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
        self.axs.Transform(trf)
        self.surf.Move(TopLoc_Location(trf))

    def SetSurf_XY(self, deg=0.0, axs="x"):
        dx, dy = get_deg(self.axs, dir_to_vec(self.rot.Direction()))
        if axs == "x":
            ax1 = gp_Ax1(self.rot.Location(), self.rot.XDirection())
            dg0 = np.rad2deg(dy)
        elif axs == "y":
            ax1 = gp_Ax1(self.rot.Location(), self.rot.YDirection())
            dg0 = np.rad2deg(dx)
        else:
            ax1 = gp_Ax1(self.rot.Location(), self.rot.Direction())
            dg0 = np.rad2deg(dx)
        trf = gp_Trsf()
        trf.SetRotation(ax1, np.deg2rad(-dg0 + deg))
        self.axs.Transform(trf)
        self.surf.Move(TopLoc_Location(trf))
        get_deg(self.rot, dir_to_vec(self.axs.Direction()))

    def MovXYZSurf(self, dst=0.0, axs="z"):
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
        # self.rot.Transform(trf)
        self.axs.Transform(trf)
        self.surf.Move(TopLoc_Location(trf))

    def MovVecSurf(self, vec=gp_Vec(0, 0, 1)):
        trf = gp_Trsf()
        trf.SetTranslation(vec)
        self.rot.Transform(trf)
        self.axs.Transform(trf)
        self.surf.Move(TopLoc_Location(trf))

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

    def MoveRel(self, trf=gp_Trsf()):
        self.axs.Transform(trf)
        self.face.Move(TopLoc_Location(trf))

    def SurfCurvature(self, nxy=[200, 200], lxy=[450, 450], rxy=[700, 0], sxy=[0, 0]):
        px = np.linspace(-1, 1, int(nxy[0])) * lxy[0] / 2
        py = np.linspace(-1, 1, int(nxy[1])) * lxy[1] / 2
        mesh = np.meshgrid(px, py)
        surf_x = curvature(mesh[0], r=rxy[0], s=sxy[0])
        surf_y = curvature(mesh[1], r=rxy[1], s=sxy[1])
        data = surf_x + surf_y
        self.surf, self.surf_pts = surf_spl_pcd(*mesh, data)

    def SurfCurvature_Loc(self, nxy=[200, 200], lxy=[450, 450], rxy=[700, 0], sxy=[0, 0]):
        px = np.linspace(-1, 1, int(nxy[0])) * lxy[0] / 2
        py = np.linspace(-1, 1, int(nxy[1])) * lxy[1] / 2
        mesh = np.meshgrid(px, py)
        surf_x = curvature(mesh[0], r=rxy[0], s=sxy[0])
        surf_y = curvature(mesh[1], r=rxy[1], s=sxy[1])
        data = surf_x + surf_y
        self.surf, self.surf_pts = surf_spl_pcd(*mesh, data)
        trf = gp_Trsf()
        trf.SetTransformation(self.axs, gp_Ax3())
        self.surf.Location(TopLoc_Location(trf))

    def get_surf_uvpnt(self, uv=[0, 0]):
        surf = BRep_Tool.Surface(self.surf)
        pnt = gp_Pnt()
        surf.D0(uv[0], uv[1], pnt)
        return pnt

    def load_rim(self, rimfile="../ticra/input/surf/mou.rim"):
        data = np.loadtxt(rimfile, skiprows=2)
        pts = []
        for xy in data:
            pts.append(gp_Pnt(*xy, 0))
        self.rim = make_polygon(pts, closed=True)

    def load_mat(self, sfcfile="../ticra/input/surf/pln_mat.sfc"):
        xs, ys, xe, ye = [float(v) for v in getline(sfcfile, 2).split()]
        nx, ny = [int(v) for v in getline(sfcfile, 3).split()]
        px = np.linspace(xs, xe, nx)
        py = np.linspace(ys, ye, ny)
        mesh = np.meshgrid(px, py)
        data = np.loadtxt(sfcfile, skiprows=3).T
        self.surf, self.surf_pts = surf_spl_pcd(*mesh, data)

    def export_rim_2d(self, rimfile="m2.rim", name="m2-rim"):
        rim_2d = dispocc.proj_rim_pln(self, self.rim, self.pln, self.axs)

        fp = open(rimfile, "w")
        fp.write(' {:s}\n'.format(name))
        fp.write('{:12d}{:12d}{:12d}\n'.format(1, 1, 1))
        rim_tmp = gp_Pnt()
        for i, e in enumerate(Topo(rim_2d).edges()):
            e_curve, u0, u1 = BRep_Tool.Curve(e)
            print(i, e, u0, u1)
            if i != 0 and rim_tmp == e_curve.Value(u0):
                u_range = np.linspace(u0, u1, 50)
                rim_tmp = e_curve.Value(u1)
                p = e_curve.Value(u0)
                p.Transform(set_trf(self.axs, gp_Ax3()))
                data = [p.X(), p.Y()]
                print(0, p, u_range[0], u_range[-1])
            elif i != 0 and rim_tmp == e_curve.Value(u1):
                u_range = np.linspace(u1, u0, 50)
                rim_tmp = e_curve.Value(u0)
                p = e_curve.Value(u1)
                p.Transform(set_trf(self.axs, gp_Ax3()))
                data = [p.X(), p.Y()]
                print(1, p, u_range[0], u_range[-1])
            else:
                u_range = np.linspace(u0, u1, 50)
                rim_tmp = e_curve.Value(u1)
                p = e_curve.Value(u0)
                p.Transform(set_trf(self.axs, gp_Ax3()))
                data = [p.X(), p.Y()]
                print(2, p, u_range[0], u_range[-1])
            fp.write(''.join([float_to_string(val) for val in data]) + '\n')
            for u in u_range[1:]:
                p = e_curve.Value(u)
                p.Transform(set_trf(self.axs, gp_Ax3()))
                data = [p.X(), p.Y()]
                fp.write(''.join([float_to_string(val)
                                  for val in data]) + '\n')
        fp.close()

    def export_sfc1_axs(self, sfcfile="m2_mat.sfc", name="M2 Mat"):
        surf = BRep_Tool.Surface(self.surf)

        trf = set_trf(self.axs, gp_Ax3())
        xy0 = dispocc.proj_pnt_pln(self, surf.Value(0, 0), self.pln, self.axs)
        xy1 = dispocc.proj_pnt_pln(self, surf.Value(1, 1), self.pln, self.axs)
        xy0.Transform(trf)
        xy1.Transform(trf)

        m2_trf = set_trf(gp_Ax3(), self.axs)
        m2_pln = BRep_Tool.Surface(self.pln)
        for px in np.linspace(-100, 100, 10):
            for py in np.linspace(-100, 100, 10):
                p0 = gp_Pnt(px, py, 0).Transformed(m2_trf)
                p1 = obj.proj_pnt_pln(p0, self.surf, self.axs)

        #ix0, ix1 = m2.surf_pts.LowerRow(), m2.surf_pts.UpperRow()
        #iy0, iy1 = m2.surf_pts.LowerCol(), m2.surf_pts.UpperCol()
        #xy0 = m2.surf_pts.Value(ix0, iy0).Transformed(trf)
        #xy1 = m2.surf_pts.Value(ix1, iy1).Transformed(trf)
        nx, ny = 200, 200
        xs, xe = xy0.X(), xy1.X()
        ys, ye = xy0.Y(), xy1.Y()
        fp = open(sfcfile, "w")
        fp.write(" {} \n".format(name))
        fp.write(" {:.2e} {:.2e} {:.2e} {:.2e}\n".format(xs, ys, xe, ye))
        fp.write(" {:d} {:d}\n".format(nx, ny))
        for ix in np.linspace(0, 1, nx):
            for iy in np.linspace(0, 1, ny):
                p0 = surf.Value(ix, iy)
                p1 = dispocc.proj_pnt_pln(self, p0, self.pln, self.axs)
                pz = p1.Transformed(trf)
                z = p0.Distance(p1)
                fp.write(" {:.5e} ".format(z))
            fp.write("\n")
        fp.close()
        print(xy0)

    def export_sfc2_axs(self, nxy=[200, 200], rx=[-250, 250], ry=[-250, 250], sfcfile="m2_mat.sfc"):
        trsf = set_trf(gp_Ax3(), self.axs)
        nx, ny = nxy
        xs, xe = rx
        ys, ye = ry
        plnx = np.linspace(xs, xe, nx)
        plny = np.linspace(ys, ye, ny)
        mesh = np.meshgrid(plnx, plny)
        data = np.zeros_like(mesh[0])
        for (ix, iy), x in np.ndenumerate(data):
            px, py = mesh[0][ix, iy], mesh[1][ix, iy]
            p0 = gp_Pnt(px, py, 0).Transformed(trsf)
            p1 = dispocc.proj_pnt_pln(self, p0, self.surf, self.axs)
            z = p0.Distance(p1)
            data[ix, iy] = z

            txt = "\r"
            txt += "{:d}, {:d} / {:d}, {:d}".format(ix, iy, ny, nx)
            sys.stdout.write(txt)
            sys.stdout.flush()
        print()
        grasp_sfc(mesh, data, sfcfile)

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
    print(argc, opt)

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
