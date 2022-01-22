import numpy as np
import matplotlib.pyplot as plt
import sys

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Lin, gp_Sphere
from OCC.Core.BRep import BRep_Tool, BRep_Builder
from OCC.Core.Geom import Geom_ToroidalSurface
from OCC.Core.Geom import Geom_Line
from OCC.Core.GeomAPI import GeomAPI_IntCS
from OCC.Core.GeomLProp import GeomLProp_SurfaceTool
from OCC.Core.TopoDS import TopoDS_Compound
from OCC.Extend.DataExchange import read_step_file, write_step_file
from OCCUtils.Construct import make_n_sided, make_n_sections
from OCCUtils.Construct import make_edge, make_polygon
from OCCUtils.Construct import vec_to_dir, dir_to_vec
from OCCUtils.Topology import Topo

from src.base import dispocc, gen_ellipsoid


class DBall (plotocc):

    def __init__(self):
        plotocc.__init__(self)
        self.b1 = self.get_face(gen_ellipsoid(rxyz=[100, 100, 105]))
        self.b2 = self.get_face(gen_ellipsoid(rxyz=[210, 210, 210]))
        self.beam = gp_Ax3(gp_Pnt(0, 150, 0), gp_Dir(1, 1.5, 0))
        print(self.b1)

        h_surf = BRep_Tool.Surface(self.b2)
        ray = Geom_Line(self.beam.Axis())
        self.RayTrace = GeomAPI_IntCS(ray, h_surf)
        print(self.RayTrace.NbPoints())
        self.num = 0
        self.pts = [self.beam.Location()]
        for i in range(5):
            self.beam = self.reflect_b1(num=1)
            self.beam = self.reflect_b1(num=2)
            self.beam = self.reflect_b2(num=1)
            self.beam = self.reflect_b2(num=2)

    def beam_run(self):
        for i in range(3):
            self.beam = self.reflect_b1()
            self.beam = self.reflect_b2()

    def get_face(self, sol):
        top_api = Topo(sol)
        print(top_api.number_of_faces())
        for face in top_api.faces():
            sol_face = face
        return sol_face

    def reflect_b1(self, num=1):
        h_surf = BRep_Tool.Surface(self.b1)
        ray = Geom_Line(self.beam.Axis())
        self.RayTrace.Perform(ray, h_surf)
        if self.RayTrace.NbPoints() == 0:
            beam = self.beam
        else:
            self.num += 1
            uvw = self.RayTrace.Parameters(num)
            u, v, w = uvw
            p1, vx, vy = gp_Pnt(), gp_Vec(), gp_Vec()
            GeomLProp_SurfaceTool.D1(h_surf, u, v, p1, vx, vy)
            vz = vx.Crossed(vy)
            norm = gp_Ax3(p1, vec_to_dir(vz), vec_to_dir(vx))
            self.show_axs_pln(norm, scale=10)
            beam = self.beam
            beam.SetLocation(p1)
            beam.SetDirection(beam.Direction().Reversed())
            beam.Mirror(norm.Ax2())
            print(self.num, self.b1, p1)
            self.pts.append(p1)
            # self.beam.XReverse()
            # self.beam.Mirror(norm.Ax2())
        return beam

    def reflect_b2(self, num=1):
        h_surf = BRep_Tool.Surface(self.b2)
        ray = Geom_Line(self.beam.Axis())
        self.RayTrace.Perform(ray, h_surf)
        if self.RayTrace.NbPoints() == 0:
            beam = self.beam
        else:
            self.num += 1
            uvw = self.RayTrace.Parameters(num)
            u, v, w = uvw
            p1, vx, vy = gp_Pnt(), gp_Vec(), gp_Vec()
            GeomLProp_SurfaceTool.D1(h_surf, u, v, p1, vx, vy)
            vz = vx.Crossed(vy).Reversed()
            norm = gp_Ax3(p1, vec_to_dir(vz), vec_to_dir(vx))
            self.show_axs_pln(norm, scale=10)
            beam = self.beam
            beam.SetLocation(p1)
            beam.SetDirection(beam.Direction().Reversed())
            beam.Mirror(norm.Axis())
            print(self.num, self.b2, p1)
            self.pts.append(p1)
            # self.beam.XReverse()
            # self.beam.Mirror(norm.Ax2())
        return beam

    # def reflect(self, p0, v0, face):
    #    h_surf = BRep_Tool.Surface(face)
    #    ray = Geom_Line(gp_Lin(p0, vec_to_dir(v0)))
    #    uvw = GeomAPI_IntCS(ray.GetHandle(), h_surf).Parameters(1)
    #    u, v, w = uvw
    #    p1, vx, vy = gp_Pnt(), gp_Vec(), gp_Vec()
    #    GeomLProp_SurfaceTool.D1(h_surf, u, v, p1, vx, vy)
    #    vz = vx.Crossed(vy)
    #    vx.Normalize()
    #    vy.Normalize()
    #    vz.Normalize()
    #    v1 = v0.Mirrored(gp_Ax2(p1, vec_to_dir(vz)))
    #    return p1, v1

    def export_file(self):
        builder = BRep_Builder()
        compound = TopoDS_Compound()
        builder.MakeCompound(compound)
        builder.Add(compound, self.b1)
        builder.Add(compound, self.b2)
        builder.Add(compound, make_polygon(self.pts))
        self.export_stp(compound)

    def display_ball(self):
        #self.show_vec(self.beam, scale=50)
        self.display.DisplayShape(self.b1, transparency=0.7, color="RED")
        self.display.DisplayShape(self.b2, transparency=0.7, color="BLUE1")
        self.display.DisplayShape(self.beam.Location())
        self.show_axs_pln(scale=100)
        self.show_axs_pln(self.beam, scale=10)
        self.show_pts(self.pts)
        self.show()
        self.export_file()


if __name__ == '__main__':
    obj = DBall()
    obj.display_ball()
