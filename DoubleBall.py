import numpy as np
import matplotlib.pyplot as plt
import sys

from OCC.Display.SimpleGui import init_display
from OCC.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.gp import gp_Lin, gp_Sphere
from OCC.BRep import BRep_Tool
from OCC.Geom import Geom_ToroidalSurface
from OCC.Geom import Geom_Line
from OCC.GeomAPI import GeomAPI_IntCS
from OCC.TopoDS import TopoDS_Shape, TopoDS_Compound
from OCCUtils.Construct import make_n_sided, make_n_sections
from OCCUtils.Construct import make_edge
from OCCUtils.Topology import Topo

from base import plotocc, gen_ellipsoid


class DBall (plotocc):

    def __init__(self):
        plotocc.__init__(self)
        self.b1 = self.get_face(gen_ellipsoid(rxyz=[100, 100, 105]))
        self.b2 = self.get_face(gen_ellipsoid(rxyz=[210, 210, 210]))
        self.beam = gp_Ax3(gp_Pnt(0, 150, 0), gp_Dir(1, 0, 0))
        print(self.b1)

        h_surf = BRep_Tool.Surface(self.b2)
        ray = Geom_Line(self.beam.Axis())
        self.RayTrace = GeomAPI_IntCS(ray.GetHandle(), h_surf)

    def get_face(self, sol):
        top_api = Topo(sol)
        print(top_api.number_of_faces())
        for face in top_api.faces():
            sol_face = face
        return sol_face

    def reflect_b1(self):
        print(self.b1)
        h_surf = BRep_Tool.Surface(self.b1)
        ray = Geom_Line(self.beam.Axis())
        self.RayTrace.Perform(ray, h_surf)

    def reflect_b2(self):
        print(self.b2)
        h_surf = BRep_Tool.Surface(self.b2)
        ray = Geom_Line(self.beam.Axis())
        self.RayTrace.Perform(ray, h_surf)

    #def reflect(self, p0, v0, face):
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

    def display_ball(self):
        self.show_vec(self.beam, scale=50)
        self.display.DisplayShape(self.b1, transparency=0.7, color="RED")
        self.display.DisplayShape(self.b2, transparency=0.7, color="BLUE")
        self.show_axs_pln(scale=100)
        self.show()


if __name__ == '__main__':
    obj = DBall()
    obj.display_ball()
