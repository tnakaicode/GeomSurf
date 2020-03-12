import numpy as np
import matplotlib.pyplot as plt
import sys

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Lin, gp_Sphere, gp_Pln
from OCC.Core.BRep import BRep_Tool, BRep_Builder, BRep_PointsOnSurface
from OCC.Core.BRepFill import BRepFill_Filling
from OCC.Core.BRepCheck import BRepCheck_Face, BRepCheck_Analyzer
from OCC.Core.BRepTools import BRepTools_MapOfVertexPnt2d
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse, BRepAlgoAPI_Section
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeWire
from OCC.Core.BRepProj import BRepProj_Projection
from OCC.Core.Geom import Geom_ToroidalSurface
from OCC.Core.Geom import Geom_Line
from OCC.Core.GeomAbs import GeomAbs_G1, GeomAbs_G2
from OCC.Core.GeomAbs import GeomAbs_C0, GeomAbs_C1, GeomAbs_C2, GeomAbs_C3
from OCC.Core.GeomAPI import GeomAPI_IntCS, GeomAPI_IntSS
from OCC.Core.GeomLProp import GeomLProp_SurfaceTool
from OCC.Core.TopoDS import TopoDS_Compound, TopoDS_Builder, TopoDS_Face, TopoDS_Shell
from OCC.Extend.DataExchange import read_step_file, write_step_file
from OCCUtils.Construct import make_n_sided, make_n_sections
from OCCUtils.Construct import make_edge, make_polygon, make_vertex
from OCCUtils.Construct import vec_to_dir, dir_to_vec
from OCCUtils.Construct import point_to_vector, vector_to_point
from OCCUtils.Topology import Topo, dumpTopology

from base import plotocc, gen_ellipsoid, set_loc, spl_face


def line_from_axs(axs=gp_Ax3(), length=100):
    vec = point_to_vector(axs.Location()) + \
        dir_to_vec(axs.Direction()) * length
    return make_edge(axs.Location(), vector_to_point(vec))


class HexPlane (plotocc):

    def __init__(self):
        plotocc.__init__(self)
        self.compound = TopoDS_Compound()
        self.builder = BRep_Builder()
        self.builder.MakeCompound(self.compound)

        self.beam = gp_Ax3()
        self.beam.SetLocation(gp_Pnt(0.5, 0.5, 0.0))
        self.beam.SetDirection(gp_Dir(0.0, 0.5, 1.0))
        self.beam_line = line_from_axs(self.beam, length=20)
        self.builder.Add(self.compound, self.beam_line)

        ax = gp_Ax3(gp_Pnt(0, 0, 10), gp_Dir(0, 0, -1))
        px = np.linspace(-1, 1, 10) * 10
        py = np.linspace(-1, 1, 10) * 10
        mesh = np.meshgrid(px, py)
        surf = mesh[0]**2 / 100 + mesh[1]**2 / 150
        self.surf = spl_face(*mesh, surf, ax)
        self.surf_bound = self.make_PolySurf(radi=5, axs=ax)

        self.beam_glin = Geom_Line(self.beam.Location(), self.beam.Direction())
        self.ics = GeomAPI_IntCS(self.beam_glin, BRep_Tool.Surface(self.surf))
        print(self.ics.NbPoints())
        # print(self.ics.Point(1))

        self.ics = GeomAPI_IntCS(
            self.beam_glin, BRep_Tool.Surface(self.surf_bound))
        print(self.ics.NbPoints())

        #self.display.DisplayShape(self.surf, transparency=0.7)
        self.display.DisplayShape(self.surf_bound, transparency=0.7)
        self.plns = TopoDS_Shell()
        self.builder.MakeShell(self.plns)
        for ix in np.linspace(0, 1, 5):
            for iy in np.linspace(0, 1, 5):
                p1, vx, vy = gp_Pnt(), gp_Vec(), gp_Vec()
                GeomLProp_SurfaceTool.D1(
                    BRep_Tool.Surface(self.surf), ix, iy, p1, vx, vy)
                vz = vx.Crossed(vy)
                axs = gp_Ax3(p1, vec_to_dir(vz), vec_to_dir(vx))
                pln = self.make_PolyPlane(axs=axs, radi=2.5, shft=15.0)
                print(pln)

                self.builder.Add(self.compound, make_vertex(p1))
                self.builder.Add(self.plns, pln)
        self.builder.Add(self.compound, self.plns)

        for face in Topo(self.plns).faces():
            self.ics.Perform(self.beam_glin, BRep_Tool.Surface(face))
            uvw = self.ics.Parameters(1)
            u, v, w = uvw
            p1, vx, vy = gp_Pnt(), gp_Vec(), gp_Vec()
            GeomLProp_SurfaceTool.D1(
                BRep_Tool.Surface(face), u, v, p1, vx, vy)
            vz = vx.Crossed(vy)
            if u > 0 and v > 0:
                print(u, v)
                print(p1)
                print(self.ics.Point(1))
                self.display.DisplayShape(p1)
                self.display.DisplayShape(face, color="BLUE")
            else:
                print(u, v)

    def make_PolyPlane(self, num=6, radi=1.0, shft=0.0, axs=gp_Ax3()):
        lxy = radi - 0.1
        pnts = []
        angl = 360 / num
        for i in range(num):
            thet = np.deg2rad(i * angl) + np.deg2rad(shft)
            x, y = radi * np.sin(thet), radi * np.cos(thet)
            pnts.append(gp_Pnt(x, y, 0))
        pnts.append(pnts[0])
        poly = make_polygon(pnts)
        brep = BRepBuilderAPI_MakeFace(gp_Pln(), poly)
        brep.Add(poly)
        face = brep.Face()
        face.Location(set_loc(gp_Ax3(), axs))
        return face

    def make_PolySurf(self, num=6, radi=1.0, shft=0.0, axs=gp_Ax3()):
        lxy = radi - 0.1
        pnts = []
        angl = 360 / num
        for i in range(num):
            thet = np.deg2rad(i * angl) + np.deg2rad(shft)
            x, y = radi * np.sin(thet), radi * np.cos(thet)
            pnts.append(gp_Pnt(x, y, 0))
        pnts.append(pnts[0])
        poly = make_polygon(pnts)
        proj = BRepProj_Projection(poly, self.surf, gp_Pnt(0, 0, 10))
        # print(proj.Current())
        brep = BRepBuilderAPI_MakeFace(self.surf, poly)
        # brep.Add(poly)
        face = brep.Face()
        face.Location(set_loc(gp_Ax3(), axs))
        return face

    def export_file(self):
        write_step_file(self.compound, self.tempname + ".stp")

    def display_Plane(self):
        self.display.DisplayShape(self.beam_line)
        self.display.DisplayShape(self.compound)
        self.show_axs_pln(scale=1.0)
        self.show()


if __name__ == '__main__':
    obj = HexPlane()
    obj.export_file()
    obj.display_Plane()
