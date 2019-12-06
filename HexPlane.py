import numpy as np
import matplotlib.pyplot as plt
import sys

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Lin, gp_Sphere
from OCC.Core.BRep import BRep_Tool, BRep_Builder
from OCC.Core.BRepFill import BRepFill_Filling
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse, BRepAlgoAPI_Section
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeWire
from OCC.Core.Geom import Geom_ToroidalSurface
from OCC.Core.Geom import Geom_Line
from OCC.Core.GeomAbs import GeomAbs_G1, GeomAbs_G2
from OCC.Core.GeomAbs import GeomAbs_C0, GeomAbs_C1, GeomAbs_C2, GeomAbs_C3
from OCC.Core.GeomAPI import GeomAPI_IntCS
from OCC.Core.GeomLProp import GeomLProp_SurfaceTool
from OCC.Core.TopoDS import TopoDS_Compound, TopoDS_Builder
from OCC.Extend.DataExchange import read_step_file, write_step_file
from OCCUtils.Construct import make_n_sided, make_n_sections
from OCCUtils.Construct import make_edge, make_polygon, make_vertex
from OCCUtils.Construct import vec_to_dir, dir_to_vec
from OCCUtils.Topology import Topo

from base import plotocc, gen_ellipsoid, set_loc, spl_face

# https://www.opencascade.com/doc/occt-7.4.0/refman/html/class_b_rep_fill___filling.html
# N-Side Filling This algorithm avoids to build a face from:
#   a set of edges defining the bounds of the face and some constraints the surface support has to satisfy
#   a set of edges and points defining some constraints the support surface has to satisfy
#       an initial surface to deform for satisfying the constraints
#
#   a set of parameters to control the constraints.
#   The support surface of the face is computed by deformation of the initial surface in order to satisfy the given constraints.
#   The set of bounding edges defines the wire of the face.
#
#   If no initial surface is given, the algorithm computes it automatically.
#   If the set of edges is not connected (Free constraint) missing edges are automatically computed.
#
#   Limitations:
#
#   If some constraints are not compatible
#       The algorithm does not take them into account.
#       So the constraints will not be satisfyed in an area containing the incompatibilitries.
#       he constraints defining the bound of the face have to be entered in order to have a continuous wire.
#
#   Other Applications:
#   Deformation of a face to satisfy internal constraints
#   Deformation of a face to improve Gi continuity with connected faces


class HexPlane (plotocc):

    def __init__(self):
        plotocc.__init__(self)
        self.compound = TopoDS_Compound()
        self.builder = BRep_Builder()
        self.builder.MakeCompound(self.compound)

        ax = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1))
        px = np.linspace(-1, 1, 10) * 10
        py = np.linspace(-1, 1, 10) * 10
        mesh = np.meshgrid(px, py)
        surf = mesh[0]**2 / 100 + mesh[1]**2 / 150
        self.surf = spl_face(*mesh, surf, ax)

        self.display.DisplayShape(self.surf, transparency=0.7)
        for ix in np.linspace(0, 1, 5):
            for iy in np.linspace(0, 1, 5):
                p1, vx, vy = gp_Pnt(), gp_Vec(), gp_Vec()
                GeomLProp_SurfaceTool.D1(
                    BRep_Tool.Surface(self.surf), ix, iy, p1, vx, vy)
                vz = vx.Crossed(vy)
                axs = gp_Ax3(p1, vec_to_dir(vz), vec_to_dir(vx))
                pln = self.make_PolyPlane(axs=axs, radi=2.5, shft=15.0)
                print(ix, iy, pln)
                self.display.DisplayShape(p1)
                self.display.DisplayShape(pln)
                self.builder.Add(pln, make_vertex(p1))
                self.builder.Add(self.compound, pln)

    def make_PolyPlane(self, num=6, radi=1.0, shft=0.0, axs=gp_Ax3()):
        pnts = []
        angl = 360 / num
        for i in range(num):
            thet = np.deg2rad(i * angl) + np.deg2rad(shft)
            x, y = radi * np.sin(thet), radi * np.cos(thet)
            pnts.append(gp_Pnt(x, y, 0))
        pnts.append(pnts[0])
        poly = make_polygon(pnts)

        # n_sided = BRepFill_Filling()
        # n_sided.SetApproxParam()
        # n_sided.SetResolParam()
        # n_sided.SetConstrParam()
        # for i, pnt in enumerate(pnts[:-1]):
        #    i0, i1 = i, i + 1
        #    n_sided.Add(pnt)
        #    n_sided.Add(make_edge(pnts[i0], pnts[i1]), GeomAbs_C0)
        # n_sided.Add(gp_Pnt(0, 0, 1))
        # n_sided.Build()
        # face = n_sided.Face()
        # print(n_sided.G0Error())

        face = BRepBuilderAPI_MakeFace(poly).Face()
        face.Location(set_loc(gp_Ax3(), axs))
        return face

    def export_file(self):
        write_step_file(self.compound, "./tmp/HexPlane.stp")

    def display_Plane(self):
        self.show_axs_pln(scale=1.0)
        self.show()


if __name__ == '__main__':
    obj = HexPlane()
    obj.export_file()
    obj.display_Plane()
