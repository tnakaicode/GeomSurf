import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import argparse
from linecache import getline, clearcache, updatecache

basename = os.path.dirname(__file__)

sys.path.append(os.path.join("./"))
from base_occ import dispocc

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Circ
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Compound
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_ThruSections
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
from OCC.Core.BRepLProp import BRepLProp_SLProps
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCCUtils.Common import vertex2pnt
from OCCUtils.Construct import make_box, make_line, make_wire, make_edge
from OCCUtils.Construct import make_plane, make_polygon
from OCCUtils.Construct import point_to_vector, vector_to_point
from OCCUtils.Construct import dir_to_vec, vec_to_dir


class BeamReflect (dispocc):

    def __init__(self, temp=True, disp=True, touch=False):
        super().__init__(temp, disp, touch)

        r0 = 10.0
        r1 = 11.0
        r2 = 15.0
        r3 = 1.0
        l0 = 0.0
        l1 = 100.0
        l2 = 150.0
        l3 = 200.0
        ax0 = gp_Ax3(gp_Pnt(0, 0, l0), gp_Dir(0, 0, 1))
        ax1 = gp_Ax3(gp_Pnt(0, 0, l1), gp_Dir(0, 0, 1))
        ax2 = gp_Ax3(gp_Pnt(0, 0, l2), gp_Dir(0, 1, 1))
        ax3 = gp_Ax3(gp_Pnt(0, 0, l3), gp_Dir(1, 0, 1))
        cr0 = make_edge(gp_Circ(ax0.Ax2(), r0))
        cr1 = self.make_ellips(ax1, [r1, 20], 30.0, 200, None)
        cr2 = self.make_ellips(ax2, [r2, 16], -5.0, 200, None)
        cr3 = self.make_ellips(ax3, [r3, 1], -15.0, 200, None)

        api = BRepOffsetAPI_ThruSections()
        api.AddWire(make_wire(cr0))
        api.AddWire(make_wire(cr1))
        api.AddWire(make_wire(cr2))
        api.AddWire(make_wire(cr3))
        api.Build()
        self.shll = api.Shape()
        self.display.DisplayShape(self.shll, transparency=0.9)
        self.face = [f for f in TopologyExplorer(self.shll).faces()][0]
        self.surf = BRepAdaptor_Surface(self.face)
        self.us, self.ue = self.surf.FirstUParameter(), self.surf.LastUParameter()
        self.vs, self.ve = self.surf.FirstVParameter(), self.surf.LastVParameter()
        print(self.surf.FirstUParameter(), self.surf.LastUParameter())
        print(self.surf.FirstVParameter(), self.surf.LastVParameter())

        self.display.DisplayShape(cr0)
        self.display.DisplayShape(cr1)
        self.display.DisplayShape(cr2)
        self.display.DisplayShape(self.face)

    def calc_uv(self, u=0.1, v=0.1):
        return self.us + u * (self.ue - self.us), self.vs + v * (self.ve - self.vs)

    def display_uv(self, u=0, v=0):
        prop = BRepLProp_SLProps(self.surf, *self.calc_uv(u, v), 2, 1.0E-9)
        p1, vx, vy = prop.Value(), prop.D1U(), prop.D1V()
        print(u, v, p1)
        self.display.DisplayShape(p1)


if __name__ == '__main__':
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument("--pxyz", dest="pxyz",
                        default=[0.0, 0.0, 0.0], type=float, nargs=3)
    opt = parser.parse_args()
    print(opt)

    obj = BeamReflect(touch=True)
    obj.display_uv(0.0, 0.0)
    obj.display_uv(0.0, 0.1)
    obj.display_uv(0.0, 0.2)
    obj.display_uv(0.0, 0.3)
    obj.display_uv(0.1, 0.0)
    obj.show_axs_pln(scale=10)

    beam0 = gp_Ax3(gp_Pnt(0, 0, 0),
                   gp_Dir(1, 1, 1))
    pts = [beam0.Location()]
    for i in range(20):
        beam0 = obj.reflect_beam(obj.shll, beam0=beam0, tr=0)
        pts.append(beam0.Location())
    obj.display.DisplayShape(make_polygon(pts))

    obj.ShowOCC()
