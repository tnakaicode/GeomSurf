import numpy as np
from numpy import sin, cos, pi
from skimage import measure

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pnt, gp_Ax3, gp_Dir
from OCC.Core.GeomAbs import GeomAbs_C0
from OCC.Core.GeomAbs import GeomAbs_Intersection, GeomAbs_Arc
from OCC.Core.BRepFill import BRepFill_Filling
from OCC.Core.BRepFill import BRepFill_CurveConstraint
from OCC.Core.BRepOffset import BRepOffset_MakeOffset, BRepOffset_Skin, BRepOffset_Interval
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakePolygon
from OCC.Extend.TopologyUtils import TopologyExplorer, WireExplorer
from OCC.Extend.ShapeFactory import make_face, make_vertex
from OCC.Extend.DataExchange import read_step_file, write_step_file
from OCCUtils.Construct import make_n_sided, make_n_sections
from OCCUtils.Construct import make_edge, make_polygon
from OCCUtils.Topology import Topo

from src.base_occ import dispocc


class MakeSoild (dispocc):

    def __init__(self):
        dispocc.__init__(self)
        axs = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1))
        shp = self.make_PolyWire(num=15, skin=None, axs=axs)
        self.display.DisplayShape(shp)
        self.export_stp(shp)

        axs = gp_Ax3(gp_Pnt(0, 0, 1), gp_Dir(0, 0, 1))
        shp = self.make_PolyWire(num=15, skin=0, axs=axs)
        self.display.DisplayShape(shp)
        self.export_stp(shp)

        axs = gp_Ax3(gp_Pnt(0, 0, 3), gp_Dir(0, 0, 1))
        shp = self.make_StarWire(num=15, skin=1, axs=axs, radi=[1.1, 1.0])
        self.display.DisplayShape(shp)
        self.export_stp(shp)

        axs = gp_Ax3(gp_Pnt(0, 0, 5), gp_Dir(0, 1, 1))
        shp = self.make_StarWire(num=15, skin=1, axs=axs, radi=[0.9, 1.0])
        self.display.DisplayShape(shp)
        self.export_stp(shp)


if __name__ == '__main__':
    obj = MakeSoild()
    obj.show_axs_pln(scale=1.0)
    obj.ShowOCC()
