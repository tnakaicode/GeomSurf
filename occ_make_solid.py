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
from OCC.Core.TopoDS import topods, TopoDS_Shape, TopoDS_Compound
from OCC.Core.TopoDS import TopoDS_Wire, TopoDS_Vertex
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Core.TopTools import TopTools_ListOfShape
from OCC.Core.GeomAbs import GeomAbs_C2
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Cut
from OCC.Core.BRepFilletAPI import BRepFilletAPI_MakeFillet
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_ThruSections, BRepOffsetAPI_MakePipeShell
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeSolid
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCCUtils.Common import vertex2pnt
from OCCUtils.Construct import ShapeToTopology
from OCCUtils.Construct import make_box, make_line, make_wire, make_edge, make_circle
from OCCUtils.Construct import make_plane, make_polygon, make_face
from OCCUtils.Construct import make_pipe, make_coons
from OCCUtils.Construct import make_solid
from OCCUtils.Construct import point_to_vector, vector_to_point
from OCCUtils.Construct import dir_to_vec, vec_to_dir


def make_loft(
    elements,
    ruled=True,
    tolerance=1.0E-6,
    continuity=GeomAbs_C2,
    check_compatibility=True,
):

    sections = BRepOffsetAPI_ThruSections(False, ruled, tolerance)
    for i in elements:
        if isinstance(i, TopoDS_Wire):
            sections.AddWire(i)
        elif isinstance(i, TopoDS_Vertex):
            sections.AddVertex(i)
        else:
            raise TypeError(
                "elements is a list of TopoDS_Wire or TopoDS_Vertex, found a %s fool"
                % i.__class__
            )

    sections.CheckCompatibility(check_compatibility)
    sections.SetContinuity(continuity)
    sections.Build()

    return make_solid(sections.Shape())


if __name__ == '__main__':
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument("--pxyz", dest="pxyz",
                        default=[0.0, 0.0, 0.0], type=float, nargs=3)
    opt = parser.parse_args()
    print(opt)

    obj = dispocc(touch=True)
    axs = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1))
    box = make_box(axs.Ax2(), 3000, 3000, 5000)

    shp1 = make_box(gp_Ax2(gp_Pnt(100, 100, 200),
                           gp_Dir(0.1, 0.1, 1.0)),
                    100, 100, 500.0)
    box = BRepAlgoAPI_Cut(box, shp1).Shape()

    # wire1 = make_wire(make_circle(gp_Pnt(1500, 1500, 3000), 500))
    # wire2 = make_wire(make_circle(gp_Pnt(1500, 1500, 3500), 500))
    # shp2 = make_loft([wire1, wire2], ruled=True)
    # print(shp2)
    # box = BRepAlgoAPI_Cut(box, shp2).Shape()

    shp2 = make_box(gp_Ax2(gp_Pnt(1000, 2000, 1500),
                           gp_Dir(-0.1, 0.5, 1.0)),
                    500, 500, 1000.0)
    edges = [e for e in TopologyExplorer(shp2).edges()]
    fill = BRepFilletAPI_MakeFillet(shp2)
    fill.Add(edges[0])
    fill.Add(edges[2])
    fill.Add(edges[3])

    for i in range(1, fill.NbContours() + 1):
        length = fill.Length(i)
        radius = 0.15 * length
        fill.SetRadius(radius, i, 1)
    shp2_fill = fill.Shape()
    box = BRepAlgoAPI_Cut(box, shp2_fill).Shape()

    wire2 = make_wire(make_edge(gp_Circ(gp_Ax2(gp_Pnt(1500, 1500, 2500),
                                               gp_Dir(0.2, 0.4, 0.5)),
                                        100),
                                -np.pi / 2, np.pi / 3)
                      )
    shell = BRepOffsetAPI_MakePipeShell(wire2)
    shell.Add(make_polygon([gp_Pnt(100, 100, 800),
                            gp_Pnt(300, 400, 800),
                            gp_Pnt(200, 800, 700)], closed=True))
    shell.Build()
    shell.MakeSolid()
    shp3 = shell.Shape()
    box = BRepAlgoAPI_Cut(box, shp3).Shape()

    obj.selected_shape = [box]
    obj.export_stp_selected()

    obj.display.DisplayShape(box, transparency=0.9)
    obj.display.DisplayShape(axs.Location())
    obj.ShowOCC()
