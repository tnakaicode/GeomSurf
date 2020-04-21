import numpy as np

from OCC.Core.gp import gp_Pln, gp_Dir, gp_Pnt
from OCC.Core.BRep import BRep_Builder
from OCC.Core.BRepTools import breptools_Read
from OCC.Core.TopoDS import TopoDS_Shape
from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRepBndLib import brepbndlib_Add
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Section
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCC.Extend.ShapeFactory import get_boundingbox, make_box
from OCC.Extend.DataExchange import write_step_file

from base import plotocc


def get_brep():
    cylinder_head = TopoDS_Shape()
    builder = BRep_Builder()
    breptools_Read(cylinder_head, './core_example/cylinder_head.brep', builder)
    return cylinder_head


if __name__ == '__main__':
    shp = get_brep()
    xyz_min_max = get_boundingbox(shp)
    p1 = gp_Pnt(*xyz_min_max[0:3])
    p2 = gp_Pnt(*xyz_min_max[3:])
    box = make_box(p1, p2)

    obj = plotocc()
    obj.display.DisplayShape(box, transparency=0.9)
    obj.display.DisplayShape(shp)
    obj.show_axs_pln(scale=75)
    obj.show()
