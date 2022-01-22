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
from OCC.Extend.ShapeFactory import get_boundingbox
from OCC.Extend.DataExchange import write_step_file
from OCCUtils.Construct import make_box

from src.base_occ import dispocc


def vectorized_slicer(li):
    # Create Plane defined by a point and the perpendicular direction
    z_values, shape = li
    _slices = []
    for z in z_values:
        # print 'slicing index:', z, 'sliced by process:', os.getpid()
        plane = gp_Pln(gp_Pnt(0., 0., z), gp_Dir(0., 0.01, 1.))
        face = BRepBuilderAPI_MakeFace(plane).Shape()
        # Computes Shape/Plane intersection
        section = BRepAlgoAPI_Section(shape, face)
        section.Build()
        if section.IsDone():
            _slices.append(section.Shape())
    return _slices


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

    obj = dispocc()
    obj.display.DisplayShape(box, transparency=0.9)
    obj.display.DisplayShape(shp)

    z_delta = abs(xyz_min_max[2] - xyz_min_max[5])
    for z in np.linspace(xyz_min_max[2], xyz_min_max[5], 5):
        print(z)
        plane = gp_Pln(gp_Pnt(0., 0., z), gp_Dir(0., 0.0, 1.))
        face = BRepBuilderAPI_MakeFace(plane).Shape()
        # Computes Shape/Plane intersection
        section = BRepAlgoAPI_Section(shp, face)
        section.Build()
        if section.IsDone():
            obj.display.DisplayShape(section.Shape(), color="BLUE1")
            obj.export_stp(section.Shape())

    obj.show_axs_pln(scale=75)
    obj.show_occ()
