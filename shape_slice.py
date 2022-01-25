"""
This script demonstrates loading 2d cross-section from iges which can be crated from 2d cad program like autocad.
"""

import OCC.Extend.DataExchange as OCC_Extend_Data_Exchange
import OCC.Extend.TopologyUtils as OCC_Extend_Topo_Utils
import OCC.Core.BRepBndLib as OCC_BREP_Bounding_Lib
import OCC.Core.Bnd as OCC_Bounding
import OCC.Core.BRepBuilderAPI as OCC_Brep_Builder
import OCC.Core.BRepAdaptor as OCC_Brep_Adaptor


def get_direction_of_wire(wire):
    """
    A method to detect 2d wire orientation is clockwise or counterclockwise.
    For counterclockwise it should be +z.
    For clockwise it should be -z.
    :param wire: Wire to be investigating
    :return: direction of the normal of wire
    """
    make_face_wire = OCC_Brep_Builder.BRepBuilderAPI_MakeFace(wire)
    make_face_wire.Build()
    surface_wire = OCC_Brep_Adaptor.BRepAdaptor_Surface(make_face_wire.Shape())
    plane_wire = surface_wire.Plane()
    plane_wire_normal = plane_wire.Axis().Direction()
    return plane_wire_normal


def load_section_from_iges(file):
    """
    A method to create a face entity to be extruded to solid or mesh from iges file

    In 2D drawings all wires should be joint and they should be polyline or composite curve.

    Method firstly calculate bbox of all wires and detect outer wire which have the largest diagonal of bbox.
    Other wires are inner wires.

    After that, orientation of the outer wire should be counterclockwise so check the normal reverse it
    if normal is not upward.
    Orientation of inner wire should be clockwise so check the normals and reverse wire if normal is not downward.

    :param file: Iges file to load
    :return: closed face with inner holes
    """
    data = OCC_Extend_Data_Exchange.read_iges_file(file)
    topo_explorer = OCC_Extend_Topo_Utils.TopologyExplorer(data)
    wires = list(topo_explorer.wires())

    wire_dict = {'outer': None, 'inner': None}

    bbox_max = 0.0
    index_max = 0
    for i, wire in enumerate(wires):
        bbox = OCC_Bounding.Bnd_Box()
        bbox.SetGap(1e-6)
        OCC_BREP_Bounding_Lib.brepbndlib_Add(wire, bbox)
        diagonal = bbox.SquareExtent()
        if diagonal > bbox_max:
            bbox_max = diagonal
            index_max = i
    wire_dict['outer'] = wires[index_max]
    wires.pop(index_max)
    wire_dict['inner'] = wires

    plane_outer_normal = get_direction_of_wire(wire_dict['outer'])

    if plane_outer_normal.Z() < 0:
        wire_dict['outer'].Reverse()
    make_face = OCC_Brep_Builder.BRepBuilderAPI_MakeFace(wire_dict['outer'])

    for wire_inner in wire_dict['inner']:
        plane_inner_normal = get_direction_of_wire(wire_inner)
        if plane_inner_normal.Z() > 0:
            wire_inner.Reverse()
        make_face.Add(wire_inner)

    make_face.Build()
    return make_face.Shape()


from src.base_occ import dispocc, gen_ellipsoid, set_loc, spl_face

from OCC.Extend.DataExchange import write_iges_file, read_iges_file
box = read_iges_file("./box.iges")
pln = load_section_from_iges("./box.iges")

obj = dispocc()
obj.show_axs_pln(scale=50)
obj.display.DisplayShape(box, transparency=0.9)
obj.display.DisplayShape(pln)
obj.ShowOCC()
