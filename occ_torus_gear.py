import numpy as np
import sys
import os

basename = os.path.dirname(__file__)

sys.path.append(os.path.join("./"))
from src.base_occ import dispocc

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeTorus
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire, BRepBuilderAPI_MakeFace
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_MakePipe
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Cut
from OCC.Core.BRepCheck import BRepCheck_Analyzer
from OCC.Core.Geom import Geom_Circle
from OCC.Extend.TopologyUtils import TopologyExplorer


def torus_surface_point(major_radius, minor_radius, u, v):
    """Return a 3D point on the torus surface for parameters u and v."""
    return gp_Pnt(
        (major_radius + minor_radius * np.cos(v)) * np.cos(u),
        (major_radius + minor_radius * np.cos(v)) * np.sin(u),
        minor_radius * np.sin(v),
    )


def torus_surface_wire(major_radius, minor_radius, steps=144, v_amplitude=0.12, tooth_count=24):
    """Build a closed wire along the torus surface for one full major loop.

    The path travels once around the major circle (u from 0 to 2π), while the
    minor parameter v oscillates according to the tooth_count to create multiple
    grooves in the small-radius direction.
    """
    pts = []
    for i in range(steps + 1):
        u = 2.0 * np.pi * i / steps
        v = v_amplitude * np.sin(tooth_count * u)
        pts.append(torus_surface_point(major_radius, minor_radius, u, v))

    wire_builder = BRepBuilderAPI_MakeWire()
    for i in range(steps):
        edge = BRepBuilderAPI_MakeEdge(pts[i], pts[i + 1]).Edge()
        wire_builder.Add(edge)
    return wire_builder.Wire()


def section_face(radius=9.0):
    """Create a planar circle face to sweep along the path."""
    circle = Geom_Circle(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1), gp_Dir(1, 0, 0)), radius)
    edge = BRepBuilderAPI_MakeEdge(circle).Edge()
    wire = BRepBuilderAPI_MakeWire(edge).Wire()
    face = BRepBuilderAPI_MakeFace(wire).Face()
    return face


def make_torus_gear(
    major_radius=120.0,
    minor_radius=40.0,
    path_steps=360,
    section_radius=25.0,
    v_amplitude=0.6,
    tooth_count=12,
):
    """Build a torus gear by sweeping a section along a torus surface wire and cutting the torus."""
    print("[DEBUG] make_torus_gear start")
    print(f"[DEBUG]   major_radius : {major_radius}")
    print(f"[DEBUG]   minor_radius : {minor_radius}")
    print(f"[DEBUG]   path_steps   : {path_steps}")
    print(f"[DEBUG]   section_radius: {section_radius}")
    print(f"[DEBUG]   v_amplitude  : {v_amplitude}")
    print(f"[DEBUG]   tooth_count  : {tooth_count}")

    torus = BRepPrimAPI_MakeTorus(major_radius, minor_radius).Shape()
    print(f"[DEBUG]   base_torus_is_null: {torus.IsNull()}")

    path = torus_surface_wire(major_radius, minor_radius, path_steps, v_amplitude, tooth_count)
    section = section_face(section_radius)
    pipe = BRepOffsetAPI_MakePipe(path, section)
    pipe.Build()
    swept = pipe.Shape()
    swept_null = swept is None or swept.IsNull()
    print(f"[DEBUG]   swept_null: {swept_null}")
    if swept_null:
        print("[DEBUG]   pipe sweep failed")
        return None

    cut_builder = BRepAlgoAPI_Cut(torus, swept)
    if not cut_builder.IsDone():
        print("[DEBUG]   cut operation not done")
        return None
    gear = cut_builder.Shape()
    gear_null = gear is None or gear.IsNull()
    print(f"[DEBUG]   gear_is_null: {gear_null}")
    if not gear_null:
        gear_valid = BRepCheck_Analyzer(gear).IsValid()
    else:
        gear_valid = False
    print(f"[DEBUG]   gear_is_valid: {gear_valid}")
    print("[DEBUG] make_torus_gear end")
    return gear, path, swept


def debug_shape_info(shape):
    """Print minimal diagnostics to confirm geometry creation state."""
    is_null = shape.IsNull()
    is_valid = False if is_null else BRepCheck_Analyzer(shape).IsValid()

    face_count = 0
    edge_count = 0
    if not is_null:
        topo = TopologyExplorer(shape)
        face_count = len(list(topo.faces()))
        edge_count = len(list(topo.edges()))

    print("[DEBUG] Torus gear build status")
    print(f"[DEBUG]   is_null   : {is_null}")
    print(f"[DEBUG]   is_valid  : {is_valid}")
    print(f"[DEBUG]   faces     : {face_count}")
    print(f"[DEBUG]   edges     : {edge_count}")

if __name__ == '__main__':
    major = 120.0
    minor = 40.0
    path_steps = 360
    section_radius = 14.0
    v_amplitude = 0.35
    turns = 12

    print("[DEBUG] Build parameters")
    print(f"[DEBUG]   major        : {major}")
    print(f"[DEBUG]   minor        : {minor}")
    print(f"[DEBUG]   path_steps   : {path_steps}")
    print(f"[DEBUG]   section_radius: {section_radius}")
    print(f"[DEBUG]   v_amplitude  : {v_amplitude}")
    print(f"[DEBUG]   tooth_count  : {turns}")

    obj = dispocc(touch=True)

    gear, path, swept = make_torus_gear(
        major_radius=major,
        minor_radius=minor,
        path_steps=path_steps,
        section_radius=section_radius,
        v_amplitude=v_amplitude,
        tooth_count=turns,
    )
    debug_shape_info(gear)
    obj.display.DisplayShape(path, color="BLUE", transparency=0.7)
    obj.display.DisplayShape(swept, color="RED", transparency=0.6)
    obj.show_axs_pln(gp_Ax3(), scale=major * 0.4)
    obj.display.DisplayShape(gear, transparency=0.7)
    obj.ShowOCC()
