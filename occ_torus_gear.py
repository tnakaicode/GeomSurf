import numpy as np
import sys
import os

basename = os.path.dirname(__file__)

sys.path.append(os.path.join("./"))
from src.base_occ import dispocc

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax2, gp_Ax3
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeTorus, BRepPrimAPI_MakeCylinder
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Cut
from OCC.Core.BRepCheck import BRepCheck_Analyzer
from OCC.Extend.TopologyUtils import TopologyExplorer


def _tooth_cutter(center, tangent, helix_dir, groove_radius, groove_length):
    """Build one cylindrical cutter used to carve a tooth groove."""
    axis_vec = tangent + helix_dir
    axis_vec.Normalize()
    offset = gp_Vec(axis_vec.X(), axis_vec.Y(), axis_vec.Z())
    offset.Multiply(0.5 * groove_length)
    start = gp_Pnt(center.X() - offset.X(), center.Y() - offset.Y(), center.Z() - offset.Z())

    # Use radial direction as X-direction for a stable local frame.
    radial = gp_Dir(center.X(), center.Y(), 0.0)
    cutter_axs = gp_Ax2(start, gp_Dir(axis_vec.X(), axis_vec.Y(), axis_vec.Z()), radial)
    return BRepPrimAPI_MakeCylinder(cutter_axs, groove_radius, groove_length).Shape()


def make_torus_gear(
    major_radius=120.0,
    minor_radius=40.0,
    tooth_count=36,
    groove_radius=9.0,
    groove_length=320.0,
    helix_angle_deg=28.0,
    double_helix=True,
):
    """Create a torus gear by cutting repeated helical grooves on a torus."""
    print("[DEBUG] make_torus_gear start")
    print(f"[DEBUG]   major_radius   : {major_radius}")
    print(f"[DEBUG]   minor_radius   : {minor_radius}")
    print(f"[DEBUG]   tooth_count    : {tooth_count}")
    print(f"[DEBUG]   groove_radius  : {groove_radius}")
    print(f"[DEBUG]   groove_length  : {groove_length}")
    print(f"[DEBUG]   helix_angle_deg: {helix_angle_deg}")
    print(f"[DEBUG]   double_helix   : {double_helix}")

    gear = BRepPrimAPI_MakeTorus(major_radius, minor_radius).Shape()
    print(f"[DEBUG]   base_torus_is_null: {gear.IsNull()}")

    helix_angle = np.deg2rad(helix_angle_deg)
    helix_scale = np.tan(helix_angle)
    print(f"[DEBUG]   helix_angle_rad: {helix_angle}")
    print(f"[DEBUG]   helix_scale    : {helix_scale}")

    for i in range(tooth_count):
        th = 2.0 * np.pi * i / tooth_count
        cth = np.cos(th)
        sth = np.sin(th)

        center = gp_Pnt(major_radius * cth, major_radius * sth, 0.0)
        tangent = gp_Vec(-sth, cth, 0.0)
        print(
            f"[DEBUG] tooth {i + 1}/{tooth_count} | "
            f"th={th:.6f} rad | center=({center.X():.3f}, {center.Y():.3f}, {center.Z():.3f})"
        )

        helix_vec = gp_Vec(0.0, 0.0, helix_scale)
        cutter_a = _tooth_cutter(center, tangent, helix_vec, groove_radius, groove_length)
        print(f"[DEBUG]   cutter_a_is_null: {cutter_a.IsNull()}")
        gear = BRepAlgoAPI_Cut(gear, cutter_a).Shape()
        cut_a_null = gear.IsNull()
        cut_a_valid = False if cut_a_null else BRepCheck_Analyzer(gear).IsValid()
        print(f"[DEBUG]   after_cut_a -> is_null={cut_a_null}, is_valid={cut_a_valid}")

        if double_helix:
            cutter_b = _tooth_cutter(center, tangent, gp_Vec(0.0, 0.0, -helix_scale), groove_radius, groove_length)
            print(f"[DEBUG]   cutter_b_is_null: {cutter_b.IsNull()}")
            gear = BRepAlgoAPI_Cut(gear, cutter_b).Shape()
            cut_b_null = gear.IsNull()
            cut_b_valid = False if cut_b_null else BRepCheck_Analyzer(gear).IsValid()
            print(f"[DEBUG]   after_cut_b -> is_null={cut_b_null}, is_valid={cut_b_valid}")

    print("[DEBUG] make_torus_gear end")
    return gear


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
    teeth = 36
    groove = 9.0
    helix = 28.0
    single = False

    print("[DEBUG] Build parameters")
    print(f"[DEBUG]   major     : {major}")
    print(f"[DEBUG]   minor     : {minor}")
    print(f"[DEBUG]   teeth     : {teeth}")
    print(f"[DEBUG]   groove    : {groove}")
    print(f"[DEBUG]   helix_deg : {helix}")
    print(f"[DEBUG]   single    : {single}")

    obj = dispocc(touch=True)

    gear = make_torus_gear(
        major_radius=major,
        minor_radius=minor,
        tooth_count=teeth,
        groove_radius=groove,
        groove_length=2.5 * (major + minor),
        helix_angle_deg=helix,
        double_helix=not single,
    )
    debug_shape_info(gear)

    obj.show_axs_pln(gp_Ax3(), scale=major * 0.4)
    obj.display.DisplayShape(gear)
    obj.ShowOCC()
