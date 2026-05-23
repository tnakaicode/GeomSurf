import numpy as np
import sys
import os
import argparse

basename = os.path.dirname(__file__)

sys.path.append(os.path.join("./"))
from src.base_occ import dispocc

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax2, gp_Ax3
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeTorus, BRepPrimAPI_MakeCylinder
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Cut


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
    gear = BRepPrimAPI_MakeTorus(major_radius, minor_radius).Shape()

    helix_angle = np.deg2rad(helix_angle_deg)
    helix_scale = np.tan(helix_angle)

    for i in range(tooth_count):
        th = 2.0 * np.pi * i / tooth_count
        cth = np.cos(th)
        sth = np.sin(th)

        center = gp_Pnt(major_radius * cth, major_radius * sth, 0.0)
        tangent = gp_Vec(-sth, cth, 0.0)

        helix_vec = gp_Vec(0.0, 0.0, helix_scale)
        cutter_a = _tooth_cutter(center, tangent, helix_vec, groove_radius, groove_length)
        gear = BRepAlgoAPI_Cut(gear, cutter_a).Shape()

        if double_helix:
            cutter_b = _tooth_cutter(center, tangent, gp_Vec(0.0, 0.0, -helix_scale), groove_radius, groove_length)
            gear = BRepAlgoAPI_Cut(gear, cutter_b).Shape()

    return gear

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument("--pxyz", dest="pxyz",
                        default=[0.0, 0.0, 0.0], type=float, nargs=3)
    parser.add_argument("--major", type=float, default=120.0)
    parser.add_argument("--minor", type=float, default=40.0)
    parser.add_argument("--teeth", type=int, default=36)
    parser.add_argument("--groove", type=float, default=9.0)
    parser.add_argument("--helix", type=float, default=28.0,
                        help="helix angle in degree")
    parser.add_argument("--single", action="store_true",
                        help="use single helix grooves (disable double helix)")
    opt = parser.parse_args()
    print(opt)

    obj = dispocc(touch=True)

    gear = make_torus_gear(
        major_radius=opt.major,
        minor_radius=opt.minor,
        tooth_count=opt.teeth,
        groove_radius=opt.groove,
        groove_length=2.5 * (opt.major + opt.minor),
        helix_angle_deg=opt.helix,
        double_helix=not opt.single,
    )

    obj.show_axs_pln(gp_Ax3(), scale=opt.major * 0.4)
    obj.display.DisplayShape(gear, color="GRAY", transparency=0.05)
    obj.ShowOCC()
