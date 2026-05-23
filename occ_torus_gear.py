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
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire


def torus_surface_point(major_radius, minor_radius, u, v):
    """Return a 3D point on the torus surface for parameters u and v."""
    return gp_Pnt(
        (major_radius + minor_radius * np.cos(v)) * np.cos(u),
        (major_radius + minor_radius * np.cos(v)) * np.sin(u),
        minor_radius * np.sin(v),
    )


def torus_surface_wire(
    major_radius,
    minor_radius,
    steps=144,
    u0=0.0,
    skew=0.08,
    v_phase=0.0,
):
    """Build a closed wire around the torus tube with a slight major-angle slant."""
    pts = []
    for i in range(steps + 1):
        t = 2.0 * np.pi * i / steps
        u = u0 + skew * np.sin(t)
        v = t + v_phase
        pts.append(torus_surface_point(major_radius, minor_radius, u, v))

    wire_builder = BRepBuilderAPI_MakeWire()
    for i in range(steps):
        edge = BRepBuilderAPI_MakeEdge(pts[i], pts[i + 1]).Edge()
        wire_builder.Add(edge)
    return wire_builder.Wire()


def make_torus_wires(
    major_radius=120.0,
    minor_radius=40.0,
    wires=5,
    steps=240,
):
    """Create multiple torus-surface wires around the minor radius."""
    print("[DEBUG] make_torus_wires start")
    print(f"[DEBUG]   major_radius: {major_radius}")
    print(f"[DEBUG]   minor_radius: {minor_radius}")
    print(f"[DEBUG]   wires       : {wires}")
    print(f"[DEBUG]   steps       : {steps}")

    torus = BRepPrimAPI_MakeTorus(major_radius, minor_radius).Shape()
    wires_list = []
    for i in range(wires):
        u0 = 2.0 * np.pi * i / wires
        v_phase = 0.0
        wire = torus_surface_wire(
            major_radius,
            minor_radius,
            steps=steps,
            u0=u0,
            v_phase=v_phase,
        )
        wires_list.append(wire)
        print(f"[DEBUG]   wire {i + 1}/{wires} u0={u0:.3f} created")

    print("[DEBUG] make_torus_wires end")
    return torus, wires_list


def debug_shape_info(shape, name="shape"):
    """Print minimal diagnostics to confirm geometry creation state."""
    print(f"[DEBUG] {name}.IsNull   : {shape.IsNull()}")
    print(f"[DEBUG] {name}.ShapeType: {shape.ShapeType()}")


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

    torus, wires = make_torus_wires(
        major_radius=major,
        minor_radius=minor,
        wires=5,
        steps=360,
    )

    debug_shape_info(torus, "torus")
    obj.display.DisplayShape(torus, color="BLUE", transparency=0.8)
    for idx, wire in enumerate(wires, start=1):
        debug_shape_info(wire, f"wire_{idx}")
        obj.display.DisplayShape(wire, color="RED" if idx % 2 == 0 else "BLUE", transparency=0.0)

    obj.show_axs_pln(gp_Ax3(), scale=major * 0.4)
    obj.ShowOCC()
