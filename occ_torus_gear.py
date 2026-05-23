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
from OCC.Core.Geom import Geom_Circle


def torus_surface_point(major_radius, minor_radius, u, v):
    """Return a 3D point on the torus surface for parameters u and v."""
    return gp_Pnt(
        (major_radius + minor_radius * np.cos(v)) * np.cos(u),
        (major_radius + minor_radius * np.cos(v)) * np.sin(u),
        minor_radius * np.sin(v),
    )


def torus_surface_normal(u, v):
    """Return the normalized torus surface normal for parameters u and v."""
    return gp_Vec(np.cos(v) * np.cos(u), np.cos(v) * np.sin(u), np.sin(v))


def torus_surface_wire(
    major_radius,
    minor_radius,
    steps=144,
    u0=0.0,
    tilt=0.08,
    v_phase=0.0,
):
    """Build a closed wire around the torus tube with a slight tilt."""
    pts = []
    for i in range(steps + 1):
        t = 2.0 * np.pi * i / steps
        u = u0 + tilt * np.sin(t)
        v = t + v_phase
        pts.append(torus_surface_point(major_radius, minor_radius, u, v))

    wire_builder = BRepBuilderAPI_MakeWire()
    for i in range(steps):
        edge = BRepBuilderAPI_MakeEdge(pts[i], pts[i + 1]).Edge()
        wire_builder.Add(edge)
    return wire_builder.Wire(), pts


def make_section_face(section_radius, p0, tangent_vec, normal_vec, bite_depth=1.5):
    """Create a planar circle face perpendicular to the path tangent and bitten into the torus."""
    center = p0.Translated(normal_vec.Reversed().Scaled(bite_depth))
    tangent_dir = gp_Dir(tangent_vec)
    alt_dir = gp_Dir(1, 0, 0) if abs(tangent_dir.Dot(gp_Dir(0, 0, 1))) > 0.9 else gp_Dir(0, 0, 1)
    ax2 = gp_Ax2(center, tangent_dir, alt_dir)
    circle = Geom_Circle(ax2, section_radius)
    edge = BRepBuilderAPI_MakeEdge(circle).Edge()
    wire = BRepBuilderAPI_MakeWire(edge).Wire()
    return BRepBuilderAPI_MakeFace(wire).Face()


def make_torus_wires(
    major_radius=120.0,
    minor_radius=40.0,
    wires=5,
    steps=240,
    section_radius=8.0,
    tilt=0.08,
):
    """Create multiple torus-surface wires and sweep section faces along them."""
    print("[DEBUG] make_torus_wires start")
    print(f"[DEBUG]   major_radius: {major_radius}")
    print(f"[DEBUG]   minor_radius: {minor_radius}")
    print(f"[DEBUG]   wires       : {wires}")
    print(f"[DEBUG]   steps       : {steps}")
    print(f"[DEBUG]   section_radius: {section_radius}")
    print(f"[DEBUG]   tilt        : {tilt}")

    torus = BRepPrimAPI_MakeTorus(major_radius, minor_radius).Shape()
    wires_list = []
    sweep_list = []
    for i in range(wires):
        u0 = 2.0 * np.pi * i / wires
        v_phase = 0.0
        wire, pts = torus_surface_wire(
            major_radius,
            minor_radius,
            steps=steps,
            u0=u0,
            tilt=tilt,
            v_phase=v_phase,
        )
        wires_list.append(wire)

        tangent_vec = gp_Vec(pts[0], pts[1])
        normal_vec = torus_surface_normal(u0, v_phase)
        section = make_section_face(section_radius, pts[0], tangent_vec, normal_vec, bite_depth=2.0)
        pipe = BRepOffsetAPI_MakePipe(wire, section)
        pipe.Build()
        sweep = pipe.Shape()
        sweep_list.append(sweep)
        print(f"[DEBUG]   wire {i + 1}/{wires} swept")

    print("[DEBUG] make_torus_wires end")
    return torus, wires_list, sweep_list


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

    torus, wires, sweeps = make_torus_wires(
        major_radius=major,
        minor_radius=minor,
        wires=5,
        steps=360,
        section_radius=8.0,
        tilt=0.08
    )

    debug_shape_info(torus, "torus")
    obj.display.DisplayShape(torus, color="BLUE", transparency=0.8)
    for idx, wire in enumerate(wires, start=1):
        debug_shape_info(wire, f"wire_{idx}")
        obj.display.DisplayShape(wire, color="RED" if idx % 2 == 0 else "BLUE", transparency=0.0)

    for idx, sweep in enumerate(sweeps, start=1):
        debug_shape_info(sweep, f"sweep_{idx}")
        obj.display.DisplayShape(sweep, transparency=0.5)

    obj.show_axs_pln(gp_Ax3(), scale=major * 0.4)
    obj.ShowOCC()
