import numpy as np
import sys
import os

basename = os.path.dirname(__file__)

sys.path.append(os.path.join("./"))
from src.base_occ import dispocc

import logging

logging.getLogger("matplotlib").setLevel(logging.ERROR)

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeTorus, BRepPrimAPI_MakePrism
from OCC.Core.BRepBuilderAPI import (
    BRepBuilderAPI_MakeEdge,
    BRepBuilderAPI_MakeWire,
    BRepBuilderAPI_MakeFace,
)
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_MakePipeShell
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
from OCC.Core.TColgp import TColgp_HArray1OfPnt
from OCC.Core.GeomAPI import GeomAPI_Interpolate


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


def torus_signed_distance(major_radius, minor_radius, point):
    """Signed distance to the torus surface: negative is inside."""
    rho = np.hypot(point.X(), point.Y())
    return (rho - major_radius) ** 2 + point.Z() ** 2 - minor_radius**2


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
    for i in range(steps):
        t = 2.0 * np.pi * i / steps
        u = u0 + tilt * np.sin(t)
        v = t + v_phase
        pts.append(torus_surface_point(major_radius, minor_radius, u, v))

    pts_array = TColgp_HArray1OfPnt(1, len(pts))
    for idx, p in enumerate(pts, start=1):
        pts_array.SetValue(idx, p)

    interp = GeomAPI_Interpolate(pts_array, True, 1e-6)
    interp.Perform()
    curve = interp.Curve()
    edge = BRepBuilderAPI_MakeEdge(curve).Edge()
    wire = BRepBuilderAPI_MakeWire(edge).Wire()
    return wire, pts


def make_section_face(
    section_width,
    section_height,
    p0,
    major_dir,
    tangent_vec,
    normal_vec,
):
    """Create a parallelogram section whose plane is perpendicular to the wire."""
    tangent_dir = gp_Vec(tangent_vec.X(), tangent_vec.Y(), tangent_vec.Z())
    if tangent_dir.Magnitude() < 1e-9:
        tangent_dir = gp_Vec(1, 0, 0)
    tangent_dir.Normalize()

    normal_dir = gp_Vec(normal_vec.X(), normal_vec.Y(), normal_vec.Z())
    if normal_dir.Magnitude() < 1e-9:
        normal_dir = gp_Vec(0, 0, 1)
    normal_dir.Normalize()

    x_dir = normal_dir.Crossed(tangent_dir)
    if x_dir.Magnitude() < 1e-9:
        x_dir = gp_Vec(1, 0, 0)
    x_dir.Normalize()

    major_dir = gp_Vec(major_dir.X(), major_dir.Y(), major_dir.Z())
    if major_dir.Magnitude() < 1e-9:
        major_dir = gp_Vec(p0.X(), p0.Y(), 0.0)
    if x_dir.Dot(major_dir) > 0:
        x_dir.Reverse()

    half_width = section_width
    half_height = section_height

    p00 = p0.Translated(x_dir.Scaled(-half_width)).Translated(normal_dir.Scaled(-half_height))
    p10 = p0.Translated(x_dir.Scaled(half_width)).Translated(normal_dir.Scaled(-half_height))
    p11 = p0.Translated(x_dir.Scaled(half_width)).Translated(normal_dir.Scaled(half_height))
    p01 = p0.Translated(x_dir.Scaled(-half_width)).Translated(normal_dir.Scaled(half_height))

    wire_builder = BRepBuilderAPI_MakeWire()
    wire_builder.Add(BRepBuilderAPI_MakeEdge(p00, p01).Edge())
    wire_builder.Add(BRepBuilderAPI_MakeEdge(p01, p11).Edge())
    wire_builder.Add(BRepBuilderAPI_MakeEdge(p11, p10).Edge())
    wire_builder.Add(BRepBuilderAPI_MakeEdge(p10, p00).Edge())
    section_wire = wire_builder.Wire()

    return BRepBuilderAPI_MakeFace(section_wire).Face(), section_wire


def make_torus_wires(
    major_radius=120.0,
    minor_radius=40.0,
    wires=5,
    steps=240,
    section_width=16.0,
    section_height=8.0,
    tilt=0.08,
):
    """Create multiple torus-surface wires and sweep section faces along them."""
    print("[DEBUG] make_torus_wires start")
    print(f"[DEBUG]   major_radius: {major_radius}")
    print(f"[DEBUG]   minor_radius: {minor_radius}")
    print(f"[DEBUG]   wires       : {wires}")
    print(f"[DEBUG]   steps       : {steps}")
    print(f"[DEBUG]   section_width: {section_width}")
    print(f"[DEBUG]   section_height: {section_height}")
    print(f"[DEBUG]   tilt        : {tilt}")

    torus = BRepPrimAPI_MakeTorus(major_radius, minor_radius).Shape()
    sweep_list = []
    wire_list = []
    section_list = []
    common_shapes = []
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
        wire_list.append(wire)

        major_dir = gp_Vec(pts[0].X(), pts[0].Y(), 0.0)
        if major_dir.Magnitude() < 1e-9:
            major_dir = gp_Vec(np.cos(u0), np.sin(u0), 0.0)
        major_dir.Normalize()

        normal_vec = torus_surface_normal(u0, v_phase)
        tangent_vec = gp_Vec(pts[-1], pts[1])
        tangent_dir = gp_Vec(tangent_vec.X(), tangent_vec.Y(), tangent_vec.Z())
        if tangent_dir.Magnitude() < 1e-9:
            tangent_dir = gp_Vec(1, 0, 0)
        tangent_dir.Normalize()

        inward1 = normal_vec.Crossed(tangent_dir)
        inward2 = tangent_dir.Crossed(normal_vec)
        if inward1.Magnitude() < 1e-9:
            x_dir = inward2
        elif inward2.Magnitude() < 1e-9:
            x_dir = inward1
        else:
            inward1.Normalize()
            inward2.Normalize()
            test_pt = pts[0].Translated(inward1.Scaled(section_width * 0.2))
            x_dir = inward1 if torus_signed_distance(major_radius, minor_radius, test_pt) < 0 else inward2

        center = pts[0].Translated(normal_vec.Reversed().Scaled(2.0))
        sd = torus_signed_distance(major_radius, minor_radius, center)
        print(
            f"[DEBUG]   section center signed distance: {sd:.6f} ({'inside' if sd < 0 else 'outside' if sd > 0 else 'on surface'})"
        )
        section, section_wire = make_section_face(
            section_width,
            section_height,
            pts[0],
            major_dir,
            tangent_vec,
            normal_vec,
        )
        section_list.append(section)

        pipe_shell = BRepOffsetAPI_MakePipeShell(wire)
        pipe_shell.Add(section_wire, True)
        pipe_shell.Build()
        if not pipe_shell.IsDone():
            print(f"[DEBUG]   wire {i + 1}/{wires} MakePipeShell failed, fallback to prism")
            cutter = BRepPrimAPI_MakePrism(section, tangent_vec.Scaled(1000)).Shape()
        else:
            pipe_shell.MakeSolid()
            cutter = pipe_shell.Shape()
        sweep_list.append(cutter)
        print(f"[DEBUG]   wire {i + 1}/{wires} cutter built (pipe shell sweep, type={cutter.ShapeType()})")

    result = torus
    for idx, cutter in enumerate(sweep_list, start=1):
        print(
            f"[DEBUG] cut_{idx}: result type={type(result)}, cutter type={type(cutter)}, cutter ShapeType={cutter.ShapeType()}"
        )
        cut_op = BRepAlgoAPI_Cut(result, cutter)
        cut_op.Build()
        if not cut_op.IsDone():
            print(f"[DEBUG] cut_{idx}: operation failed")
            continue
        new_result = cut_op.Shape()
        if new_result is None or new_result.IsNull():
            print(f"[DEBUG] cut_{idx}: result is null")
            continue
        result = new_result
        print(f"[DEBUG] cut_{idx}: removed cutter volume")

    print("[DEBUG] make_torus_wires end")
    return torus, result, wire_list, section_list, sweep_list, common_shapes


def debug_shape_info(shape, name="shape"):
    """Print minimal diagnostics to confirm geometry creation state."""
    print(f"[DEBUG] {name}.IsNull   : {shape.IsNull()}")
    print(f"[DEBUG] {name}.ShapeType: {shape.ShapeType()}")


if __name__ == "__main__":
    major_radius = 120.0
    minor_radius = 40.0
    wire_count = 5
    wire_steps = 360
    section_width = 16.0
    section_height = 8.0
    tilt = 0.2

    print("[DEBUG] Build parameters")
    print(f"[DEBUG]   major_radius : {major_radius}")
    print(f"[DEBUG]   minor_radius : {minor_radius}")
    print(f"[DEBUG]   wire_count   : {wire_count}")
    print(f"[DEBUG]   wire_steps   : {wire_steps}")
    print(f"[DEBUG]   section_width: {section_width}")
    print(f"[DEBUG]   section_height: {section_height}")
    print(f"[DEBUG]   tilt         : {tilt}")

    obj = dispocc(touch=True)

    torus, final_shape, wire_list, section_list, sweep_list, common_shapes = make_torus_wires(
        major_radius=major_radius,
        minor_radius=minor_radius,
        wires=wire_count,
        steps=wire_steps,
        section_width=section_width,
        section_height=section_height,
        tilt=tilt
    )

    #if final_shape is not None and not final_shape.IsNull():
    #    obj.display.DisplayShape(final_shape, transparency=0.5)
    #    print("[DEBUG] final result displayed")

    debug_shape_info(final_shape, "final_shape")
    for idx, wire in enumerate(wire_list, start=1):
        obj.display.DisplayShape(wire, color="RED")
    for idx, section in enumerate(section_list, start=1):
        obj.display.DisplayShape(section, transparency=0.7, color="BLUE")
    #for idx, sweep in enumerate(sweep_list, start=1):
    #    obj.display.DisplayShape(sweep, transparency=0.7, color="GREEN")
    #for idx, common_shape in enumerate(common_shapes, start=1):
    #    obj.display.DisplayShape(common_shape, transparency=0.5, color="BLUE1")
    obj.display.DisplayShape(torus)
    obj.show_axs_pln(gp_Ax3(), scale=major_radius * 0.4)
    obj.save_view("debug")
    obj.ShowOCC()
