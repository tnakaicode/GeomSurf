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
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Common
from OCC.Core.BRepCheck import BRepCheck_Analyzer
from OCC.Core.TColgp import TColgp_HArray1OfPnt
from OCC.Core.GeomAPI import GeomAPI_Interpolate
from OCC.Extend.TopologyUtils import TopologyExplorer


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
    return (rho - major_radius) ** 2 + point.Z() ** 2 - minor_radius ** 2


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


def make_section_face(section_width, section_height, p0, tangent_vec, normal_vec, bite_depth=1.5):
    """Create a planar rectangular face perpendicular to the path tangent."""
    center = p0.Translated(normal_vec.Reversed().Scaled(bite_depth))
    tangent_dir = gp_Dir(tangent_vec)
    # choose a stable in-plane axis
    alt_dir = gp_Dir(1, 0, 0) if abs(tangent_dir.Dot(gp_Dir(0, 0, 1))) > 0.9 else gp_Dir(0, 0, 1)
    x_dir = gp_Vec(alt_dir)
    x_dir.Normalize()
    # make x_dir perpendicular to tangent_dir
    tangent_vec_norm = gp_Vec(tangent_dir)
    x_dir = x_dir.Subtracted(tangent_vec_norm.Scaled(x_dir.Dot(tangent_vec_norm)))
    x_dir.Normalize()
    x_dir_dir = gp_Dir(x_dir)
    y_dir = gp_Vec(tangent_dir.Crossed(x_dir_dir))

    hx = section_width
    hy = section_height
    p00 = center.Translated(x_dir.Scaled(-hx)).Translated(y_dir.Scaled(-hy))
    p01 = center.Translated(x_dir.Scaled(hx)).Translated(y_dir.Scaled(-hy))
    p11 = center.Translated(x_dir.Scaled(hx)).Translated(y_dir.Scaled(hy))
    p10 = center.Translated(x_dir.Scaled(-hx)).Translated(y_dir.Scaled(hy))

    wire_builder = BRepBuilderAPI_MakeWire()
    wire_builder.Add(BRepBuilderAPI_MakeEdge(p00, p01).Edge())
    wire_builder.Add(BRepBuilderAPI_MakeEdge(p01, p11).Edge())
    wire_builder.Add(BRepBuilderAPI_MakeEdge(p11, p10).Edge())
    wire_builder.Add(BRepBuilderAPI_MakeEdge(p10, p00).Edge())
    return BRepBuilderAPI_MakeFace(wire_builder.Wire()).Face()


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
        center = pts[0].Translated(normal_vec.Reversed().Scaled(2.0))
        sd = torus_signed_distance(major_radius, minor_radius, center)
        print(f"[DEBUG]   section center signed distance: {sd:.6f} ({'inside' if sd < 0 else 'outside' if sd > 0 else 'on surface'})")
        section = make_section_face(section_width, section_height, pts[0], tangent_vec, normal_vec, bite_depth=2.0)
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
        section_width=16.0,
        section_height=8.0,
        tilt=0.2,
    )

    debug_shape_info(torus, "torus")
    obj.display.DisplayShape(torus, color="BLUE1", transparency=0.5)
    for idx, wire in enumerate(wires, start=1):
        debug_shape_info(wire, f"wire_{idx}")
        obj.display.DisplayShape(wire, color="RED", transparency=0.0)

    common_shapes = []
    for idx, sweep in enumerate(sweeps, start=1):
        common = BRepAlgoAPI_Common(torus, sweep)
        if not common.IsDone():
            print(f"[DEBUG] common_{idx}: operation failed")
            continue
        common_shape = common.Shape()
        if common_shape.IsNull():
            print(f"[DEBUG] common_{idx}: no intersection")
            continue
        print(f"[DEBUG] common_{idx}: intersection found (type={common_shape.ShapeType()})")
        common_shapes.append(common_shape)

    if common_shapes:
        print(f"[DEBUG] Displaying {len(common_shapes)} common shape(s)")
        for idx, common_shape in enumerate(common_shapes, start=1):
            debug_shape_info(common_shape, f"common_{idx}")
            obj.display.DisplayShape(common_shape, color="GREEN", transparency=0.5)
    else:
        print("[DEBUG] No common intersection found for any sweep.")

    obj.show_axs_pln(gp_Ax3(), scale=major * 0.4)
    obj.ShowOCC()
