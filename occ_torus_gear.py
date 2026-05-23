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
    BRepBuilderAPI_MakeSolid,
    BRepBuilderAPI_Sewing,
)
from OCC.Core.TopoDS import topods
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_MakePipe, BRepOffsetAPI_MakePipeShell, BRepOffsetAPI_ThruSections
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Cut, BRepAlgoAPI_Common, BRepAlgoAPI_Fuse
from OCC.Core.ShapeFix import ShapeFix_Solid, ShapeFix_Shape
from OCC.Core.BRepGProp import brepgprop
from OCC.Core.GProp import GProp_GProps
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
    return wire, pts, curve


def shape_fix_shape(shape):
    if shape is None or shape.IsNull():
        return shape
    if shape.ShapeType() == 2:
        fixer = ShapeFix_Solid()
        fixer.Init(shape)
        fixer.Perform()
        fixed = fixer.Shape()
    else:
        fixer = ShapeFix_Shape()
        fixer.Init(shape)
        fixer.Perform()
        fixed = fixer.Shape()
    if fixed is None or fixed.IsNull():
        return shape
    return fixed


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

    y_dir = normal_dir

    half_width = section_width * 0.5
    half_height = section_height * 0.5
    # Center section on the torus surface: outer edge on surface, inner edge penetrates by section_height
    bite_depth = half_height

    center = p0.Translated(normal_dir.Reversed().Scaled(bite_depth))
    p00 = center.Translated(x_dir.Scaled(-half_width)).Translated(y_dir.Scaled(-half_height))
    p10 = center.Translated(x_dir.Scaled(half_width)).Translated(y_dir.Scaled(-half_height))
    p11 = center.Translated(x_dir.Scaled(half_width)).Translated(y_dir.Scaled(half_height))
    p01 = center.Translated(x_dir.Scaled(-half_width)).Translated(y_dir.Scaled(half_height))

    wire_builder = BRepBuilderAPI_MakeWire()
    wire_builder.Add(BRepBuilderAPI_MakeEdge(p00, p01).Edge())
    wire_builder.Add(BRepBuilderAPI_MakeEdge(p01, p11).Edge())
    wire_builder.Add(BRepBuilderAPI_MakeEdge(p11, p10).Edge())
    wire_builder.Add(BRepBuilderAPI_MakeEdge(p10, p00).Edge())
    section_wire = wire_builder.Wire()

    return BRepBuilderAPI_MakeFace(section_wire).Face(), section_wire


def count_solids(shape):
    return sum(1 for _ in TopologyExplorer(shape).solids())


def shape_volume(shape):
    props = GProp_GProps()
    try:
        brepgprop.VolumeProperties(shape, props)
        return props.Mass()
    except Exception:
        return 0.0


def choose_largest_solid(shape):
    solids = list(TopologyExplorer(shape).solids())
    if not solids:
        return None
    if len(solids) == 1:
        return solids[0]
    best_solid = None
    best_abs_volume = -1.0
    for solid in solids:
        props = GProp_GProps()
        brepgprop.VolumeProperties(solid, props)
        volume = props.Mass()
        abs_volume = abs(volume)
        print(f"[DEBUG]   candidate solid volume={volume:.6f}, abs={abs_volume:.6f}")
        if abs_volume > best_abs_volume:
            best_abs_volume = abs_volume
            best_solid = solid
    print(f"[DEBUG]   chosen largest solid abs volume={best_abs_volume:.6f}")
    return best_solid


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
    cutters_for_cut = []
    prism_cutters = []
    wire_list = []
    section_list = []
    common_shapes = []
    common_volumes = []

    for i in range(wires):
        u0 = 2.0 * np.pi * i / wires
        v_phase = 0.0
        wire, pts, curve = torus_surface_wire(
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
        tangent_vec = gp_Vec(1, 0, 0)
        try:
            p0 = gp_Pnt()
            d1 = gp_Vec()
            curve.D1(curve.FirstParameter(), p0, d1)
            tangent_vec = gp_Vec(d1.X(), d1.Y(), d1.Z())
        except Exception:
            tangent_vec = gp_Vec(pts[0], pts[1])
            if tangent_vec.Magnitude() < 1e-9:
                tangent_vec = gp_Vec(pts[-1], pts[0])
        tangent_dir = gp_Vec(tangent_vec.X(), tangent_vec.Y(), tangent_vec.Z())
        if tangent_dir.Magnitude() < 1e-9:
            tangent_dir = gp_Vec(1, 0, 0)
        tangent_dir.Normalize()
        if tangent_dir.Magnitude() < 1e-9:
            tangent_dir = gp_Vec(1, 0, 0)
        tangent_dir.Normalize()

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

        # Build groove cutter as a ring of discrete cells (ThruSections approach).
        # Section extends outside torus surface to avoid tangency issues in boolean ops.
        _n_cells = 16
        _half_w = section_width * 0.5
        _inner_ext = section_height          # depth inside torus (groove depth)
        _outer_ext = section_height * 0.5   # extension OUTSIDE torus surface
        _half_h_n = (_inner_ext + _outer_ext) * 0.5  # half total height along normal
        # center offset from surface: negative = inside torus
        _center_from_surf = (_outer_ext - _inner_ext) * 0.5  # = -section_height/4

        def _make_rect_wire(pk, nk, tk, u_k):
            # Azimuthal direction: always perpendicular to torus normal, no singularity
            az = gp_Vec(-np.sin(u_k), np.cos(u_k), 0.0)
            dot_at = az.Dot(tk)
            xk = gp_Vec(az.X() - dot_at * tk.X(),
                        az.Y() - dot_at * tk.Y(),
                        az.Z() - dot_at * tk.Z())
            if xk.Magnitude() < 1e-9:
                xk = az
            xk.Normalize()
            # Center: slightly inside torus surface (outer edge is outside surface)
            ck = pk.Translated(nk.Scaled(_center_from_surf))
            q00 = ck.Translated(xk.Scaled(-_half_w)).Translated(nk.Scaled(-_half_h_n))  # inner edge
            q10 = ck.Translated(xk.Scaled( _half_w)).Translated(nk.Scaled(-_half_h_n))  # inner edge
            q11 = ck.Translated(xk.Scaled( _half_w)).Translated(nk.Scaled( _half_h_n))  # outer edge
            q01 = ck.Translated(xk.Scaled(-_half_w)).Translated(nk.Scaled( _half_h_n))  # outer edge
            wb = BRepBuilderAPI_MakeWire()
            wb.Add(BRepBuilderAPI_MakeEdge(q00, q01).Edge())
            wb.Add(BRepBuilderAPI_MakeEdge(q01, q11).Edge())
            wb.Add(BRepBuilderAPI_MakeEdge(q11, q10).Edge())
            wb.Add(BRepBuilderAPI_MakeEdge(q10, q00).Edge())
            return wb.Wire() if wb.IsDone() else None

        def _torus_tangent(u_k, v_k, dt=1e-5):
            pn = torus_surface_point(major_radius, minor_radius,
                                     u0 + tilt * np.sin(v_k + dt), v_k + dt)
            pk = torus_surface_point(major_radius, minor_radius, u_k, v_k)
            tv = gp_Vec(pk, pn)
            if tv.Magnitude() < 1e-9:
                tv = gp_Vec(0, 0, 1)
            tv.Normalize()
            return tv

        _cells = []
        for _k in range(_n_cells):
            _t0 = 2.0 * np.pi * _k / _n_cells
            _t1 = 2.0 * np.pi * (_k + 1) / _n_cells
            _u0k = u0 + tilt * np.sin(_t0)
            _u1k = u0 + tilt * np.sin(_t1)
            _p0 = torus_surface_point(major_radius, minor_radius, _u0k, _t0)
            _p1 = torus_surface_point(major_radius, minor_radius, _u1k, _t1)
            _n0 = torus_surface_normal(_u0k, _t0); _n0.Normalize()
            _n1 = torus_surface_normal(_u1k, _t1); _n1.Normalize()
            _tv0 = _torus_tangent(_u0k, _t0)
            _tv1 = _torus_tangent(_u1k, _t1)
            _w0 = _make_rect_wire(_p0, _n0, _tv0, _u0k)
            _w1 = _make_rect_wire(_p1, _n1, _tv1, _u1k)
            if _w0 is None or _w1 is None:
                continue
            _loft = BRepOffsetAPI_ThruSections(True, True)  # solid, ruled
            _loft.AddWire(_w0)
            _loft.AddWire(_w1)
            _loft.Build()
            if _loft.IsDone():
                _cell = shape_fix_shape(_loft.Shape())
                if _cell and not _cell.IsNull():
                    _cells.append(_cell)

        # Store individual cells for sequential cutting (ring fusion is unreliable)
        _cutter_cells = []
        if _cells:
            _cutter_cells = _cells
            # Diagnostic: use first cell as representative for display
            cutter = shape_fix_shape(_cells[0])
            total_cell_vol = sum(shape_volume(c) for c in _cells)
        else:
            print(f"[DEBUG]   wire {i + 1}/{wires} cell-loft failed, fallback to prism")
            cutter = BRepPrimAPI_MakePrism(section, normal_vec.Reversed().Scaled(-section_height * 4)).Shape()
            cutter = shape_fix_shape(cutter)
            _cutter_cells = [cutter] if cutter and not cutter.IsNull() else []
            total_cell_vol = shape_volume(cutter)

        prism_cutter = BRepPrimAPI_MakePrism(section, tangent_vec.Scaled(1000)).Shape()
        prism_cutter = shape_fix_shape(prism_cutter)
        prism_cutters.append(prism_cutter)

        cutters_for_cut.append(_cutter_cells)  # list of cells

        common_op = BRepAlgoAPI_Common(torus, cutter)
        common_op.Build()
        common_shape = common_op.Shape() if common_op.IsDone() else None
        common_volume = 0.0
        if common_shape is not None and not common_shape.IsNull():
            common_shape = shape_fix_shape(common_shape)
            common_volume = shape_volume(common_shape)
            common_shapes.append(common_shape)
        common_volumes.append(common_volume)

        cutter_volume = shape_volume(cutter)
        print(
            f"[DEBUG]   wire {i + 1}/{wires} cutter built ({len(_cutter_cells)} cells, type={cutter.ShapeType()}), total_cell_vol={total_cell_vol:.3f}, common volume={common_volume:.6f}"
        )
        sweep_list.append(cutter)

    result = torus
    print(f"[DEBUG] before cuts: torus solids={count_solids(result)}, ShapeType={result.ShapeType()}")
    for idx, cells in enumerate(cutters_for_cut, start=1):
        n_cell_ok = 0
        for k, cell in enumerate(cells):
            cut_op = BRepAlgoAPI_Cut(result, cell)
            cut_op.Build()
            if not cut_op.IsDone():
                continue
            new_result = cut_op.Shape()
            if new_result is None or new_result.IsNull():
                continue
            new_result = shape_fix_shape(new_result)
            n_solids = count_solids(new_result)
            if n_solids == 0 and new_result.ShapeType() == 0:
                continue
            if n_solids == 1 and new_result.ShapeType() == 0:
                single = next(TopologyExplorer(new_result).solids(), None)
                if single:
                    result = shape_fix_shape(single)
                    n_cell_ok += 1
                    continue
            if n_solids > 1:
                best = choose_largest_solid(new_result)
                result = shape_fix_shape(best) if best else new_result
                n_cell_ok += 1
                continue
            result = new_result
            n_cell_ok += 1
        print(f"[DEBUG] wire_{idx}: {n_cell_ok}/{len(cells)} cells cut, result solids={count_solids(result)}, ShapeType={result.ShapeType()}")

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
    torus_vol = shape_volume(torus)
    final_vol = shape_volume(final_shape)
    print(f"[DEBUG] torus_volume : {torus_vol:.1f}")
    print(f"[DEBUG] final_volume : {final_vol:.1f}")
    print(f"[DEBUG] removed vol  : {torus_vol - final_vol:.1f} ({100*(torus_vol-final_vol)/torus_vol:.2f}%)")
    for idx, wire in enumerate(wire_list, start=1):
        obj.display.DisplayShape(wire, color="RED")
    for idx, section in enumerate(section_list, start=1):
        obj.display.DisplayShape(section, transparency=0.7, color="BLUE1")
    for idx, sweep in enumerate(sweep_list, start=1):
        obj.display.DisplayShape(sweep, transparency=0.2, color="GREEN")
    #for idx, common_shape in enumerate(common_shapes, start=1):
    #    obj.display.DisplayShape(common_shape, transparency=0.5, color="BLUE1")
    obj.display.DisplayShape(torus, transparency=0.8, color="BLUE1")
    obj.display.DisplayShape(final_shape)
    obj.show_axs_pln(gp_Ax3(), scale=major_radius * 0.4)
    obj.save_view("debug")
    obj.ShowOCC()
