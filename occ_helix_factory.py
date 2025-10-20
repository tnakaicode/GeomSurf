#!/usr/bin/env python
# coding: utf-8

"""Helix factory utility for pythonOCC.

Provides a make_helix function that builds a helical edge (TopoDS_Edge)
from simple parameters: length, number of turns and radius.

The implementation samples points along a right-handed helix and fits
an interpolating Geom_BSplineCurve through them, then converts it to
a TopoDS_Edge for display or further modeling.

Example:
    from helix_factory import make_helix
    edge = make_helix(10.0, 5.0, 1.0)

"""

from math import pi, cos, sin, acos, asin, floor, ceil
from typing import Optional
import faulthandler
import traceback

from OCC.Core.gp import gp_Pnt
from OCC.Core.TColgp import TColgp_Array1OfPnt
from OCC.Core.GeomAPI import GeomAPI_PointsToBSpline
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge
from OCC.Display.SimpleGui import init_display
from OCC.Core.Geom import Geom_CylindricalSurface
from OCC.Core.GCE2d import GCE2d_MakeSegment
from OCC.Core.gp import gp_Ax3, gp_Pnt2d, gp_Lin2d, gp_Dir2d
from OCCUtils.Construct import make_box
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Section
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_VERTEX
from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeVertex

# enable faulthandler to get C-level tracebacks on crashes
faulthandler.enable()


def make_helix(length: float, turns: float, radius: float = 1.0, num_points: int = 100):
    """Create a helical edge.

    Args:
        length: total height (along Z) of the helix.
        turns: number of full revolutions (can be float).
        radius: radius of the helix (default 1.0).
        num_points: number of sample points used to build the spline.

    Returns:
        TopoDS_Edge representing the helix curve.

    Notes:
        - The helix is centered on the Z axis, starting at z=0 and
          ascending to z=length.
        - The function builds a BSpline through sampled points; for
          production geometry you may want to tune interpolation
          parameters or construct an exact helix via a threaded BSpline.
    """
    if num_points < 2:
        raise ValueError("num_points must be >= 2")

    pts = TColgp_Array1OfPnt(1, num_points)
    for i in range(num_points):
        t = i / (num_points - 1)
        theta = 2.0 * pi * turns * t
        x = radius * cos(theta)
        y = radius * sin(theta)
        z = length * t
        pts.SetValue(i + 1, gp_Pnt(x, y, z))

    spline_builder = GeomAPI_PointsToBSpline(pts)
    spline = spline_builder.Curve()

    edge = BRepBuilderAPI_MakeEdge(spline).Edge()
    return edge


def make_helix_on_cylinder(
    length: float, turns: float, radius: float = 1.0, start_angle: float = 0.0
):
    """Create an exact helix edge by mapping a 2D segment onto a cylindrical surface.

    This follows the approach used in `core_helix.py`: build a Geom_CylindricalSurface,
    create a 2D segment in the cylinder's parametric (u,v) space that spans the
    angular range corresponding to the requested turns, then make an edge on the
    surface between the parametric bounds. The resulting edge lies exactly on
    the cylinder (no spline interpolation).
    """
    # Create a cylindrical surface centered/oriented with default Ax3
    cyl = Geom_CylindricalSurface(gp_Ax3(), radius)

    # Build a 2D straight segment in the cylinder's parameter space where
    # the first coordinate is the angular parameter (u) and the second is z (v).
    # This decouples angle and height: u runs from start_u..end_u, v runs 0..length.
    start_u = start_angle
    end_u = start_angle + 2.0 * pi * turns
    p1 = gp_Pnt2d(start_u, 0.0)
    p2 = gp_Pnt2d(end_u, length)
    seg = GCE2d_MakeSegment(p1, p2)

    # Create an edge by mapping the 2D segment onto the cylindrical surface.
    # Use the overload that accepts (Geom2d_Curve, Geom_Surface).
    edge = BRepBuilderAPI_MakeEdge(seg.Value(), cyl).Edge()
    return edge


# compute intersections between helix edge and box
def intersect_edge_with_shape(edge, shape, tol: float = 1e-7):
    print("DEBUG: about to construct BRepAlgoAPI_Section")
    try:
        sec = BRepAlgoAPI_Section(edge, shape)
    except Exception as e:
        print("BRepAlgoAPI_Section constructor failed:", e)
        traceback.print_exc()
        return []

    print(
        "DEBUG: constructed Section object; calling ComputePCurveOn1/Approximation/SetFuzzyValue/Build"
    )
    try:
        sec.ComputePCurveOn1(True)
        sec.Approximation(True)
        sec.SetFuzzyValue(tol)
        sec.Build()
    except Exception as e:
        print("BRepAlgoAPI_Section.Build() failed:", e)
        traceback.print_exc()
        return []

    try:
        if not sec.IsDone():
            print("BRepAlgoAPI_Section reported not done")
            return []

        res = sec.Shape()
    except Exception as e:
        print("Error retrieving section shape:", e)
        return []

    pts = []
    exp = TopExp_Explorer(res, TopAbs_VERTEX)
    while exp.More():
        v = exp.Current()
        try:
            p = BRep_Tool.Pnt(v)
            pts.append(p)
        except Exception:
            # skip invalid vertex
            pass
        exp.Next()

    return pts


def intersect_helix_with_box(
    length: float,
    turns: float,
    radius: float,
    start_angle: float,
    box_min,
    box_max,
    tol=1e-9,
):
    """Analytically compute intersection points between a helix on a cylinder
    (as produced by make_helix_on_cylinder) and an axis-aligned box.

    This avoids native intersection routines by solving for parameter u (angle)
    where the helix meets box faces (x = const, y = const, z = const).

    box_min, box_max: sequences or gp_Pnt giving (xmin,ymin,zmin) and (xmax,ymax,zmax)
    Returns a list of tuples (gp_Pnt, face_id) where face_id is one of
    'x0','x1','y0','y1','z0','z1'.
    """
    xmin, ymin, zmin = (
        (box_min.X(), box_min.Y(), box_min.Z())
        if hasattr(box_min, "X")
        else tuple(box_min)
    )
    xmax, ymax, zmax = (
        (box_max.X(), box_max.Y(), box_max.Z())
        if hasattr(box_max, "X")
        else tuple(box_max)
    )

    u0 = start_angle
    u1 = start_angle + 2.0 * pi * turns
    results = []

    def in_urange(u):
        return (u + 1e-12) >= u0 and (u - 1e-12) <= u1

    def z_of_u(u):
        if abs(u1 - u0) < 1e-12:
            return 0.0
        return length * (u - u0) / (u1 - u0)

    # helper to add candidate u (normalized by 2pi shifts)
    def add_u_candidate(u_base, face_id, x_check=None, y_check=None, z_check=None):
        # generate k such that u = u_base + 2pi*k is in [u0,u1]
        two_pi = 2.0 * pi
        kmin = int(ceil((u0 - u_base) / two_pi))
        kmax = int(floor((u1 - u_base) / two_pi))
        for k in range(kmin, kmax + 1):
            u = u_base + two_pi * k
            if not in_urange(u):
                continue
            z = z_of_u(u)
            x = radius * cos(u)
            y = radius * sin(u)
            ok = True
            if x_check is not None:
                if not (x_check[0] - tol <= x <= x_check[1] + tol):
                    ok = False
            if y_check is not None:
                if not (y_check[0] - tol <= y <= y_check[1] + tol):
                    ok = False
            if z_check is not None:
                if not (z_check[0] - tol <= z <= z_check[1] + tol):
                    ok = False
            if ok:
                results.append((gp_Pnt(x, y, z), face_id))

    # intersections with x = xmin and x = xmax (if within radius)
    for x_val, face in ((xmin, "x0"), (xmax, "x1")):
        if abs(x_val) <= abs(radius) + 1e-12:
            c = x_val / radius
            if c >= -1.0 - 1e-12 and c <= 1.0 + 1e-12:
                try:
                    a = acos(max(-1.0, min(1.0, c)))
                except Exception:
                    continue
                # two base solutions: +a and -a
                add_u_candidate(
                    a,
                    face,
                    x_check=(x_val, x_val),
                    y_check=(ymin, ymax),
                    z_check=(zmin, zmax),
                )
                add_u_candidate(
                    -a,
                    face,
                    x_check=(x_val, x_val),
                    y_check=(ymin, ymax),
                    z_check=(zmin, zmax),
                )

    # intersections with y = ymin/ymax via asin
    for y_val, face in ((ymin, "y0"), (ymax, "y1")):
        if abs(y_val) <= abs(radius) + 1e-12:
            s = y_val / radius
            if s >= -1.0 - 1e-12 and s <= 1.0 + 1e-12:
                try:
                    a = asin(max(-1.0, min(1.0, s)))
                except Exception:
                    continue
                add_u_candidate(
                    a,
                    face,
                    y_check=(y_val, y_val),
                    x_check=(xmin, xmax),
                    z_check=(zmin, zmax),
                )
                add_u_candidate(
                    pi - a,
                    face,
                    y_check=(y_val, y_val),
                    x_check=(xmin, xmax),
                    z_check=(zmin, zmax),
                )

    # intersections with z = zmin/zmax: solve for u from linear relation
    if abs(length) > 1e-12:
        for z_val, face in ((zmin, "z0"), (zmax, "z1")):
            # u = u0 + (u1-u0) * (z_val/length)
            frac = z_val / length
            u = u0 + (u1 - u0) * frac
            if in_urange(u):
                x = radius * cos(u)
                y = radius * sin(u)
                if (xmin - tol <= x <= xmax + tol) and (ymin - tol <= y <= ymax + tol):
                    results.append((gp_Pnt(x, y, z_val), face))

    return results


if __name__ == "__main__":
    display, start_display, _, _ = init_display()
    length = 10.0
    turns = 100.1
    radius = 5.0
    start_angle = 0.0
    e = make_helix_on_cylinder(length, turns, radius=radius, start_angle=start_angle)
    box = make_box(gp_Pnt(0, 0, 0), gp_Pnt(10, 10, 10))
    display.DisplayShape(e, update=True)
    display.DisplayShape(box, transparency=0.9)

    # analytic intersection (helix on cylinder vs axis-aligned box)
    hits = intersect_helix_with_box(
        length, turns, radius, start_angle, gp_Pnt(0, 0, 0), gp_Pnt(10, 10, 10)
    )
    print(f"analytic hits: {len(hits)}")
    for p, face in hits:
        display.DisplayShape(p)
    display.FitAll()
    start_display()
