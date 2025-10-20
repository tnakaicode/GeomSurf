"""
Helpers for creating a B-Spline surface with pythonocc and projecting a Wire (or point sequence)
onto it. Provides utilities to close the projected wire and to clamp/replace out-of-bounds
projected points using the surface boundary when necessary.

Notes:
- This module assumes pythonocc-core is installed and available as OCP/occ libraries.
- The projection strategy: for each point on the input wire we compute the closest point on
  the surface using GeomAPI_ProjectPointOnSurf. If the projected (u,v) parameter is outside
  the surface parametric domain, we clamp to the domain and optionally interpolate along
  the surface boundary to keep the resulting curve closed.

Functions:
- make_bspline_surface(points_u, points_v, degree_u=3, degree_v=3): create a surface
- project_wire_to_surface(wire, surface, tol=1e-6): project edges/points of a wire
- close_projected_wire(wire_on_surface, surface, tol=1e-6): ensure the projected wire is closed

"""

# ---------------------------------------------------------------------------
# Module overview
#
# English:
# This module provides utilities to build a B-spline surface and to project a
# TopoDS_Wire (or sampled points from a wire) onto that surface. It includes
# helpers to clamp projections to the surface domain, to close projected wires
# (with simple heuristics), and some utility routines to safely create edges
# and to assemble hybrid wires (polygon + arc).
#
# The typical flow is:
#  1. Build a Geom_BSplineSurface from a grid of 3D points (make_bspline_surface).
#  2. For each edge of an input wire, sample points along the edge and project
#     those samples onto the surface (project_wire_to_surface). Sampling avoids
#     "tabbing" artifacts caused when only vertices are projected and long
#     straight connectors are created.
#  3. If the projected wire is not geometrically closed, attempt to close it
#     by interpolating along a surface boundary when the endpoints map to the
#     same parametric boundary (close_projected_wire), otherwise connect with
#     a straight segment as a fallback.
#  4. Helpers include safe_make_edge (avoid degenerate edges), polyline
#     self-intersection detection, arc sampling, and functions to insert
#     arc endpoints into polygon point lists so a single closed hybrid wire can
#     be assembled.
#
# How to resolve duplicates:
# - Frequently duplication arises from multiple fallback implementations or
#   repeated import lines. To clean up:
#   1) Consolidate imports at the top of the file and remove duplicates. Keep
#      only one canonical import per symbol (e.g. one GeomAPI_ProjectPointOnSurf).
#   2) If the file contains duplicated functions (earlier and later versions),
#      keep the most robust variant and remove the older one. Use unit tests or
#      a quick smoke run to ensure the kept variant passes.
#   3) Factor repeated logic into small helpers (e.g. a single safe_make_edge,
#      a single projection wrapper) and call these helpers instead of repeating
#      inline code blocks.
#   4) For API-version differences, centralize the fallbacks into small
#      adapter helpers (e.g. get_surface_uv_bounds) instead of scattering
#      try/excepts throughout the code.
#
# Self-intersection avoidance Strategies to avoid/repair:
# - Prevention (primary):
#   * Use dense, adaptive sampling on input edges before projection. Increase
#     sampling where curvature or param change is large so the projected polyline
#     follows the surface mapping and does not produce long chords that cut in.
#   * Preserve sampling order and, where possible, use the surface's (u,v)
#     parameter continuity to detect and unwrap domain wrap-around (periodic
#     surfaces). Sorting or reordering by parameter can avoid unnatural jumps.
#   * If an edge maps partially outside the domain, consider projecting to the
#     boundary and then sampling along the boundary curve (iso-curve) instead
#     of clamping isolated points; this produces a topologically consistent
#     seam.
#
# - Detection & Repair (secondary):
#   * Use the provided polyline_self_intersects() helper (or an external 2D
#     geometry library like shapely) to detect intersections in candidate
#     point sequences. When detected:
#       - Try alternate splice directions (as implemented: candidate A/B) and
#         pick the non-intersecting one.
#       - If both candidates self-intersect, attempt local refinement: insert
#         additional samples around the problematic region and re-check.
#       - As a last resort, run a wire-fixing routine such as OCC's
#         ShapeFix_Wire (if available) or convert the polyline to a trimmed
#         face and re-extract its outer wire.
#
# - Stronger approach (recommended for robust workflows):
#   * Instead of connecting projected 3D sample points by straight edges, build
#     a 2D curve in UV-space: project input edge samples to (u,v) pairs, build a
#     smooth 2D interpolating curve (Geom2dAPI_Interpolate) in the parametric
#     domain, and then create a PCurve on the surface. This keeps the topology
#     on the surface and avoids many self-intersection issues that arise from
#     naive 3D straight-line stitching.
#   * If the surface boundary must be used to close a wire, extract the
#     boundary curve(s) with ShapeAnalysis or Surface methods and splice using
#     parametric interpolation rather than straight chords.
#
# Practical suggestions:
# - Add small diagnostic logging: write per-sample (u,v) sequences and 3D
#   distances to a CSV for problematic wires so you can see where large jumps
#   occur and tune sampling thresholds.
# - Add unit tests: create small shapes where you know the expected projected
#   result, then assert wires are closed and non-self-intersecting.
# - Keep fallbacks centralized and remove duplicated logic after verification.
#
# ---------------------------------------------------------------------------

from __future__ import annotations

# Try to import pythonocc (OCC.Core or OCP). Support both common namespaces.
from OCC.Core.gp import gp_Pnt, gp_Vec
from OCC.Core.TColgp import TColgp_Array2OfPnt
from OCC.Core.GeomAPI import GeomAPI_PointsToBSplineSurface, GeomAPI_ProjectPointOnSurf
from OCC.Core.TColStd import TColStd_Array1OfReal
from OCC.Core.GeomAPI import GeomAPI_ProjectPointOnSurf
from OCC.Core.GeomAPI import GeomAPI_ExtremaCurveSurface
from OCC.Core.Geom import Geom_BSplineSurface
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeWire, BRepBuilderAPI_MakeEdge
from OCC.Core.TopoDS import TopoDS_Wire, topods
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_VERTEX, TopAbs_EDGE
from OCC.Core.BRep import BRep_Tool
from OCC.Core.ShapeAnalysis import ShapeAnalysis_Surface
from OCC.Core.GeomLProp import GeomLProp_SLProps
from OCC.Core.Geom2dAPI import Geom2dAPI_Interpolate

from typing import List, Tuple, Optional
import math


def make_bspline_surface(
    points_u: List[List[Tuple[float, float, float]]],
    degree_u: int = 3,
    degree_v: int = 3,
) -> Geom_BSplineSurface:
    """
    Create a B-spline surface from a 2D grid of points.

    points_u: list of rows (v direction), each row is a list of (x,y,z) points along u.
    degree_u, degree_v: polynomial degrees in U and V.

    Returns a Geom_BSplineSurface.
    """
    rows = len(points_u)
    if rows == 0:
        raise ValueError("points_u must contain at least one row")
    cols = len(points_u[0])
    arr = TColgp_Array2OfPnt(1, rows, 1, cols)
    for i in range(rows):
        if len(points_u[i]) != cols:
            raise ValueError("All rows must have the same number of columns")
        for j in range(cols):
            x, y, z = points_u[i][j]
            arr.SetValue(i + 1, j + 1, gp_Pnt(x, y, z))

    # GeomAPI_PointsToBSplineSurface constructor signatures vary between versions;
    # use the common (array, degU, degV) form and then call Surface().
    try:
        builder = GeomAPI_PointsToBSplineSurface(arr, degree_u, degree_v)
    except TypeError:
        # older/newer variants may accept different arg orders; fall back to single-arg style
        builder = GeomAPI_PointsToBSplineSurface(arr)
    # GeomAPI_PointsToBSplineSurface provides Surface() returning Geom_Surface
    surf = builder.Surface()
    if not isinstance(surf, Geom_BSplineSurface):
        # Attempt to downcast if necessary
        surf = Geom_BSplineSurface.DownCast(surf)
    return surf


def project_point_to_surface(
    p: gp_Pnt, surface: Geom_BSplineSurface, tol: float = 1e-6
) -> Tuple[gp_Pnt, float, float]:
    """
    Project a single 3D point onto the surface. Returns the projected 3D point and the (u, v)
    parameters on the surface. If projection falls outside param domain, parameter values may be
    outside the range — caller should clamp if needed.
    """
    projector = GeomAPI_ProjectPointOnSurf(p, surface)
    if not projector.IsDone() or projector.NbPoints() == 0:
        raise RuntimeError("Projection failed")
    uv_pnt = projector.Point(1)
    # Use the SWIG-exposed Parameters method which returns (u, v)
    params = projector.Parameters(1)
    # Parameters may return a tuple-like sequence
    try:
        u, v = params[0], params[1]
    except Exception:
        # If it's already two values
        u, v = params
    return uv_pnt, u, v


def clamp_uv(u: float, v: float, surface: Geom_BSplineSurface) -> Tuple[float, float]:
    """
    Clamp (u,v) to the surface parametric domain.
    """
    u1, u2, v1, v2 = get_surface_uv_bounds(surface)
    cu = max(min(u, u2), u1)
    cv = max(min(v, v2), v1)
    return cu, cv


def get_surface_uv_bounds(surface) -> Tuple[float, float, float, float]:
    """
    Robustly obtain the surface parametric bounds (u1, u2, v1, v2).
    Different pythonocc/OCCT bindings expose different method names; try several fallbacks.
    """
    # common method names in various bindings
    candidates = [
        ("U1", "U2", "V1", "V2"),
        ("FirstUParameter", "LastUParameter", "FirstVParameter", "LastVParameter"),
        ("Bounds",),  # may return (u1,u2,v1,v2) or fill refs
    ]

    for cand in candidates:
        try:
            if cand[0] == "Bounds":
                # try calling Bounds()
                b = surface.Bounds()
                # Bounds might return a tuple
                if isinstance(b, (tuple, list)) and len(b) >= 4:
                    return float(b[0]), float(b[1]), float(b[2]), float(b[3])
                # some wrappers may fill by reference via arguments - skip here
                continue
            u1 = getattr(surface, cand[0])()
            u2 = getattr(surface, cand[1])()
            v1 = getattr(surface, cand[2])()
            v2 = getattr(surface, cand[3])()
            return float(u1), float(u2), float(v1), float(v2)
        except Exception:
            continue

    # Fallback: use ShapeAnalysis_Surface to estimate param bounds by sampling
    try:
        sas = ShapeAnalysis_Surface(surface)
        # ShapeAnalysis_Surface has methods to get parameter range in some bindings
        getu1 = getattr(sas, "UFirst", None) or getattr(sas, "FirstU", None)
        getu2 = getattr(sas, "ULast", None) or getattr(sas, "LastU", None)
        getv1 = getattr(sas, "VFirst", None) or getattr(sas, "FirstV", None)
        getv2 = getattr(sas, "VLast", None) or getattr(sas, "LastV", None)
        if callable(getu1) and callable(getu2) and callable(getv1) and callable(getv2):
            return float(getu1()), float(getu2()), float(getv1()), float(getv2())
    except Exception:
        pass
    # As ultimate fallback return a large domain to avoid crashes; caller should notice if wrong
    return -1e6, 1e6, -1e6, 1e6


def project_wire_to_surface(
    wire: TopoDS_Wire, surface: Geom_BSplineSurface, tol: float = 1e-6
) -> TopoDS_Wire:
    """
    Project every vertex of a TopoDS_Wire onto the surface and return a new wire made of
    edges connecting the projected points (in the same order). If a projected point is outside
    the parametric domain, it will be clamped to the boundary.
    """
    # Iterate edges and sample each edge, projecting sampled points to the surface.
    exp_e = TopExp_Explorer(wire, TopAbs_EDGE)
    all_pts: List[gp_Pnt] = []
    sas = ShapeAnalysis_Surface(surface)
    max_samples_per_edge = 200
    while exp_e.More():
        edge = topods.Edge(exp_e.Current())
        # try to get underlying 3d curve for the edge
        try:
            curve_data = BRep_Tool.Curve(edge)
        except Exception:
            curve_data = None

        pts_on_edge: List[gp_Pnt] = []
        if curve_data:
            try:
                curve, f, l = curve_data
                # estimate length by sampling few points
                sample_len = 10
                prev_p = gp_Pnt()
                curve.D0(f, prev_p)
                est_len = 0.0
                for k in range(1, sample_len + 1):
                    t = f + (l - f) * k / sample_len
                    ptmp = gp_Pnt()
                    curve.D0(t, ptmp)
                    est_len += prev_p.Distance(ptmp)
                    prev_p = ptmp
                # choose number of samples proportional to length (one sample per ~3.0 units)
                target = max(2, min(max_samples_per_edge, int(est_len / 3.0) + 1))
                for k in range(target + 1):
                    t = f + (l - f) * k / target
                    ptmp = gp_Pnt()
                    curve.D0(t, ptmp)
                    pts_on_edge.append(ptmp)
            except Exception:
                pts_on_edge = []
        # fallback: if no curve, sample edge vertices
        if not pts_on_edge:
            exp_v = TopExp_Explorer(edge, TopAbs_VERTEX)
            while exp_v.More():
                v = topods.Vertex(exp_v.Current())
                pts_on_edge.append(BRep_Tool.Pnt(v))
                exp_v.Next()

        # project each sampled point on this edge to surface
        for p in pts_on_edge:
            try:
                proj = GeomAPI_ProjectPointOnSurf(p, surface)
                if not proj.IsDone() or proj.NbPoints() == 0:
                    continue
                # get parameters robustly
                params = proj.Parameters(1)
                try:
                    u, vpar = params[0], params[1]
                except Exception:
                    u, vpar = params
                cu, cv = clamp_uv(u, vpar, surface)
                surf_p = sas.Value(cu, cv)
                # avoid adding nearly-duplicate points
                if not all_pts or surf_p.Distance(all_pts[-1]) > tol:
                    all_pts.append(surf_p)
            except Exception:
                continue

        exp_e.Next()

    # Build wire from projected points (connect sequentially)
    mk = BRepBuilderAPI_MakeWire()
    for i in range(len(all_pts)):
        a = all_pts[i]
        b = all_pts[(i + 1) % len(all_pts)]
        e = safe_make_edge(a, b, tol=1e-6)
        if e is not None:
            try:
                mk.Add(e)
            except Exception:
                pass
    return mk.Wire()


def close_projected_wire(
    wire_on_surface: TopoDS_Wire, surface: Geom_BSplineSurface, tol: float = 1e-6
) -> TopoDS_Wire:
    """
    Ensure the given wire (assumed to lie on or near the surface) is geometrically closed.
    If endpoints are separated beyond tol, attempt to connect them along the nearest surface
    boundary: sample the surface boundary near the endpoints and interpolate a closing edge.

    This is a heuristic function: it will try to detect whether the wire wraps around the
    parametric domain edges and will connect endpoints along that boundary if possible.
    """
    # Collect vertices
    exp = TopExp_Explorer(wire_on_surface, TopAbs_VERTEX)
    verts: List[gp_Pnt] = []
    while exp.More():
        v = topods.Vertex(exp.Current())
        p = BRep_Tool.Pnt(v)
        verts.append(p)
        exp.Next()
    if len(verts) == 0:
        raise ValueError("Empty wire")

    # Check closure
    start = verts[0]
    end = verts[-1]
    d = start.Distance(end)
    if d <= tol:
        return wire_on_surface
    # Compute projections of start/end on surface to get their (u,v)
    proj = GeomAPI_ProjectPointOnSurf
    p1proj = proj(start, surface)
    p2proj = proj(end, surface)
    if (
        not p1proj.IsDone()
        or not p2proj.IsDone()
        or p1proj.NbPoints() == 0
        or p2proj.NbPoints() == 0
    ):
        # fallback: copy existing edges and connect with straight edge
        mk = BRepBuilderAPI_MakeWire()
        exp_e = TopExp_Explorer(wire_on_surface, TopAbs_EDGE)
        while exp_e.More():
            e = topods.Edge(exp_e.Current())
            mk.Add(e)
            exp_e.Next()
        e = safe_make_edge(start, end)
        if e is not None:
            mk.Add(e)
        return mk.Wire()

    params1 = p1proj.Parameters(1)
    params2 = p2proj.Parameters(1)
    try:
        u1, v1 = params1[0], params1[1]
    except Exception:
        u1, v1 = params1
    try:
        u2, v2 = params2[0], params2[1]
    except Exception:
        u2, v2 = params2

    u1c, v1c = clamp_uv(u1, v1, surface)
    u2c, v2c = clamp_uv(u2, v2, surface)

    umin, umax, vmin, vmax = get_surface_uv_bounds(surface)

    def on_u_boundary(u, v):
        return abs(u - umin) < 1e-9 or abs(u - umax) < 1e-9

    def on_v_boundary(u, v):
        return abs(v - vmin) < 1e-9 or abs(v - vmax) < 1e-9

    # Build the resulting wire by copying existing edges
    mk = BRepBuilderAPI_MakeWire()
    exp_e = TopExp_Explorer(wire_on_surface, TopAbs_EDGE)
    while exp_e.More():
        e = topods.Edge(exp_e.Current())
        mk.Add(e)
        exp_e.Next()

    # If both endpoints lie on the same parametric boundary (u=const or v=const), interpolate
    # points along that boundary between the two parameter positions.
    if on_u_boundary(u1c, v1c) and on_u_boundary(u2c, v2c) and abs(u1c - u2c) < 1e-9:
        uconst = u1c
        t1 = v1c
        t2 = v2c
        # sample N points along boundary
        N = max(4, int(abs(t2 - t1) / ((vmax - vmin) / 20)))
        pts = []
        for i in range(N + 1):
            t = t1 + (t2 - t1) * i / N
            p = surface.Value(uconst, t)
            pts.append(p)
        # create edges along pts
        for i in range(len(pts) - 1):
            e = safe_make_edge(pts[i], pts[i + 1], tol=1e-9)
            if e is not None:
                mk.Add(e)
        return mk.Wire()

    if on_v_boundary(u1c, v1c) and on_v_boundary(u2c, v2c) and abs(v1c - v2c) < 1e-9:
        vconst = v1c
        t1 = u1c
        t2 = u2c
        N = max(4, int(abs(t2 - t1) / ((umax - umin) / 20)))
        pts = []
        for i in range(N + 1):
            t = t1 + (t2 - t1) * i / N
            p = surface.Value(t, vconst)
            pts.append(p)
        for i in range(len(pts) - 1):
            e = safe_make_edge(pts[i], pts[i + 1], tol=1e-9)
            if e is not None:
                mk.Add(e)
        return mk.Wire()

    # Otherwise, if endpoints are on different boundaries or interior, just add straight edge
    e = safe_make_edge(start, end)
    if e is not None:
        mk.Add(e)
    return mk.Wire()


def safe_make_edge(a: gp_Pnt, b: gp_Pnt, tol: float = 1e-9):
    """
    Try to make an edge between two gp_Pnt safely. Returns TopoDS_Edge or None.
    Avoid making edges when points are coincident or extremely close (which can raise StdFail_NotDoneBRep_API).
    """
    try:
        if a.Distance(b) <= tol:
            return None
        mk = BRepBuilderAPI_MakeEdge(a, b)
        if not mk.IsDone():
            return None
        return mk.Edge()
    except Exception:
        return None


def insert_point_into_poly_pts(poly_pts: List[gp_Pnt], pt: gp_Pnt, tol: float = 1e-6):
    """Insert pt into poly_pts list at the appropriate segment if it projects onto a segment
    (within tol). Returns (new_list, inserted_index).
    If no suitable segment projection, insert after nearest vertex.
    """
    best_seg = None
    best_t = None
    best_dist = float("inf")
    for i in range(len(poly_pts)):
        a = poly_pts[i]
        b = poly_pts[(i + 1) % len(poly_pts)]
        # vector from a to b
        vx = b.X() - a.X()
        vy = b.Y() - a.Y()
        vz = b.Z() - a.Z()
        wx = pt.X() - a.X()
        wy = pt.Y() - a.Y()
        wz = pt.Z() - a.Z()
        denom = vx * vx + vy * vy + vz * vz
        if denom == 0:
            continue
        t = (vx * wx + vy * wy + vz * wz) / denom
        if t < 0.0 or t > 1.0:
            continue
        # closest point on segment
        cx = a.X() + t * vx
        cy = a.Y() + t * vy
        cz = a.Z() + t * vz
        dx = pt.X() - cx
        dy = pt.Y() - cy
        dz = pt.Z() - cz
        dist2 = dx * dx + dy * dy + dz * dz
        if dist2 < best_dist:
            best_dist = dist2
            best_seg = i
            best_t = t
    if best_seg is not None and best_dist <= tol * tol:
        # insert new gp_Pnt at position best_seg+1
        new_pt = gp_Pnt(pt.X(), pt.Y(), pt.Z())
        new_list = poly_pts[: best_seg + 1] + [new_pt] + poly_pts[best_seg + 1 :]
        return new_list, best_seg + 1
    # fallback: insert after nearest vertex
    best_v = min(range(len(poly_pts)), key=lambda i: poly_pts[i].Distance(pt))
    new_pt = gp_Pnt(pt.X(), pt.Y(), pt.Z())
    new_list = poly_pts[: best_v + 1] + [new_pt] + poly_pts[best_v + 1 :]
    return new_list, best_v + 1


def sample_arc_edge_to_points(arc_edge, samples: int = 24):
    """Return a list of gp_Pnt sampling the Geom curve underlying arc_edge."""
    pts = []
    try:
        # try to get curve and params
        curve, f, l = BRep_Tool.Curve(topods.Edge(arc_edge))
        if curve is None:
            return pts
        for i in range(samples + 1):
            t = f + (l - f) * i / samples
            p = gp_Pnt()
            curve.D0(t, p)
            pts.append(p)
    except Exception:
        # fallback: sample by vertices of edge
        exp = TopExp_Explorer(arc_edge, TopAbs_VERTEX)
        while exp.More():
            v = topods.Vertex(exp.Current())
            pts.append(BRep_Tool.Pnt(v))
            exp.Next()
    return pts


def polyline_self_intersects(points: List[gp_Pnt]) -> bool:
    """Detect self-intersection in closed polyline (points in order). Uses 2D (X,Y).
    Excludes adjacent segments from checks.
    """

    def seg_intersect(a1, a2, b1, b2):
        x1, y1 = a1.X(), a1.Y()
        x2, y2 = a2.X(), a2.Y()
        x3, y3 = b1.X(), b1.Y()
        x4, y4 = b2.X(), b2.Y()
        denom = (y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1)
        if abs(denom) < 1e-12:
            return False
        ua = ((x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3)) / denom
        ub = ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)) / denom
        return 0.0 < ua < 1.0 and 0.0 < ub < 1.0

    m = len(points)
    if m < 4:
        return False
    for i in range(m):
        a1 = points[i]
        a2 = points[(i + 1) % m]
        for j in range(i + 2, m):
            # skip adjacent and same segment wrap
            if j == i or (j + 1) % m == i:
                continue
            b1 = points[j]
            b2 = points[(j + 1) % m]
            if seg_intersect(a1, a2, b1, b2):
                return True
    return False


# Minimal usage example in module docstring
if __name__ == "__main__":
    # Runnable example using pythonocc display
    from OCC.Display.SimpleGui import init_display
    from OCC.Core.gp import gp_Pnt, gp_Ax2, gp_Dir, gp_Circ
    from OCC.Core.GC import GC_MakeArcOfCircle
    import math as _math

    # build a simple grid of points (rows in V, cols in U) forming a slightly curved patch
    rows = 9
    cols = 8
    grid = []
    for i in range(rows):
        row = []
        for j in range(cols):
            x = j * (rows + 2)
            y = i * (cols + 3)
            z = math.sin(j * 0.45) * math.cos(i * 0.25) * 4.0
            row.append((x, y, z))
        grid.append(row)

    surf = make_bspline_surface(grid, degree_u=3, degree_v=3)

    # Helper: make a polygonal ring (outer polygon) with center offset
    def make_polygon(center_x, center_y, radii, z=0.0):
        pts = []
        n = len(radii)
        for i in range(n):
            ang = 2.0 * _math.pi * i / n
            r = radii[i]
            x = center_x + r * _math.cos(ang)
            y = center_y + r * _math.sin(ang)
            pts.append(gp_Pnt(x, y, z))
        mk = BRepBuilderAPI_MakeWire()
        for i in range(len(pts)):
            a = pts[i]
            b = pts[(i + 1) % len(pts)]
            e = safe_make_edge(a, b)
            if e is not None:
                mk.Add(e)
        return mk.Wire(), pts

    # Helper: make a circular arc as an edge (from angle1 to angle2)
    def make_arc(center_x, center_y, radius, angle1, angle2, z=0.0):
        c = gp_Circ(gp_Ax2(gp_Pnt(center_x, center_y, z), gp_Dir(0, 0, 1)), radius)
        # Use overload (gp_Circ, Standard_Real, Standard_Real, Standard_Boolean)
        try:
            arc_maker = GC_MakeArcOfCircle(c, angle1, angle2, True)
            if not arc_maker.IsDone():
                raise Exception("arc maker failed")
            curve = arc_maker.Value()
        except Exception:
            # Fallback: build arc from three points on the circle
            p1 = gp_Pnt(
                center_x + radius * _math.cos(angle1),
                center_y + radius * _math.sin(angle1),
                z,
            )
            pm = gp_Pnt(
                center_x + radius * _math.cos((angle1 + angle2) / 2.0),
                center_y + radius * _math.sin((angle1 + angle2) / 2.0),
                z,
            )
            p2 = gp_Pnt(
                center_x + radius * _math.cos(angle2),
                center_y + radius * _math.sin(angle2),
                z,
            )
            try:
                arc_maker2 = GC_MakeArcOfCircle(p1, pm, p2)
                if not arc_maker2.IsDone():
                    return None
                curve = arc_maker2.Value()
            except Exception:
                return None
        # create edge from curve
        try:
            from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge

            ek = BRepBuilderAPI_MakeEdge(curve)
            if ek.IsDone():
                return ek.Edge()
        except Exception:
            return None
        return None

    # Create two hybrid wires: one fully inside the surface region, one partially outside
    cx = 24.0
    cy = 18.0

    # inside: polygon + small arc inside
    poly_inside, poly_inside_pts = make_polygon(cx, cy, [12, 11, 12, 11, 12, 11])
    arc_edge = make_arc(cx, cy, 6.0, 0.0, _math.pi)

    # assemble hybrid: replace polygon segment between nearest polygon vertices with the arc
    def assemble_hybrid(poly_wire, poly_pts, arc_e):
        # If no arc or no polygon points, fallback
        if arc_e is None or not poly_pts:
            return poly_wire

        n = len(poly_pts)
        # sample arc to points
        arc_sample = sample_arc_edge_to_points(arc_e, samples=48)
        if not arc_sample:
            return poly_wire

        # insert arc endpoints into polygon vertices (split edges where appropriate)
        pts_with_first, idx_first = insert_point_into_poly_pts(poly_pts, arc_sample[0])
        # after inserting first, indices may shift; insert second into updated list
        pts_with_both, idx_second = insert_point_into_poly_pts(
            pts_with_first, arc_sample[-1]
        )

        # find indices in final list matching the inserted endpoints
        idx0 = idx_first
        idx1 = idx_second
        n2 = len(pts_with_both)

        # build helper to get keep sequence
        def build_keep_sequence(start, end):
            seq = []
            i = start
            while True:
                seq.append(i)
                if i == end:
                    break
                i = (i + 1) % n2
            return seq

        forward = (idx1 - idx0) % n2
        backward = (idx0 - idx1) % n2

        # candidate A: keep indices from idx1 -> idx0 and append arc_sample
        keepA = build_keep_sequence(idx1, idx0)
        ptsA = [pts_with_both[i] for i in keepA] + arc_sample

        # candidate B: keep indices from idx0 -> idx1 and append reversed arc_sample
        keepB = build_keep_sequence(idx0, idx1)
        ptsB = [pts_with_both[i] for i in keepB] + list(reversed(arc_sample))

        # prefer the candidate that does not self-intersect; if both OK pick the one that removes shorter segment
        okA = not polyline_self_intersects(ptsA)
        okB = not polyline_self_intersects(ptsB)

        candidate = None
        if okA and not okB:
            candidate = ptsA
        elif okB and not okA:
            candidate = ptsB
        elif okA and okB:
            # choose based on shorter removed segment
            if forward <= backward:
                candidate = ptsA
            else:
                candidate = ptsB
        else:
            # both intersect; pick the one that removes the shorter segment
            candidate = ptsA if forward <= backward else ptsB

        # build wire from candidate points
        mk = BRepBuilderAPI_MakeWire()
        for i in range(len(candidate)):
            a = candidate[i]
            b = candidate[(i + 1) % len(candidate)]
            e = safe_make_edge(a, b)
            if e is not None:
                try:
                    mk.Add(e)
                except Exception:
                    pass
        try:
            return mk.Wire()
        except Exception:
            return poly_wire

    hybrid_inside = assemble_hybrid(poly_inside, poly_inside_pts, arc_edge)

    # outside: larger polygon that extends beyond surface extents + arc
    poly_out, poly_out_pts = make_polygon(cx + 8.0, cy - 6.0, [40, 35, 40, 35, 40, 35])
    arc_edge2 = make_arc(cx + 8.0, cy - 6.0, 30.0, -0.5 * _math.pi, 0.5 * _math.pi)
    hybrid_out = assemble_hybrid(poly_out, poly_out_pts, arc_edge2)

    # Project and close both wires
    projected_in = project_wire_to_surface(hybrid_inside, surf)
    closed_in = close_projected_wire(projected_in, surf)

    projected_out = project_wire_to_surface(hybrid_out, surf)
    closed_out = close_projected_wire(projected_out, surf)

    # display
    display, start_display, add_menu, add_function_to_menu = init_display()
    display.DisplayShape(surf, update=True)
    display.DisplayShape(hybrid_inside, color="BLACK", update=True)
    display.DisplayShape(projected_in, color="BLUE1", update=True)
    display.DisplayShape(closed_in, color="RED", update=True)

    display.DisplayShape(hybrid_out, color="BLACK", update=True)
    display.DisplayShape(projected_out, color="GREEN1", update=True)
    display.DisplayShape(closed_out, color="YELLOW", update=True)

    print(
        "Displaying: surface, two hybrid input wires (black), projected (blue/green), closed (red/yellow)"
    )
    start_display()
