"""
plot_gaussian_slice.py

Create a Wire that represents a z-direction cross-section of z = Ggaussian(x, y)
and display it with pythonocc-core.

"""

import math
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import random

from OCC.Core.gp import gp_Pnt, gp_Pln, gp_Ax3, gp_Dir, gp_Vec
from OCC.Core.Geom import Geom_Plane
from OCC.Core.BRepBuilderAPI import (
    BRepBuilderAPI_MakePolygon,
    BRepBuilderAPI_MakeVertex,
)
from OCC.Display.SimpleGui import init_display
from OCC.Core.TColgp import TColgp_Array1OfPnt
from OCC.Core.GeomAPI import GeomAPI_PointsToBSpline
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire
from OCCUtils.Construct import make_polygon, make_edge, make_wire, make_face
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_ThruSections
from OCC.Core.GeomAPI import GeomAPI_Interpolate
from OCC.Core.TColgp import (
    TColgp_Array1OfPnt,
    TColgp_Array2OfPnt,
    TColgp_HArray1OfPnt,
    TColgp_HArray2OfPnt,
)
from OCC.Core.BRep import BRep_Builder
from OCC.Core.TopoDS import TopoDS_Compound
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_VERTEX, TopAbs_IN, TopAbs_WIRE
from OCC.Core.BRepTools import breptools
from OCC.Core.BRepTopAdaptor import BRepTopAdaptor_FClass2d
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
from OCC.Core.gp import gp_Pnt2d
from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepIntCurveSurface import BRepIntCurveSurface_Inter
from OCC.Core.BRepLProp import BRepLProp_SLProps
from OCC.Core.gp import gp_Lin
from OCC.Core.BRepProj import BRepProj_Projection


def Ggaussian(x, y, A=1.0, sx=1.0, sy=1.0, cx=0.0, cy=0.0, alpha=1.0):
    """A 2D Gaussian centered at (cx, cy).

    z = A * exp(-(((x-cx)/sx)**2 + ((y-cy)/sy)**2))
    """
    return A * math.exp(-alpha * (((x - cx) / sx) ** 2 + ((y - cy) / sy) ** 2))


def make_slice_wire_xz(y0, x_range, n, smooth=False):
    """Return a Wire following (x, y0, Ggaussian(x,y0)) sampled along x.

    If smooth=True a B-spline curve is fitted through the points and returned as a Wire.
    """
    xs = np.linspace(x_range[0], x_range[1], n)
    points = [
        gp_Pnt(float(x), float(y0), float(Ggaussian(float(x), float(y0)))) for x in xs
    ]
    if not smooth:
        poly = BRepBuilderAPI_MakePolygon()
        for p in points:
            poly.Add(p)
        return poly.Wire()

    # build a TColgp_Array1OfPnt for the BSpline fitter (1-based indexing)
    arr = TColgp_Array1OfPnt(1, len(points))
    for i, p in enumerate(points, start=1):
        arr.SetValue(i, p)
    bspline = GeomAPI_PointsToBSpline(arr).Curve()
    edge = BRepBuilderAPI_MakeEdge(bspline).Edge()
    wire = BRepBuilderAPI_MakeWire(edge).Wire()
    return wire


def make_slice_wire_yz(x0, y_range, n, smooth=False):
    """Return a Wire following (x0, y, Ggaussian(x0,y)) sampled along y.

    If smooth=True a B-spline curve is fitted through the points and returned as a Wire.
    """
    ys = np.linspace(y_range[0], y_range[1], n)
    points = [
        gp_Pnt(float(x0), float(y), float(Ggaussian(float(x0), float(y)))) for y in ys
    ]
    if not smooth:
        poly = BRepBuilderAPI_MakePolygon()
        for p in points:
            poly.Add(p)
        return poly.Wire()

    arr = TColgp_Array1OfPnt(1, len(points))
    for i, p in enumerate(points, start=1):
        arr.SetValue(i, p)
    bspline = GeomAPI_PointsToBSpline(arr).Curve()
    edge = BRepBuilderAPI_MakeEdge(bspline).Edge()
    wire = BRepBuilderAPI_MakeWire(edge).Wire()
    return wire


def make_isoline_wire_z(z0, n=360, A=1.0, sx=1.0, sy=1.0, center=(0.0, 0.0)):
    """Create a closed Wire for the isoline Ggaussian(x,y) == z0 (plane z=z0).

    For a Gaussian A*exp(-((x/sx)**2 + (y/sy)**2)) the level set is an ellipse when
    0 < z0 <= A: (x-cx)^2/(rx^2) + (y-cy)^2/(ry^2) = 1 where
      rx = sx * sqrt(-ln(z0/A)), ry = sy * sqrt(-ln(z0/A)).

    Returns a closed Wire (polygon approximating the ellipse).
    """
    if z0 <= 0.0 or z0 > A:
        raise ValueError(f"z0 must satisfy 0 < z0 <= A (A={A}). Got z0={z0}")

    # reuse isoline_points_z to generate gp_Pnt vertices
    pts = isoline_points_z(z0, n=n, A=A, sx=sx, sy=sy, center=center)
    poly = BRepBuilderAPI_MakePolygon()
    for p in pts:
        poly.Add(p)
    poly.Close()
    return poly.Wire()


def isoline_points_z(z0, n=360, A=1.0, wxy=(1.0, 1.0), sxy=(0.0, 0.0), alpha=1.0):
    """Return a list of gp_Pnt points approximating the isoline Ggaussian(x,y)==z0.

    Points are ordered around the closed contour; they are gp_Pnt instances as requested.
    """
    if z0 <= 0.0 or z0 > A:
        raise ValueError(f"z0 must satisfy 0 < z0 <= A (A={A}). Got z0={z0}")

    sx, sy = sxy
    wx, wy = wxy
    # solve A * exp(-alpha * ((x-sx)^2/wx^2 + (y-sy)^2/wy^2)) = z0
    # => (x-sx)^2/(wx^2) + (y-sy)^2/(wy^2) = -ln(z0/A)/alpha
    r = math.sqrt(-math.log(z0 / A) / alpha)
    rx = wx * r
    ry = wy * r

    pts = TColgp_HArray1OfPnt(1, n)
    for i in range(n):
        t = 2.0 * math.pi * i / (n + 1)
        x = sx + rx * math.cos(t)
        y = sy + ry * math.sin(t)
        pts.SetValue(i + 1, gp_Pnt(float(x), float(y), float(z0)))
    return pts


def make_multiple_isolines_z(
    z_min,
    z_max,
    count=10,
    n_per_contour=180,
    A=1.0,
    wxy=(1.0, 1.0),
    sxy=(0.0, 0.0),
    alpha=1.0,
):
    """Create multiple isoline Wires between z_min and z_max (inclusive).

    Returns a tuple (wires, points_lists) where points_lists[i] is the list of gp_Pnt
    that were used to construct wires[i].
    """
    if count <= 0:
        raise ValueError("count must be >= 1")
    # ensure z_min/z_max are in valid range (0, A]
    eps = 1e-12
    z_min_clamped = max(z_min, eps)
    z_max_clamped = min(z_max, A)
    if z_min_clamped >= z_max_clamped:
        raise ValueError(f"z_min ({z_min}) must be < z_max ({z_max}) and within (0,A]")

    zs = np.linspace(z_min_clamped, z_max_clamped, count)
    wires = []
    for z0 in zs:
        pts = isoline_points_z(z0, n=n_per_contour, A=A, wxy=wxy, sxy=sxy, alpha=alpha)
        api = GeomAPI_Interpolate(pts, True, 0.1e-3)
        api.Perform()
        wires.append(make_wire(make_edge(api.Curve())))

    return wires


def make_surface_thru_sections(wirep, wires, peak_vertex):
    """Build a surface (shell/face) passing through the given wires using
    BRepOffsetAPI_ThruSections and return the resulting shape.
    """

    # surf1: sections through wires + (optional) peak vertex
    builder1 = BRepOffsetAPI_ThruSections(False, False, 1e-6)
    for w in wires:
        builder1.AddWire(w)
    builder1.AddVertex(peak_vertex)
    builder1.Build()
    surf1 = builder1.Shape()

    # surf2: sections between wirep and the first wire
    builder2 = BRepOffsetAPI_ThruSections(False, False, 1e-6)
    builder2.AddWire(wirep)
    builder2.AddWire(wires[0])
    builder2.Build()
    surf2 = builder2.Shape()

    # surf3: a face made directly from wirep
    surf3 = make_face(wirep)

    # assemble into a compound and return
    comp = TopoDS_Compound()
    bldr = BRep_Builder()
    bldr.MakeCompound(comp)
    bldr.Add(comp, surf1)
    bldr.Add(comp, surf2)
    return comp, surf3


def make_vertex_shapes_from_points(points_lists):
    """Given an iterable of lists of gp_Pnt, return a flat list of TopoDS_Vertex shapes."""
    verts = []
    for pts in points_lists:
        for p in pts:
            v = BRepBuilderAPI_MakeVertex(p).Vertex()
            verts.append(v)
    return verts


def extract_polygon_2d_from_wire(wire):
    """Extract an ordered list of (x,y) points from a wire by iterating its vertices.

    Note: this assumes the wire is planar and represented by a polygonal wire (vertices at corners).
    """
    pts = []
    exp = TopExp_Explorer(wire, TopAbs_VERTEX)
    while exp.More():
        v = exp.Current()
        p = BRep_Tool.Pnt(v)
        pts.append((p.X(), p.Y()))
        exp.Next()
    # If the wire is closed, ensure the polygon is in order and first point != last
    if pts and pts[0] == pts[-1]:
        pts = pts[:-1]
    return pts


def point_in_polygon(x, y, poly):
    """Return True if point (x,y) is inside polygon poly (list of (x,y)).
    Ray-casting algorithm (even-odd rule).
    """
    inside = False
    n = len(poly)
    if n < 3:
        return False
    j = n - 1
    for i in range(n):
        xi, yi = poly[i]
        xj, yj = poly[j]
        intersect = ((yi > y) != (yj > y)) and (
            x < (xj - xi) * (y - yi) / (yj - yi + 1e-30) + xi
        )
        if intersect:
            inside = not inside
        j = i
    return inside


def sample_points_inside_wire(wire, nx=40, ny=40, zfunc=None, pad=0.0):
    """Sample a regular grid inside the 2D polygon defined by wire.

    - wire: TopoDS_Wire defining the planar polygon region
    - nx, ny: grid resolution in x and y
    - zfunc: optional callable zfunc(x,y) -> z; if None, z=0 is used
    - pad: padding fraction inside bbox (0..0.5)

    Returns a list of gp_Pnt
    """
    poly2d = extract_polygon_2d_from_wire(wire)
    if not poly2d:
        return []
    xs = [p[0] for p in poly2d]
    ys = [p[1] for p in poly2d]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    dx = xmax - xmin
    dy = ymax - ymin
    xmin += dx * pad
    xmax -= dx * pad
    ymin += dy * pad
    ymax -= dy * pad

    pts = []
    if nx <= 0 or ny <= 0:
        return pts
    for i in range(nx):
        x = xmin + (i + 0.5) * (xmax - xmin) / nx
        for j in range(ny):
            y = ymin + (j + 0.5) * (ymax - ymin) / ny
            if point_in_polygon(x, y, poly2d):
                z = zfunc(x, y) if zfunc is not None else 0.0
                pts.append(gp_Pnt(float(x), float(y), float(z)))
    return pts


def make_random_ax3_at_point(p):
    """Create a gp_Ax3 at gp_Pnt p with a random gp_Dir as Z direction."""

    # random direction
    rx, ry, rz = random.uniform(-1, 1), random.uniform(-1, 1), random.uniform(-1, 1)
    # avoid zero vector
    if abs(rx) < 1e-6 and abs(ry) < 1e-6 and abs(rz) < 1e-6:
        rz = 1.0
    zdir = gp_Dir(rx, ry, rz)
    # ensure Z component is positive
    if zdir.Z() < 0:
        zdir = gp_Dir(-zdir.X(), -zdir.Y(), -zdir.Z())

    # pick a helper axis to construct an orthonormal X direction
    helper = gp_Vec(1, 0, 0)
    zv = gp_Vec(zdir.X(), zdir.Y(), zdir.Z())
    xvec = helper.Crossed(zv)
    if xvec.Magnitude() < 1e-6:
        helper = gp_Vec(0, 1, 0)
        xvec = helper.Crossed(zv)
    xdir = gp_Dir(xvec)

    return gp_Ax3(p, zdir, xdir)


def sample_n_points_inside_wire(wire, n=100, max_attempts=5000, with_axes=False):
    """Sample exactly n points inside the planar polygon defined by wire.

    - wire: TopoDS_Wire (planar)
    - n: desired number of points (default 100)
    - max_attempts: maximum random trials before falling back to grid fill

    Returns a list of gp_Pnt. The z coordinate for all points is taken as the
    average Z of the wire's vertices (i.e. lying on the wire's plane).
    """
    # collect 3D vertex coordinates from the wire
    verts = []
    exp = TopExp_Explorer(wire, TopAbs_VERTEX)
    while exp.More():
        v = exp.Current()
        p = BRep_Tool.Pnt(v)
        verts.append((float(p.X()), float(p.Y()), float(p.Z())))
        exp.Next()
    if not verts:
        return []

    # 2D polygon (XY) and plane z
    poly2d = [(x, y) for (x, y, z) in verts]
    if poly2d and poly2d[0] == poly2d[-1]:
        poly2d = poly2d[:-1]
    xs = [p[0] for p in poly2d]
    ys = [p[1] for p in poly2d]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    z_plane = sum(z for (_x, _y, z) in verts) / len(verts)

    pts = []
    axes = []
    attempts = 0
    while len(pts) < n and attempts < max_attempts:
        attempts += 1
        x = random.uniform(xmin, xmax)
        y = random.uniform(ymin, ymax)
        if point_in_polygon(x, y, poly2d):
            p = gp_Pnt(float(x), float(y), float(z_plane))
            pts.append(p)
            if with_axes:
                axes.append(make_random_ax3_at_point(p))

    # If random sampling didn't find enough points, fill using a grid
    if len(pts) < n:
        need = n - len(pts)
        # choose grid size to likely provide enough candidates
        g = int(math.ceil(math.sqrt(need * 4)))
        for i in range(g):
            if len(pts) >= n:
                break
            for j in range(g):
                if len(pts) >= n:
                    break
                x = xmin + (i + 0.5) * (xmax - xmin) / g
                y = ymin + (j + 0.5) * (ymax - ymin) / g
                if point_in_polygon(x, y, poly2d):
                    p = gp_Pnt(float(x), float(y), float(z_plane))
                    pts.append(p)
                    if with_axes:
                        axes.append(make_random_ax3_at_point(p))
                    if len(pts) >= n:
                        break

    if with_axes:
        return pts[:n], axes[:n]
    return pts[:n]


def sample_n_points_on_face(face, n=100, max_attempts=5000, tol=1e-9):
    """Sample exactly n points inside the (possibly trimmed) face using UV-space sampling.

    - face: TopoDS_Face
    - n: number of points to return
    - max_attempts: how many random UV trials before falling back to UV-grid
    """
    try:
        umin, umax, vmin, vmax = breptools.UVBounds(face)
    except Exception:
        return []

    classifier = BRepTopAdaptor_FClass2d(face, tol)
    srf_adapt = BRepAdaptor_Surface(face)

    axes = []
    attempts = 0
    while len(axes) < n and attempts < max_attempts:
        attempts += 1
        u = random.uniform(umin, umax)
        v = random.uniform(vmin, vmax)
        uv = gp_Pnt2d(u, v)
        if classifier.Perform(uv) == TopAbs_IN:
            p = srf_adapt.Value(u, v)
            gp = gp_Pnt(float(p.X()), float(p.Y()), float(p.Z()))
            axes.append(make_random_ax3_at_point(gp))

    # fallback to uniform UV-grid sampling
    if len(axes) < n:
        need = n - len(axes)
        g = int(math.ceil(math.sqrt(need * 4)))
        for i in range(g):
            if len(axes) >= n:
                break
            for j in range(g):
                if len(axes) >= n:
                    break
                u = umin + (i + 0.5) * (umax - umin) / g
                v = vmin + (j + 0.5) * (vmax - vmin) / g
                uv = gp_Pnt2d(u, v)
                if classifier.Perform(uv) == TopAbs_IN:
                    p = srf_adapt.Value(u, v)
                    gp = gp_Pnt(float(p.X()), float(p.Y()), float(p.Z()))
                    axes.append(make_random_ax3_at_point(gp))
                    if len(axes) >= n:
                        break

    return axes[:n]


def sample_n_points_on_face_uniform(
    face,
    n=100,
    grid_u=60,
    grid_v=60,
    tol=1e-9,
    density_func=None,
    use_area_weight=True,
    region_func=None,
    unique=False,
):
    """Sample n points on a (possibly trimmed) face by (area * spatial density).

    Parameters:
    - face: TopoDS_Face
    - n: desired number of sampled points
    - grid_u, grid_v: UV grid resolution used to discretize the surface for weighting
    - tol: tolerance for classifiers / surface props
    - density_func: optional callable f(x,y)->weight (spatial density). If None, treated as 1.
    - use_area_weight: when True multiply spatial density by |D1U x D1V|; when False use only density_func
    - region_func: optional callable g(x,y)->bool. If provided only points where g(x,y) is True are accepted.
    - unique: if True, attempt to avoid near-duplicate points by rounding coordinates.

    Approach:
    - Build a UV grid and evaluate per-cell weight = (area density) * density_func(x,y) or density_func alone.
    - Repeatedly select cells by weight and sample uniformly inside the chosen cell. Reject points that lie
      outside the trimmed face or that fail region_func. Retry until n points are collected or a maximum
      attempt limit is reached.
    """
    try:
        umin, umax, vmin, vmax = breptools.UVBounds(face)
    except Exception:
        return []

    srf_adapt = BRepAdaptor_Surface(face)
    classifier = BRepTopAdaptor_FClass2d(face, tol)

    # grid edges and cell centers
    us = np.linspace(umin, umax, grid_u + 1)
    vs = np.linspace(vmin, vmax, grid_v + 1)
    u_centers = 0.5 * (us[:-1] + us[1:])
    v_centers = 0.5 * (vs[:-1] + vs[1:])

    weights = np.zeros((len(u_centers), len(v_centers)), dtype=float)
    # evaluate area density at each cell center
    for i, u in enumerate(u_centers):
        for j, v in enumerate(v_centers):
            try:
                prop = BRepLProp_SLProps(srf_adapt, u, v, 1, tol)
                du = prop.D1U()
                dv = prop.D1V()
                area_dens = du.Crossed(dv).Magnitude()
                # evaluate surface point to apply optional spatial density weighting
                p = srf_adapt.Value(u, v)
                x, y = float(p.X()), float(p.Y())
                spatial_w = 1.0
                if density_func is not None:
                    try:
                        spatial_w = float(density_func(x, y))
                    except Exception:
                        spatial_w = 1.0
                if use_area_weight:
                    weights[i, j] = max(area_dens * max(spatial_w, 0.0), 0.0)
                else:
                    weights[i, j] = max(max(spatial_w, 0.0), 0.0)
            except Exception:
                weights[i, j] = 0.0

    total = weights.sum()
    if total <= 0.0:
        # fallback to existing random uv sampling
        return sample_n_points_on_face(face, n=n, max_attempts=5000, tol=tol)

    flat_weights = weights.ravel() / total

    axes = []
    attempts = 0
    max_try = max(1000, n * 20)
    # sample until we have n valid points or we exhaust attempts
    while len(axes) < n and attempts < max_try:
        attempts += 1
        k = np.random.choice(flat_weights.size, p=flat_weights)
        i = k // weights.shape[1]
        j = k % weights.shape[1]
        u0, u1 = us[i], us[i + 1]
        v0, v1 = vs[j], vs[j + 1]
        # sample uniformly within the chosen cell
        u = float(np.random.uniform(u0, u1))
        v = float(np.random.uniform(v0, v1))
        uv = gp_Pnt2d(u, v)
        # ensure point is inside trimmed face
        if classifier.Perform(uv) != TopAbs_IN:
            # try a few times locally
            ok = False
            for _ in range(5):
                u = float(np.random.uniform(u0, u1))
                v = float(np.random.uniform(v0, v1))
                uv = gp_Pnt2d(u, v)
                if classifier.Perform(uv) == TopAbs_IN:
                    ok = True
                    break
            if not ok:
                continue
        p = srf_adapt.Value(u, v)
        gp_pt = gp_Pnt(float(p.X()), float(p.Y()), float(p.Z()))
        axes.append(make_random_ax3_at_point(gp_pt))

    # if we sampled fewer due to trimmed cells, fill the rest with fallback
    if len(axes) < n:
        need = n - len(axes)
        fallback = sample_n_points_on_face(face, n=need, max_attempts=5000, tol=tol)
        axes.extend(fallback)

    return axes[:n]


def find_first_intersection(shape, p0, dir_vec, tol=1.0e-9):
    """Find the nearest intersection of a ray (p0, dir_vec) with shape.

    Returns tuple (p_hit, face, u, v, w) or None if no hit.
    """
    try:
        ddir = gp_Dir(dir_vec.X(), dir_vec.Y(), dir_vec.Z())
    except Exception:
        # fallback: dir_vec may be tuple-like
        ddir = gp_Dir(float(dir_vec[0]), float(dir_vec[1]), float(dir_vec[2]))

    lin = gp_Lin(p0, ddir)
    api = BRepIntCurveSurface_Inter()
    api.Init(shape, lin, tol)
    best = None
    best_dist = float("inf")
    while api.More():
        p_hit = api.Pnt()
        d = p0.Distance(p_hit)
        if d < best_dist and api.W() > 1.0e-9:
            best_dist = d
            best = (p_hit, api.Face(), api.U(), api.V(), api.W())
        api.Next()
    return best


def reflect_beam_on_face(face, incoming_ax3, u, v, eps=1.0e-6):
    """Reflect an incoming gp_Ax3 off the given face at param (u,v).

    - face: TopoDS_Face where the hit occurred
    - p_hit: gp_Pnt hit location (may be ignored in favor of surface-evaluated point)
    - incoming_ax3: gp_Ax3 representing the incoming ray (location + direction)
    - u, v: surface parameters at hit (used to compute normal)
    - eps: small offset to move the returned origin along the reflected direction

    Returns a new gp_Ax3 for the reflected beam, or None on failure.
    """
    try:
        surf = BRepAdaptor_Surface(face)
        prop = BRepLProp_SLProps(surf, u, v, 2, 1.0e-9)
        p_surf = prop.Value()
        du = prop.D1U()
        dv = prop.D1V()
    except Exception:
        return None

    # compute surface normal (du x dv)
    vz = du.Crossed(dv)
    if vz.Magnitude() < 1e-12:
        return None

    # ensure normal points against the incoming direction
    inc_dir = incoming_ax3.Direction()
    inc_vec = gp_Vec(inc_dir.X(), inc_dir.Y(), inc_dir.Z())
    if vz.Dot(inc_vec) > 0:
        vz.Reverse()

    du.Normalize()
    dv.Normalize()
    vz.Normalize()

    # construct a local axis describing the surface normal at the hit

    norm_ax3 = gp_Ax3(
        p_surf, gp_Dir(vz.X(), vz.Y(), vz.Z()), gp_Dir(du.X(), du.Y(), du.Z())
    )

    # build an axis for the incoming beam positioned at the surface point
    in_ax3 = gp_Ax3(p_surf, incoming_ax3.Direction(), incoming_ax3.XDirection())

    # reflect the incoming axis about the surface normal plane
    try:
        in_ax3.Mirror(norm_ax3.Ax2())
        # ensure the reflected direction points away from the surface normal
        if in_ax3.Direction().Dot(norm_ax3.Direction()) < 0:
            in_ax3.ZReverse()
    except Exception:
        return None

    # offset the origin slightly along the reflected direction to avoid self-intersection
    rd = in_ax3.Direction()
    new_origin = gp_Pnt(
        p_surf.X() + rd.X() * eps, p_surf.Y() + rd.Y() * eps, p_surf.Z() + rd.Z() * eps
    )

    return gp_Ax3(new_origin, in_ax3.Direction(), in_ax3.XDirection())


def simulate_reflected_returns(
    plan_face,
    surf_shape,
    start_axes=None,
    max_bounces=5,
    require_reflections=None,
    require_at_least=False,
    total_reflections=None,
    stop_on_first_plan_hit=True,
):
    """Shoot rays from points on plan_face toward surf_shape, reflect, and
    collect points on plan_face where rays return after >=1 reflections.
    """
    # use provided start points/axes if available, otherwise sample on the plan face
    results = []
    # trajectories: list of lists of gp_Pnt describing each ray path (origin, hit1, hit2, ...)
    trajectories = []
    # per-ray statistics aligned with start_axes indices
    stats = []

    # iterate over provided axes (each is a gp_Ax3)
    for i, axs in enumerate(start_axes):
        # if an axis/direction was provided for this start point, use it
        dir_vec = gp_Vec(axs.Direction())
        if dir_vec.Magnitude() < 1e-9:
            continue
        dir_vec.Normalize()
        origin = axs.Location()
        bounced = 0
        traj_pts = [axs.Location()]

        returned_here = False
        plan_hits = 0
        # Unified attempt-based loop (option A): perform reflections until the beam
        # hits the plan (stop on first plan hit) or until a safety cap is reached.
        # If total_reflections (threshold) is provided, it acts as a soft threshold
        # (we still stop on the first plan hit whenever it occurs). The safety cap
        # prevents infinite loops in pathological cases.
        threshold = int(total_reflections) if total_reflections is not None else None
        safety_cap = max(200, int(max_bounces) * 10)
        # hard cap on loop iterations (protects against infinite loops)
        iterations = 0

        while iterations < safety_cap and bounced < max_bounces * 100:
            iterations += 1
            # check both surface and plan for nearest intersection from origin
            hit_surf = find_first_intersection(surf_shape, origin, dir_vec)
            hit_plan = find_first_intersection(plan_face, origin, dir_vec)

            # if neither hit, ray escapes
            if hit_surf is None and hit_plan is None:
                # no intersection with either shape -> ray escapes
                break

            # helper to get distance (or inf if missing)
            def _dist(hit):
                return origin.Distance(hit[0]) if hit is not None else float("inf")

            d_surf = _dist(hit_surf)
            d_plan = _dist(hit_plan)

            # choose nearest hit (tie: prefer plan)
            hit_first = None
            hit_type = None
            if d_plan <= d_surf:
                hit_first = hit_plan
                hit_type = "plan"
            else:
                hit_first = hit_surf
                hit_type = "surf"

            if hit_type == "surf":
                p_hit, face_hit, u, v, w = hit_first
                # record surf hit
                traj_pts.append(p_hit)

                # reflect on surface
                reflected_beam = reflect_beam_on_face(
                    face_hit,
                    gp_Ax3(origin, gp_Dir(dir_vec.X(), dir_vec.Y(), dir_vec.Z())),
                    u,
                    v,
                )
                if reflected_beam is None:
                    # reflection computation failed; stop this ray
                    break

                origin = reflected_beam.Location()
                dir_vec = gp_Vec(reflected_beam.Direction())
                bounced += 1
                # continue tracing from reflected ray
                continue

            else:
                # hit_type == 'plan'
                p_plan, f_plan, uu, vv, ww = hit_first
                traj_pts.append(p_plan)
                plan_hits += 1

            meets_req = True
            if require_reflections is not None:
                if require_at_least:
                    meets_req = bounced >= int(require_reflections)
                else:
                    meets_req = bounced == int(require_reflections)
            if meets_req:
                results.append(p_plan)
            returned_here = True

            # By default, stop tracing this ray after the first plan hit to avoid
            # multiple returns being recorded per start axis. If caller explicitly
            # sets stop_on_first_plan_hit=False, reflect and continue.
            if stop_on_first_plan_hit:
                break

            # reflect off the plan face and continue tracing if possible
            reflected_plan_beam = reflect_beam_on_face(
                f_plan,
                gp_Ax3(origin, gp_Dir(dir_vec.X(), dir_vec.Y(), dir_vec.Z())),
                uu,
                vv,
            )
            if reflected_plan_beam is None:
                break

            origin = reflected_plan_beam.Location()
            dir_vec = gp_Vec(reflected_plan_beam.Direction())
            bounced += 1
            # continue tracing
            continue

        # store trajectory for this ray (may be only origin if no hits)
        trajectories.append(traj_pts)

        # compute statistics for this ray efficiently using returned_here flag
        ret_flag = returned_here
        reflections = None
        path_length = None
        start_to_return = None
        if ret_flag:
            reflections = bounced
            # trajectory path length
            plen = 0.0
            for a in range(len(traj_pts) - 1):
                plen += traj_pts[a].Distance(traj_pts[a + 1])
            path_length = plen
            # compute start-to-return using the start axis location
            try:
                start_pnt = axs.Location()
                last_pt = traj_pts[-1]
                start_to_return = start_pnt.Distance(last_pt)
            except Exception:
                start_to_return = None

        stats.append(
            {
                "index": i,
                "returned": ret_flag,
                "reflections": reflections,
                "path_length": path_length,
                "start_to_return": start_to_return,
                "plan_hits": plan_hits,
                "meets_requirement": (
                    ret_flag
                    and (
                        require_reflections is None
                        or (
                            require_at_least
                            and reflections is not None
                            and reflections >= int(require_reflections)
                        )
                        or (
                            (not require_at_least)
                            and reflections is not None
                            and reflections == int(require_reflections)
                        )
                    )
                ),
            }
        )

    return results, trajectories, stats


if __name__ == "__main__":
    # Create 20 isolines between z_min and z_max (parallel to XY plane).
    z_min = 0.1
    z_max = 1.99
    count = 20
    n_per_contour = 200
    A = 2.0

    sxy = (0.0, 0.0)
    wxy = (2.0, 2.5)
    wires = make_multiple_isolines_z(
        z_min,
        z_max,
        count=count,
        n_per_contour=n_per_contour,
        A=A,
        wxy=wxy,
        sxy=sxy,
    )

    axs = gp_Ax3()
    pln = make_face(gp_Pln(axs), -10, 10, -10, 10)
    proj = BRepProj_Projection(wires[0], pln, axs.Direction())
    wirep = proj.Current()

    display, start_display, add_menu, add_function_to_menu = init_display()

    # display Gaussian peak vertex at the center (z = A)
    cx, cy = sxy
    peak_pnt = gp_Pnt(float(cx), float(cy), float(A))
    peak_v = BRepBuilderAPI_MakeVertex(peak_pnt).Vertex()

    # Build a surface through the wires using ThruSections
    surf, plan = make_surface_thru_sections(wirep, wires, peak_v)
    display.DisplayShape(surf, transparency=0.7)
    display.DisplayShape(peak_v)
    display.DisplayShape(wires[0])
    display.DisplayShape(wires[1])
    display.DisplayShape(plan, color="RED", transparency=0.5)

    # Build a Face from the projected outer wire and sample directly on the Face
    # multiple density functions and region constraints (pick via selectors)
    sx, sy = wxy
    cx, cy = sxy

    def density_uniform(x, y):
        return 1.0

    def density_gaussian(x, y, A=A, sx=sx, sy=sy, cx=cx, cy=cy, alpha=1.0):
        return float(Ggaussian(x, y, A=A, sx=sx, sy=sy, cx=cx, cy=cy, alpha=alpha))

    def density_ring(x, y, r0=0.6, sigma=0.12, A=1.0):
        r = math.hypot(x - cx, y - cy)
        return float(A * math.exp(-((r - r0) ** 2) / (2.0 * sigma**2)))

    def density_bimodal(x, y, dx=0.6, dy=0.0, A=1.0, sigma=0.25):
        r1 = math.hypot(x - (cx - dx), y - (cy - dy))
        r2 = math.hypot(x - (cx + dx), y - (cy + dy))
        return float(
            A
            * (math.exp(-0.5 * (r1 / sigma) ** 2) + math.exp(-0.5 * (r2 / sigma) ** 2))
        )

    # region constraints
    def region_circle(x, y, radius=1.2):
        return (x - cx) ** 2 + (y - cy) ** 2 <= radius * radius

    def region_rect(x, y, xmin=-0.25, xmax=0.25, ymin=-0.5, ymax=0.5):
        return xmin <= x <= xmax and ymin <= y <= ymax

    def region_donut(x, y, r_in=0.3, r_out=1.2):
        r2 = (x - cx) ** 2 + (y - cy) ** 2
        return r_in * r_in <= r2 <= r_out * r_out

    def region_paraboloid(x, y, z_thresh=0.2):
        # use the Gaussian height as simple bowl height and threshold
        z = float(Ggaussian(x, y, A=A, sx=sx, sy=sy, cx=cx, cy=cy, alpha=1.0))
        return z >= z_thresh

    density_options = {
        "uniform": density_uniform,
        "gaussian": density_gaussian,
        "ring": density_ring,
        "bimodal": density_bimodal,
    }

    region_options = {
        "circle": region_circle,
        "rect": region_rect,
        "donut": region_donut,
        "paraboloid": region_paraboloid,
    }

    # Select which density and region to use (change these strings to try others)
    chosen_density = "gaussian"
    chosen_region = "rect"

    density_func = density_options[chosen_density]
    region_func = region_options[chosen_region]

    print("Density options:", ",".join(density_options.keys()))
    print("Region options:", ",".join(region_options.keys()))
    print(f"Using density={chosen_density}, region={chosen_region}")

    sampled_axes = sample_n_points_on_face_uniform(
        plan,
        n=1000,
        grid_u=80,
        grid_v=40,
        density_func=density_func,
        use_area_weight=False,
        region_func=region_func,
    )
    print(f"Sampled {len(sampled_axes)} axes on plan")

    # --- Option blocks (1..4) -------------------------------------------------
    # Prepare four helper features below. Only option 1 is enabled by default.

    # 1) Display small spheres at each sampled start point (enabled)
    # from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeSphere
    #
    # sphere_radius = 0.02
    # sample_spheres = []
    # for ax in sampled_axes:
    #    p = ax.Location()
    #    try:
    #        s = BRepPrimAPI_MakeSphere(p, sphere_radius).Shape()
    #        sample_spheres.append(s)
    #        display.DisplayShape(s, color="BLUE1", transparency=0.0)
    #    except Exception:
    #        pass
    #
    # 2) Create a bowl geometry (paraboloid-like) and display it (commented)
    # from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeCone
    # bowl = BRepPrimAPI_MakeCone(gp_Pnt(0,0,-0.5), 1.5, 0.05, 0.8).Shape()
    # display.DisplayShape(bowl, color="BROWN", transparency=0.6)

    # 3) Create a matplotlib contour/mask from the Gaussian density and show on plots (commented)
    # xs = np.linspace(-1.5, 1.5, 120)
    # ys = np.linspace(-1.0, 1.0, 120)
    # Xg, Yg = np.meshgrid(xs, ys)
    # Zg = np.vectorize(lambda x,y: _gauss_density(x,y))(Xg, Yg)
    # # will be displayed later in matplotlib axes as contourf

    # 4) Unique-point filter visualization (commented)
    # def unique_by_round(pts, decimals=3):
    #     seen = set()
    #     out = []
    #     for p in pts:
    #         key = (round(p.X(), decimals), round(p.Y(), decimals))
    #         if key in seen:
    #             continue
    #         seen.add(key)
    #         out.append(p)
    #     return out

    # ---------------------------------------------------------------------------

    # Display axis directions as short edges (Z direction of each gp_Ax3)
    # determine a reasonable scale from sampled points bbox
    scale = 0.2

    # Simulate reflected returns: shoot beams toward surf and collect points that return to plan
    returns, trajectories, stats = simulate_reflected_returns(
        plan, surf, start_axes=sampled_axes, max_bounces=50
    )
    print(f"Reflected returns to plan: {len(returns)}")
    ret_verts = [BRepBuilderAPI_MakeVertex(p).Vertex() for p in returns]
    for v in ret_verts:
        display.DisplayShape(v, color="GREEN")

    # draw trajectories as short orange edges between successive points
    for traj in trajectories:
        if len(traj) < 2:
            continue
        for i in range(len(traj) - 1):
            e = BRepBuilderAPI_MakeEdge(traj[i], traj[i + 1]).Edge()
            display.DisplayShape(e, color="YELLOW")

    display.FitAll()

    # --- plot statistics with matplotlib ---
    # filter stats for returned rays
    def get_stats(key):
        return [s[key] for s in stats if s.get(key) and s.get(key) is not None]

    returned_stats = [s for s in stats if s.get("returned")]
    ref_counts = get_stats("reflections")
    path_lengths = get_stats("path_length")
    start_to_returns = get_stats("start_to_return")
    ret_xy = [(p.X(), p.Y()) for p in returns]
    plan_hits_all = [s.get("plan_hits", 0) for s in stats]
    plan_hits_nonzero = [p for p in plan_hits_all if p is not None]
    path_lengths_all = [s.get("path_length") for s in stats]

    # collect start locations for returned rays using stats.index
    start_xy = []
    for s in returned_stats:
        idx = s.get("index")
        try:
            sp = sampled_axes[idx].Location()
            start_xy.append((sp.X(), sp.Y()))
        except Exception:
            pass

    plt.figure(figsize=(15, 4))
    plt.suptitle("Returned hit locations (plan) and start locations")

    plt.subplot(1, 4, 1)
    if ref_counts:
        plt.hist(
            ref_counts, bins=range(0, max(ref_counts) + 2), align="left", color="C0"
        )
    else:
        # no returned rays: display empty plot with message
        plt.text(0.5, 0.5, "No returned rays", ha="center", va="center")
    plt.xlabel("Reflections until return")
    plt.ylabel("Count")

    plt.subplot(1, 4, 2)
    if path_lengths:
        plt.hist(path_lengths, bins=20, color="C1")
    else:
        plt.text(0.5, 0.5, "No path length data", ha="center", va="center")
    plt.xlabel("Trajectory path length")

    def plot_hist2d(axs, rxy, bins=(40, 40), range_xy=None, vmin=None, vmax=None):
        """Plot a 2D histogram as a contourf on `axs`.

        Returns the QuadContourSet (or None if rxy is empty) and the used range as
        (xmin, xmax, ymin, ymax).
        """
        if not rxy:
            return None, None
        rx, ry = zip(*rxy)
        # determine plotting range
        if range_xy is None:
            xmin = min(rx)
            xmax = max(rx)
            ymin = min(ry)
            ymax = max(ry)
            # add a tiny margin so points on the edge are visible
            dx = (xmax - xmin) * 0.02 if xmax > xmin else 0.1
            dy = (ymax - ymin) * 0.02 if ymax > ymin else 0.1
            xmin -= dx
            xmax += dx
            ymin -= dy
            ymax += dy
        else:
            xmin, xmax, ymin, ymax = range_xy

        # compute 2D histogram density with the specified bins and range
        H, xedges, yedges = np.histogram2d(
            rx, ry, bins=bins, range=[[xmin, xmax], [ymin, ymax]]
        )
        # meshgrid of bin centers
        xcenters = (xedges[:-1] + xedges[1:]) / 2.0
        ycenters = (yedges[:-1] + yedges[1:]) / 2.0
        Xc, Yc = np.meshgrid(xcenters, ycenters)
        # draw contourf; norm can be passed to share a color scale across axes
        if vmin is None or vmax is None:
            cf = axs.contourf(Xc, Yc, H.T, cmap="jet", alpha=0.7)
        else:
            norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
            cf = axs.contourf(Xc, Yc, H.T, cmap="jet", alpha=0.7, norm=norm)
        # overlay returned points for clarity
        axs.scatter(rx, ry, c="red", s=10, edgecolors="k", linewidths=0.2)
        return cf, (xmin, xmax, ymin, ymax)

    ax3 = plt.subplot(1, 4, 3)
    ax4 = plt.subplot(1, 4, 4)

    # Determine a common plotting range and vmax so both histograms share a colorbar
    bins = (40, 40)
    combined_x = []
    combined_y = []
    if start_xy:
        sx, sy = zip(*start_xy)
        combined_x.extend(sx)
        combined_y.extend(sy)
    if ret_xy:
        rx_vals, ry_vals = zip(*ret_xy)
        combined_x.extend(rx_vals)
        combined_y.extend(ry_vals)

    if combined_x and combined_y:
        xmin = min(combined_x)
        xmax = max(combined_x)
        ymin = min(combined_y)
        ymax = max(combined_y)
        # small margin
        dx = (xmax - xmin) * 0.02 if xmax > xmin else 0.1
        dy = (ymax - ymin) * 0.02 if ymax > ymin else 0.1
        plot_range = (xmin - dx, xmax + dx, ymin - dy, ymax + dy)

        # compute histograms to determine a common vmax
        H_start = None
        H_ret = None
        if start_xy:
            sx, sy = zip(*start_xy)
            H_start, _, _ = np.histogram2d(
                sx,
                sy,
                bins=bins,
                range=[[plot_range[0], plot_range[1]], [plot_range[2], plot_range[3]]],
            )
        if ret_xy:
            rx_vals, ry_vals = zip(*ret_xy)
            H_ret, _, _ = np.histogram2d(
                rx_vals,
                ry_vals,
                bins=bins,
                range=[[plot_range[0], plot_range[1]], [plot_range[2], plot_range[3]]],
            )

        vmax = 0
        if H_start is not None:
            vmax = max(vmax, float(H_start.max()))
        if H_ret is not None:
            vmax = max(vmax, float(H_ret.max()))

        # create a shared Normalize and ScalarMappable for a common colorbar
        norm = mpl.colors.Normalize(vmin=0, vmax=vmax)
        mappable = mpl.cm.ScalarMappable(norm=norm, cmap="jet")
        # draw contourf using the same norm so both plots share color scaling
        if start_xy:
            plot_hist2d(
                ax3, start_xy, bins=bins, range_xy=plot_range, vmin=0, vmax=vmax
            )
        if ret_xy:
            plot_hist2d(ax4, ret_xy, bins=bins, range_xy=plot_range, vmin=0, vmax=vmax)
        # attach single colorbar using the shared mappable
        fig = plt.gcf()
        # create a colorbar axis to the right of ax4 so it doesn't overlap
        try:
            divider = make_axes_locatable(ax4)
            cax = divider.append_axes("right", size="3%", pad=0.12)
            fig.colorbar(mappable, cax=cax, label="counts")
        except Exception:
            # fallback to default placement
            fig.colorbar(mappable, ax=[ax3, ax4], label="counts")
    else:
        # fallback: plot what exists without a shared colorbar
        if start_xy:
            plot_hist2d(ax3, start_xy)
        if ret_xy:
            plot_hist2d(ax4, ret_xy)

    ax3.set_xlabel("X")
    ax3.set_ylabel("Y")
    ax4.set_xlabel("X")
    ax4.set_ylabel("Y")
    plt.tight_layout()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    # ax1: histogram of plan hit counts per ray
    if plan_hits_nonzero:
        ax1.hist(
            plan_hits_nonzero,
            bins=range(0, max(plan_hits_nonzero) + 2),
            color="C3",
            align="left",
        )
    ax1.set_xlabel("Plan hits per ray")
    ax1.set_ylabel("Count")

    # ax2: scatter plan_hits vs trajectory path length (for rays that returned)
    xs = []
    ys = []
    for s in stats:
        ph = s.get("plan_hits")
        pl = s.get("path_length")
        if ph is not None and pl is not None:
            xs.append(ph)
            ys.append(pl)
    if xs and ys:
        ax2.scatter(xs, ys, c="C4")
    ax2.set_xlabel("Plan hits")
    ax2.set_ylabel("Path length")
    ax2.set_title("Plan hits vs path length")
    plt.tight_layout()
    plt.show()
    start_display()
