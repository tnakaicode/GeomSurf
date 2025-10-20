"""Fit a 3D line to a point cloud and return an OCC Edge (gp_Lin -> TopoDS_Edge).

Provides:
    fit_line_edge(points, extend=0.0, return_line=False)

Arguments:
    points: array-like of shape (N,3) or iterable of OCC gp_Pnt
    extend: float, additional length to extend both ends of the fitted segment (same units as points)
    return_line: if True, also return the gp_Lin object

Dependencies:
    numpy, pythonocc-core (OCC)

This module is a small utility intended for use in scripts and notebooks where
pythonocc is available. The function performs a PCA-based line fit (principal
component) and builds a bounded edge along the principal direction covering the
projection range of the input points, optionally extended.
"""

import numpy as np

from OCC.Core.gp import gp_Pnt, gp_Dir, gp_Lin
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge
from OCC.Core.TopoDS import TopoDS_Edge
from OCC.Core.GeomAPI import GeomAPI_ProjectPointOnCurve
from OCC.Core.BRep import BRep_Tool, BRep_Builder
from OCC.Core.TopoDS import TopoDS_Compound
from OCC.Core.Geom import Geom_Line


def _to_numpy(points):
    """Convert input points to Nx3 numpy array.

    Accepts sequences of (x,y,z) or OCC gp_Pnt objects.
    """
    pts = list(points)
    if len(pts) == 0:
        return np.zeros((0, 3), dtype=float)

    # Assume each item provides X(), Y(), Z() (e.g. OCC gp_Pnt). No type checks.
    arr = np.array([[p.X(), p.Y(), p.Z()] for p in pts], dtype=float)
    return arr


def gp_pnt_to_array(p):
    """Convert a gp_Pnt to a numpy array [x,y,z]."""
    return np.array([p.X(), p.Y(), p.Z()], dtype=float)


def point_to_segment_distance(pa, pb, p):
    """Return distance from point p to segment [pa,pb]. pa,pb,p are numpy arrays."""
    v = pb - pa
    w = p - pa
    c1 = np.dot(w, v)
    if c1 <= 0:
        return np.linalg.norm(p - pa)
    c2 = np.dot(v, v)
    if c2 <= c1:
        return np.linalg.norm(p - pb)
    b = c1 / c2
    pb_proj = pa + b * v
    return np.linalg.norm(p - pb_proj)


def dist_point_to_curve_gp(p, geom_curve, origin=None):
    """Distance from gp_Pnt p to Geom_Curve geom_curve using GeomAPI; fallback to origin if needed."""
    projp = GeomAPI_ProjectPointOnCurve(p, geom_curve)
    if projp.NbPoints() > 0:
        nearest = projp.NearestPoint()
        return p.Distance(nearest)
    if origin is not None:
        arr = gp_pnt_to_array(p)
        return float(np.linalg.norm(arr - origin))
    return float(np.linalg.norm(gp_pnt_to_array(p)))


def dist_point_to_polyline(p, sampled_pts):
    """Distance from gp_Pnt or numpy point to sampled polyline points (numpy array Nx3)."""
    if not hasattr(p, "X"):
        p_arr = np.asarray(p, dtype=float)
    else:
        p_arr = gp_pnt_to_array(p)
    min_d = float("inf")
    for i in range(len(sampled_pts) - 1):
        d = point_to_segment_distance(sampled_pts[i], sampled_pts[i + 1], p_arr)
        if d < min_d:
            min_d = d
    return min_d


def fit_line_edge(points, extend=0.0, poly_degree=1, n_samples=200):
    """Fit a 3D line to points and return a bounded TopoDS_Edge along that line.

    The line direction is obtained from the principal component (PCA). The
    returned edge is a straight edge on the fitted gp_Lin, bounded between the
    minimum and maximum projected scalars along the line. The bounds are
    optionally extended by `extend` on both ends.

    Returns:
      TopoDS_Edge or (TopoDS_Edge, gp_Lin) when return_line is True.
    """
    pts = _to_numpy(points)
    if pts.shape[0] < 2:
        raise ValueError("At least two points are required to fit a line")

    # Compute centroid and principal direction via PCA (covariance)
    centroid = pts.mean(axis=0)
    X = pts - centroid
    cov = np.dot(X.T, X) / float(max(1, pts.shape[0] - 1))
    eigvals, eigvecs = np.linalg.eigh(cov)
    # principal eigenvector corresponds to largest eigenvalue
    principal = eigvecs[:, np.argmax(eigvals)]
    norm = np.linalg.norm(principal)
    if norm == 0:
        raise ValueError("Points are degenerate; cannot determine a direction")
    direction = principal / norm

    # Project points onto the principal axis to get parameter range
    ts = np.dot(X, direction)
    tmin, tmax = float(np.min(ts)), float(np.max(ts))
    # extend both ends by given absolute amount
    if extend and float(extend) > 0.0:
        tmin -= extend
        tmax += extend

    # Build OCC objects
    p_origin = gp_Pnt(float(centroid[0]), float(centroid[1]), float(centroid[2]))
    p_dir = gp_Dir(float(direction[0]), float(direction[1]), float(direction[2]))
    line = gp_Lin(p_origin, p_dir)

    # Create an edge along the line between parameters tmin..tmax
    maker = BRepBuilderAPI_MakeEdge(line, tmin, tmax)
    edge = maker.Edge()

    # Determine if polynomial fitting requested via poly_degree
    poly_coeffs = None
    poly_compound = None
    sampled = None
    if poly_degree is not None:
        # Fit polynomials x(t), y(t), z(t) where t is projection along principal axis
        t_vals = ts
        coords = pts  # pts is already numpy array in this scope
        coeff_x = np.polyfit(t_vals, coords[:, 0], poly_degree)
        coeff_y = np.polyfit(t_vals, coords[:, 1], poly_degree)
        coeff_z = np.polyfit(t_vals, coords[:, 2], poly_degree)
        poly_coeffs = (coeff_x, coeff_y, coeff_z)

        ts_sample = np.linspace(tmin, tmax, n_samples)
        sampled = np.column_stack(
            [
                np.polyval(coeff_x, ts_sample),
                np.polyval(coeff_y, ts_sample),
                np.polyval(coeff_z, ts_sample),
            ]
        )

        # Build a compound of straight edges between consecutive sampled points
        builder = BRep_Builder()
        comp = TopoDS_Compound()
        builder.MakeCompound(comp)
        last_p = None
        for s in sampled:
            gp = gp_Pnt(float(s[0]), float(s[1]), float(s[2]))
            if last_p is not None:
                e = BRepBuilderAPI_MakeEdge(last_p, gp).Edge()
                builder.Add(comp, e)
            last_p = gp
        poly_compound = comp

    # compute errors (always)
    pts_arr = pts  # already numpy
    origin = np.array(
        [float(centroid[0]), float(centroid[1]), float(centroid[2])], dtype=float
    )
    dir_vec = np.array(
        [float(direction[0]), float(direction[1]), float(direction[2])], dtype=float
    )
    dir_vec = dir_vec / np.linalg.norm(dir_vec)

    vecs = pts_arr - origin
    proj = np.dot(vecs, dir_vec)
    proj_points = np.outer(proj, dir_vec) + origin
    residuals_line = np.linalg.norm(pts_arr - proj_points, axis=1)

    # distances to bounded edge (or polyline compound)
    if poly_compound is not None:
        # compute distance to polyline (min distance to each segment)
        dist_edge = np.array(
            [dist_point_to_polyline(p, sampled) for p in points], dtype=float
        )
    else:
        geom_curve, first_param, last_param = BRep_Tool.Curve(edge)
        dist_edge = np.array(
            [dist_point_to_curve_gp(p, geom_curve, origin) for p in points], dtype=float
        )

    # distances to infinite Geom_Line (OCC) for comparison
    geom_line = Geom_Line(line)
    dist_geom_line = np.array(
        [dist_point_to_curve_gp(p, geom_line, origin) for p in points], dtype=float
    )

    def summarize(a):
        return {
            "mean": float(np.mean(a)),
            "rms": float(np.sqrt(np.mean(a**2))),
            "max": float(np.max(a)),
        }

    errors = {
        "per_point": {
            "line": residuals_line,
            "edge": dist_edge,
            "geom_line": dist_geom_line,
        },
        "summary": {
            "line": summarize(residuals_line),
            "edge": summarize(dist_edge),
            "geom_line": summarize(dist_geom_line),
        },
    }

    # Return a consistent tuple: (shape, gp_Lin, errors, poly_coeffs)
    shape = poly_compound if poly_compound is not None else edge
    return shape, line, errors, poly_coeffs


if __name__ == "__main__":
    # Show the edge in an interactive pythonocc viewer
    from OCC.Display.SimpleGui import init_display

    display, start_display, add_menu, add_function_to_menu = init_display()

    # Simple smoke test: generate noisy points along a known line
    np.random.seed(1)
    base = np.array([1.0, 2.0, -0.5])
    dir_vec = np.array([0.4, 0.2, 0.9])
    dir_vec = dir_vec / np.linalg.norm(dir_vec)
    ts = np.linspace(-5.0, 5.0, 21)
    # Create gp_Pnt list (assume gp_Pnt is available)
    pts = []
    noise_scale = 0.02
    for t in ts:
        p = base + t * dir_vec
        p = p + np.random.normal(scale=noise_scale, size=3)
        pts.append(gp_Pnt(float(p[0]), float(p[1]), float(p[2])))

    shape, line, errors, poly_coeffs = fit_line_edge(pts, extend=0.5, poly_degree=3)
    # Print basic information
    print("Fitted shape type:", type(shape))
    loc = line.Location()
    dir_ = line.Direction()
    print(
        "Line origin:",
        loc.X(),
        loc.Y(),
        loc.Z(),
        "direction:",
        dir_.X(),
        dir_.Y(),
        dir_.Z(),
    )

    # Print summaries for all three error evaluations
    s_line = errors["summary"]["line"]
    s_edge = errors["summary"]["edge"]
    s_geom = errors["summary"]["geom_line"]
    print("\nResiduals (mean / RMS / max):")
    print(
        " infinite PCA-line:",
        f"{s_line['mean']:.6g} / {s_line['rms']:.6g} / {s_line['max']:.6g}",
    )
    print(
        " bounded edge    :",
        f"{s_edge['mean']:.6g} / {s_edge['rms']:.6g} / {s_edge['max']:.6g}",
    )
    print(
        " infinite Geom_Line:",
        f"{s_geom['mean']:.6g} / {s_geom['rms']:.6g} / {s_geom['max']:.6g}",
    )

    # Print per-point arrays and comparison
    per_line = errors["per_point"]["line"]
    per_edge = errors["per_point"]["edge"]
    per_geom = errors["per_point"]["geom_line"]
    print("\nper-point line array:", per_line)
    print("per-point edge array:", per_edge)
    print("per-point geom_line array:", per_geom)

    # print("\nIndex, line, edge, geom_line, diff(edge-line), diff(geom-line)")
    # for i in range(len(per_line)):
    #    d1 = per_edge[i] - per_line[i]
    #    d2 = per_geom[i] - per_line[i]
    #    print(
    #        f"{i:3d}, {per_line[i]:.6g}, {per_edge[i]:.6g}, {per_geom[i]:.6g}, {d1:.6g}, {d2:.6g}"
    #    )
    for p in pts:
        display.DisplayShape(p)
    display.DisplayShape(shape)
    display.FitAll()
    start_display()
