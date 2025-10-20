"""
Scherk's second (singly-periodic) minimal surface visualization using pythonOCC

Creates a BSpline surface by sampling the parametric form
    z = log(cos(v) / cos(u))
on a rectangular domain (u,v) and displays it with the project's `dispocc` helper.

Usage:
    python scherk_surface_occ.py

You can edit domain and resolution at the top of the file or pass via command-line args.
"""

import sys
import os
from math import pi
import argparse
import numpy as np

# bring project helper (same helper used in occ_minimal_surface.py)
sys.path.append(os.path.join("./"))
from base_occ import dispocc

from OCC.Core.gp import gp_Pnt, gp_Ax3, gp_Trsf, gp_Vec
from OCC.Core.TColgp import TColgp_Array2OfPnt
from OCC.Core.GeomAPI import GeomAPI_PointsToBSplineSurface
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.Core.TopLoc import TopLoc_Location


def _bspline_face_from_grid(points_grid, ax=None):
    """Build a BSpline face from a 2D Python list/array of gp_Pnt.

    points_grid should be indexed as points_grid[i][j] with i in [0..nu-1], j in [0..nv-1]
    """
    nu = len(points_grid)
    nv = len(points_grid[0]) if nu > 0 else 0
    arr = TColgp_Array2OfPnt(1, nu, 1, nv)
    for i in range(1, nu + 1):
        for j in range(1, nv + 1):
            arr.SetValue(i, j, points_grid[i - 1][j - 1])
    builder = GeomAPI_PointsToBSplineSurface(arr, 3, 3, 1e-6)
    bs = builder.Surface()
    face = BRepBuilderAPI_MakeFace(bs, 1e-6).Face()
    if ax is not None:
        trf = gp_Trsf()
        try:
            trf.SetTransformation(gp_Ax3(), ax)
        except Exception:
            # fallback: translate to ax origin
            p = ax.Location()
            trf.SetTranslation(gp_Vec(p.X(), p.Y(), p.Z()))
        face = face.Moved(TopLoc_Location(trf))
    return face


def make_scherk_weierstrass(nr=120, ntheta=160, r_max=0.95, ax=None, debug=False):
    """Scherk's second surface via the Wikipedia Weierstrass parametrization.

    Returns a single face built from the (r,theta) grid mapping.
    """
    rs = np.linspace(1e-4, r_max, nr)
    # avoid duplicating the 0 and 2*pi points which creates a degenerate column
    thetas = np.linspace(0, 2 * pi, ntheta, endpoint=False)
    grid = []
    for i, r in enumerate(rs):
        row = []
        for j, th in enumerate(thetas):
            # x = ln((1+r^2+2r cos θ)/(1+r^2-2r cos θ))
            numx = 1 + r * r + 2 * r * np.cos(th)
            denx = 1 + r * r - 2 * r * np.cos(th)
            x = float(np.log(numx / denx))
            # y = ln((1+r^2-2r sin θ)/(1+r^2+2r sin θ))
            numy = 1 + r * r - 2 * r * np.sin(th)
            deny = 1 + r * r + 2 * r * np.sin(th)
            y = float(np.log(numy / deny))
            # z = 2 * atan( (2 r^2 sin 2θ) / (r^4 - 1) )
            denom = r**4 - 1.0
            numer = 2.0 * r**2 * np.sin(2 * th)
            z = float(2.0 * np.arctan2(numer, denom))
            row.append(gp_Pnt(x, y, z))
        grid.append(row)
    face = _bspline_face_from_grid(grid, ax=ax)
    if debug:
        return face, grid
    return face


def make_scherk_implicit_patch(
    x_range=(-1.0, 1.0), y_range=(-1.0, 1.0), nx=120, ny=120, ax=None, t_max=0.95
):
    """Build two patches z = asin(sinh(x)*sinh(y)) and its negative where defined.

    This uses the implicit form sin(z) - sinh(x)*sinh(y) = 0 and solves for z.
    The domain is limited to where |sinh(x)*sinh(y)| <= 1.
    Returns a tuple of faces (top_face, bottom_face) where available.
    """
    xs = np.linspace(x_range[0], x_range[1], nx)
    ys = np.linspace(y_range[0], y_range[1], ny)

    # Compute validity mask where a real solution exists and is away from singularity (|t| <= t_max)
    T = np.empty((len(xs), len(ys)), dtype=float)
    valid = np.zeros_like(T, dtype=bool)
    for i, xi in enumerate(xs):
        for j, yj in enumerate(ys):
            t = np.sinh(xi) * np.sinh(yj)
            T[i, j] = t
            valid[i, j] = abs(t) <= t_max

    # Find largest rectangular subgrid that contains valid points only by trimming
    rows_any = valid.any(axis=1)
    cols_any = valid.any(axis=0)
    if not rows_any.any() or not cols_any.any():
        # No valid region at all, return tiny flat patches at z=0
        face0 = _bspline_face_from_grid([[gp_Pnt(0, 0, 0)]], ax=ax)
        return face0, face0

    i0 = int(np.argmax(rows_any))
    i1 = len(rows_any) - int(np.argmax(rows_any[::-1]))
    j0 = int(np.argmax(cols_any))
    j1 = len(cols_any) - int(np.argmax(cols_any[::-1]))

    xs_sub = xs[i0:i1]
    ys_sub = ys[j0:j1]

    # build point grids only inside valid rectangular region
    grid_top = []
    grid_bot = []
    for xi in xs_sub:
        row_t = []
        row_b = []
        for yj in ys_sub:
            t = np.sinh(xi) * np.sinh(yj)
            # clip to safe range inside [-t_max, t_max] for stability
            t = np.clip(t, -t_max + 1e-9, t_max - 1e-9)
            z = float(np.arcsin(t))
            row_t.append(gp_Pnt(xi, yj, z))
            row_b.append(gp_Pnt(xi, yj, -z))
        grid_top.append(row_t)
        grid_bot.append(row_b)

    top_face = _bspline_face_from_grid(grid_top, ax=ax)
    bot_face = _bspline_face_from_grid(grid_bot, ax=ax)
    return top_face, bot_face


def make_scherk_surface(
    u_min=-1.0, u_max=1.0, v_min=-1.0, v_max=1.0, nu=80, nv=80, scale_xy=1.0, ax=None
):
    """Build a BSpline approximation of Scherk's second surface.

    Parametric formula used:
        z = log(cos(v) / cos(u))

    The domain for u and v must lie inside (-pi/2, pi/2) to avoid cos=0.

    Returns: (bsurf, face)
    """
    # small clamp to avoid cos() == 0
    clamp = 1e-3
    u_lo = max(u_min, -pi / 2 + clamp)
    u_hi = min(u_max, pi / 2 - clamp)
    v_lo = max(v_min, -pi / 2 + clamp)
    v_hi = min(v_max, pi / 2 - clamp)

    arr = TColgp_Array2OfPnt(1, nu, 1, nv)
    for i in range(1, nu + 1):
        u = u_lo + (i - 1) / (nu - 1) * (u_hi - u_lo)
        for j in range(1, nv + 1):
            v = v_lo + (j - 1) / (nv - 1) * (v_hi - v_lo)
            # compute height
            z = float(np.log(np.cos(v) / np.cos(u)))
            x = u * scale_xy
            y = v * scale_xy
            arr.SetValue(i, j, gp_Pnt(x, y, z))

    # Create a BSpline surface that interpolates the points
    builder = GeomAPI_PointsToBSplineSurface(arr, 3, 3, 1e-6)
    bs = builder.Surface()
    # Use the overload that accepts (handle<Geom_Surface>, tolerance)
    face = BRepBuilderAPI_MakeFace(bs, 1e-6).Face()
    if ax is not None:
        trf = gp_Trsf()
        try:
            trf.SetTransformation(gp_Ax3(), ax)
        except Exception:
            p = ax.Location()
            trf.SetTranslation(gp_Vec(p.X(), p.Y(), p.Z()))
        face = face.Moved(TopLoc_Location(trf))
    return bs, face


if __name__ == "__main__":
    # --- 決め打ちパラメータ (ここを編集すれば固定値を変更できます) ---
    U_RANGE = (-1.0, 1.0)
    V_RANGE = (-1.0, 1.0)
    NU = 120  # u方向のサンプリング数 (高くすると滑らかになります)
    NV = 120  # v方向のサンプリング数
    SCALE = 1.0

    # Create viewer
    obj = dispocc(touch=True)

    # Demonstrate placement using gp_Ax3: create three Ax3 and display each example at a different location
    from OCC.Core.gp import gp_Ax3, gp_Pnt, gp_Dir

    ax_graph = gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1))
    ax_weier = gp_Ax3(gp_Pnt(6.0, 0, 0), gp_Dir(0, 0, 1))
    ax_implicit = gp_Ax3(gp_Pnt(0, -6.0, 0), gp_Dir(0, 0, 1))

    # 1) Graph form
    bs, face1 = make_scherk_surface(
        u_min=U_RANGE[0],
        u_max=U_RANGE[1],
        v_min=V_RANGE[0],
        v_max=V_RANGE[1],
        nu=NU,
        nv=NV,
        scale_xy=SCALE,
        ax=ax_graph,
    )
    obj.display.DisplayShape(face1)

    # 2) Weierstrass parameterisation placed at ax_weier
    face2, grid2 = make_scherk_weierstrass(
        nr=120, ntheta=200, r_max=0.95, ax=ax_weier, debug=True
    )
    obj.display.DisplayShape(face2, color="BLUE1", transparency=0.5)
    # display a sparse subset of sample points to inspect ordering
    # for i in range(0, len(grid2), max(1, len(grid2) // 12)):
    #    for j in range(0, len(grid2[0]), max(1, len(grid2[0]) // 24)):
    #        p = grid2[i][j]
    #        obj.display.DisplayShape(p, color="BLUE1")

    # 3) Implicit patch (two sheets) placed at ax_implicit
    top, bot = make_scherk_implicit_patch(
        x_range=(-2.0, 2.0),
        y_range=(-2.0, 2.0),
        nx=100,
        ny=100,
        ax=ax_implicit,
        t_max=0.99,
    )
    obj.display.DisplayShape(top, color="RED", transparency=0.5)
    obj.display.DisplayShape(bot, color="GREEN", transparency=0.5)

    obj.display.FitAll()
    obj.ShowOCC()
