#!/usr/bin/env python3
"""
occ_3d_polygon.py

Generate ~10 points (circle-like), give some non-zero Z, triangulate and display using pythonocc.

Features:
- generate ~10 3D points on an approximate circle
- triangulate by proximity (Delaunay if available, ear-clipping fallback)
- approximate geodesic Voronoi on the triangle mesh via multi-source Dijkstra
- display the mesh and Voronoi regions using pythonocc (no plotly/shapely for display)
"""

import numpy as np
import random
import heapq
import colorsys
import time
import os
import time
from math import cos, sin, pi
from collections import defaultdict
from scipy.spatial import Delaunay, cKDTree

# Optional dependency: shapely for robust polygon merging (nicer Voronoi boundaries)
from shapely.geometry import Polygon, MultiPolygon
from shapely.ops import unary_union
from OCC.Core.gp import gp_Pnt
from OCC.Core.BRepBuilderAPI import (
    BRepBuilderAPI_MakePolygon,
    BRepBuilderAPI_MakeFace,
    BRepBuilderAPI_MakeVertex,
    BRepBuilderAPI_Sewing,
)
from OCC.Display.SimpleGui import init_display
from OCC.Core.BRep import BRep_Builder
from OCC.Core.TopoDS import TopoDS_Compound
from OCC.Core.Quantity import Quantity_Color, Quantity_TOC_RGB

# small palette of colors (strings accepted by DisplayShape)
PALETTE = [
    "RED",
    "GREEN",
    "BLUE1",
    "YELLOW",
    "BROWN",
    "ORANGE",
    "PINK",
    "CYAN",
    "MAGENTA",
    "DARKGREEN",
]

# Diagnostics: record bisector solve failures (for edges where denom ~ 0)
bisector_failures = []


def generate_points(n=10, radius=50.0, z_variation=10.0, nonzero_prob=0.4, seed=None):
    if seed is not None:
        random.seed(seed)
    pts = []
    for i in range(n):
        a = 2.0 * pi * i / n
        r = radius * (1.0 + random.uniform(-0.2, 0.2))
        x = r * cos(a)
        y = r * sin(a)
        z = 0.0
        if random.random() < nonzero_prob:
            z = random.uniform(-z_variation, z_variation)
        pts.append((x, y, z))
    return pts


def make_outline_wire(points):
    poly = BRepBuilderAPI_MakePolygon()
    for p in points:
        poly.Add(gp_Pnt(*p))
    poly.Close()
    return poly.Wire()


def make_tri_faces_from_centroid(points):
    faces = []
    n = len(points)
    if n < 3:
        return faces
    cx = sum(p[0] for p in points) / n
    cy = sum(p[1] for p in points) / n
    cz = sum(p[2] for p in points) / n
    center = (cx, cy, cz)
    for i in range(n):
        a = points[i]
        b = points[(i + 1) % n]
        poly = BRepBuilderAPI_MakePolygon()
        poly.Add(gp_Pnt(*center))
        poly.Add(gp_Pnt(*a))
        poly.Add(gp_Pnt(*b))
        poly.Close()
        face = BRepBuilderAPI_MakeFace(poly.Wire()).Face()
        faces.append(face)
    return faces


def make_tri_faces_by_proximity(points):
    """Triangulate by proximity: try scipy.spatial.Delaunay on XY projection,
    fall back to an ear-clipping triangulation on the 2D polygon order.
    Returns list of TopoDS_Face objects.
    """
    # Project to XY
    xy = [(p[0], p[1]) for p in points]

    # Try Delaunay if available
    pts2 = np.array(xy)
    tri = Delaunay(pts2)
    faces = []
    for simplex in tri.simplices:
        i, j, k = int(simplex[0]), int(simplex[1]), int(simplex[2])
        poly = BRepBuilderAPI_MakePolygon()
        poly.Add(gp_Pnt(*points[i]))
        poly.Add(gp_Pnt(*points[j]))
        poly.Add(gp_Pnt(*points[k]))
        poly.Close()
        faces.append(BRepBuilderAPI_MakeFace(poly.Wire()).Face())
    return faces


def build_mesh_delaunay(points_xy, points3d=None):
    """Return vertices list and triangle index list using Delaunay on XY.

    points_xy: list of (x,y) or if points3d provided, points_xy may be derived.
    points3d: optional list of 3D coords corresponding to points_xy
    """
    if Delaunay is None or np is None:
        raise RuntimeError("scipy.spatial.Delaunay required for build_mesh_delaunay")
    pts2 = np.array(points_xy)
    tri = Delaunay(pts2)
    vertices = (
        points3d
        if points3d is not None
        else [(float(x), float(y), 0.0) for x, y in points_xy]
    )
    triangles = [tuple(map(int, s)) for s in tri.simplices]
    return vertices, triangles


def geodesic_voronoi(vertices, triangles, seed_vertex_indices):
    """Approximate geodesic Voronoi by Dijkstra on vertex graph.

    vertices: list of (x,y,z)
    triangles: list of (i,j,k) vertex indices
    seed_vertex_indices: list of vertex indices used as seeds

    Returns: label per vertex (list), distances per vertex (list)
    """
    n = len(vertices)
    # build adjacency
    adj = [[] for _ in range(n)]
    for i, j, k in triangles:
        for a, b in ((i, j), (j, k), (k, i)):
            # add undirected edge with length
            dx = vertices[a][0] - vertices[b][0]
            dy = vertices[a][1] - vertices[b][1]
            dz = vertices[a][2] - vertices[b][2]
            w = (dx * dx + dy * dy + dz * dz) ** 0.5
            adj[a].append((b, w))
            adj[b].append((a, w))

    dist = [float("inf")] * n
    label = [-1] * n
    pq = []
    # initialize seeds
    for s in seed_vertex_indices:
        if 0 <= s < n:
            dist[s] = 0.0
            label[s] = s
            heapq.heappush(pq, (0.0, s, s))  # (dist, vertex, seed_label)

    while pq:
        d, u, seed_lab = heapq.heappop(pq)
        if d > dist[u]:
            continue
        # propagate
        for v, w in adj[u]:
            nd = d + w
            if nd < dist[v]:
                dist[v] = nd
                label[v] = seed_lab
                heapq.heappush(pq, (nd, v, seed_lab))

    return label, dist


def sample_points_on_mesh(vertices, triangles, m=200, seed=None):
    """Sample m points uniformly on the triangle mesh (area-weighted).

    vertices: list of (x,y,z)
    triangles: list of (i,j,k)
    Returns: list of tuples ((x,y,z), triangle_index)
    """
    if seed is not None:
        random.seed(seed)
    # compute triangle areas
    areas = []
    tri_verts = []
    for i, j, k in triangles:
        vi = vertices[i]
        vj = vertices[j]
        vk = vertices[k]
        # vectors
        ux, uy, uz = vj[0] - vi[0], vj[1] - vi[1], vj[2] - vi[2]
        vx, vy, vz = vk[0] - vi[0], vk[1] - vi[1], vk[2] - vi[2]
        # cross product magnitude / 2
        cx = uy * vz - uz * vy
        cy = uz * vx - ux * vz
        cz = ux * vy - uy * vx
        area = 0.5 * (cx * cx + cy * cy + cz * cz) ** 0.5
        areas.append(area)
        tri_verts.append((vi, vj, vk))
    total = sum(areas)
    if total <= 0:
        return []
    # cumulative distribution
    cum = []
    s = 0.0
    for a in areas:
        s += a
        cum.append(s / total)

    samples = []
    for _ in range(m):
        r = random.random()
        # find triangle index
        idx = 0
        while idx < len(cum) and r > cum[idx]:
            idx += 1
        if idx >= len(triangles):
            idx = len(triangles) - 1
        vi, vj, vk = tri_verts[idx]
        # sample barycentric coordinates
        u = random.random()
        v = random.random()
        if u + v > 1.0:
            u = 1.0 - u
            v = 1.0 - v
        w = 1.0 - u - v
        x = u * vi[0] + v * vj[0] + w * vk[0]
        y = u * vi[1] + v * vj[1] + w * vk[1]
        z = u * vi[2] + v * vj[2] + w * vk[2]
        samples.append(((x, y, z), idx))
    return samples


def map_points_to_nearest_vertex(points, vertices):
    """Return list of nearest vertex indices for each point.

    Uses cKDTree if available for speed, otherwise brute-force.
    """
    if cKDTree is not None and np is not None:
        coords = np.array(vertices)
        tree = cKDTree(coords)
        pts = np.array(points)
        dists, idxs = tree.query(pts)
        return [int(i) for i in idxs]
    # brute force
    out = []
    for p in points:
        best = 0
        bd = None
        for i, v in enumerate(vertices):
            dx = p[0] - v[0]
            dy = p[1] - v[1]
            dz = p[2] - v[2]
            d = dx * dx + dy * dy + dz * dz
            if bd is None or d < bd:
                bd = d
                best = i
        out.append(best)
    return out


def triangles_by_label(triangles, vertex_labels, vertices=None):
    """Assign each triangle a single label to avoid overlapping regions.

    Strategy:
    - Prefer majority vote of the three vertex labels (keeps assignment on mesh).
    - If all three vertex labels differ, fall back to centroid->nearest-vertex mapping.
    This guarantees each triangle is uniquely assigned to exactly one region and
    keeps region boundaries along original mesh edges.
    """
    regions = defaultdict(list)
    for tri in triangles:
        i, j, k = tri
        labs = [vertex_labels[i], vertex_labels[j], vertex_labels[k]]
        # majority vote: if at least two vertices share a label, use that
        if labs[0] == labs[1] or labs[0] == labs[2]:
            lab = labs[0]
        elif labs[1] == labs[2]:
            lab = labs[1]
        else:
            # all three labels different: fall back to centroid->nearest-vertex mapping
            if vertices is None:
                # fallback to first label
                lab = labs[0]
            else:
                vi = vertices[i]
                vj = vertices[j]
                vk = vertices[k]
                cx = (vi[0] + vj[0] + vk[0]) / 3.0
                cy = (vi[1] + vj[1] + vk[1]) / 3.0
                cz = (vi[2] + vj[2] + vk[2]) / 3.0
                nearest_vs = map_points_to_nearest_vertex([(cx, cy, cz)], vertices)
                lab = vertex_labels[nearest_vs[0]]
        regions[lab].append(tri)

    return regions


def boundary_loops_from_triangles(vertices, tri_list):
    """Compute boundary loops (as lists of 3D points) from a list of triangles.

    This uses only original mesh edges: an undirected edge that appears exactly
    once in the tri_list is considered a boundary edge. The function returns a
    list of closed loops (each loop is a list of 3D points). No new vertices
    are created and nothing is re-triangulated.
    """
    if not tri_list:
        return []
    edge_count = {}
    edge_map = {}
    for i, j, k in tri_list:
        for a, b in ((i, j), (j, k), (k, i)):
            e = (min(a, b), max(a, b))
            edge_count[e] = edge_count.get(e, 0) + 1
            # remember the directed orientation seen last (used to draw)
            edge_map[e] = (a, b)

    # boundary edges are those with count == 1
    boundary_edges = [edge_map[e] for e, c in edge_count.items() if c == 1]

    # build adjacency for traversal
    adj = {}
    for a, b in boundary_edges:
        ai = int(a)
        bi = int(b)
        adj.setdefault(ai, []).append(bi)
        adj.setdefault(bi, []).append(ai)

    loops = []
    visited = set()
    for start in list(adj.keys()):
        if start in visited:
            continue
        loop = [start]
        visited.add(start)
        cur = start
        prev = None
        # walk the boundary greedily
        while True:
            nbrs = [n for n in adj.get(cur, []) if n != prev]
            if not nbrs:
                break
            nxt = nbrs[0]
            if nxt == start:
                break
            loop.append(nxt)
            visited.add(nxt)
            prev, cur = cur, nxt
        if len(loop) >= 3:
            coords = [vertices[v] for v in loop]
            loops.append(coords)
    return loops


def merge_triangles_to_polygons_shapely(vertices, tri_list):
    """Merge a set of triangles (given by vertex indices) into planar polygons using Shapely.

    Returns list of loops where each loop is a list of 3D points (ordered exterior coords).
    Raises if Shapely not available.
    """
    if not tri_list:
        return []
    # collect unique vertices used by these triangles
    used_idx = set()
    for i, j, k in tri_list:
        used_idx.update((i, j, k))
    pts = [vertices[i] for i in sorted(used_idx)]
    # best-fit plane for this label's points
    bf = _best_fit_plane(pts)
    origin = None
    u = None
    v = None
    if bf is not None:
        origin, u, v, normal, max_res = bf
    else:
        # fallback to XY plane
        origin = (0.0, 0.0, 0.0)
        u = (1.0, 0.0, 0.0)
        v = (0.0, 1.0, 0.0)

    # build shapely polygons in plane coordinates
    poly_list = []
    for i, j, k in tri_list:
        tri_pts = [vertices[i], vertices[j], vertices[k]]
        coords2 = []
        for p in tri_pts:
            vec = np.array(p) - np.array(origin)
            x = float(np.dot(vec, np.array(u)))
            y = float(np.dot(vec, np.array(v)))
            coords2.append((x, y))
        try:
            shp = Polygon(coords2)
            if shp.is_valid and not shp.is_empty and shp.area > 0:
                poly_list.append(shp)
        except Exception:
            # skip malformed triangle
            continue

    if not poly_list:
        return []

    merged = unary_union(poly_list)
    results = []
    if isinstance(merged, Polygon):
        geoms = [merged]
    elif isinstance(merged, MultiPolygon):
        geoms = list(merged.geoms)
    else:
        return []

    for g in geoms:
        exterior = list(g.exterior.coords)
        loop3 = []
        for x, y in exterior:
            p3 = (
                origin[0] + x * u[0] + y * v[0],
                origin[1] + x * u[1] + y * v[1],
                origin[2] + x * u[2] + y * v[2],
            )
            loop3.append(p3)
        results.append(loop3)
    return results


def _interp_edge_bisector(pa, pb, sa, sb, eps=1e-12):
    """Compute point along edge pa->pb where Euclidean distance to seed sa and sb are equal.

    Returns point (x,y,z) or None if denominator is too small. Uses linear solve derived from
    |pa + t*(pb-pa) - sa|^2 == |pa + t*(pb-pa) - sb|^2.
    """
    vx = pb[0] - pa[0]
    vy = pb[1] - pa[1]
    vz = pb[2] - pa[2]
    d0x = pa[0] - sa[0]
    d0y = pa[1] - sa[1]
    d0z = pa[2] - sa[2]
    d1x = pa[0] - sb[0]
    d1y = pa[1] - sb[1]
    d1z = pa[2] - sb[2]
    # numerator = d1^2 - d0^2
    n = (d1x * d1x + d1y * d1y + d1z * d1z) - (d0x * d0x + d0y * d0y + d0z * d0z)
    # denom = 2 * v . (d0 - d1)
    ddx = d0x - d1x
    ddy = d0y - d1y
    ddz = d0z - d1z
    denom = 2.0 * (vx * ddx + vy * ddy + vz * ddz)
    if abs(denom) < eps:
        # record the failure (store XY endpoints) and return None to indicate instability
        try:
            bisector_failures.append(((pa[0], pa[1]), (pb[0], pb[1])))
        except Exception:
            bisector_failures.append((None, None))
        # do NOT return midpoint here; let caller decide fallback/absorption
        return None
    t = n / denom
    if t < -1e-6 or t > 1.0 + 1e-6:
        # intersection outside edge (shouldn't happen often); clamp
        t = max(0.0, min(1.0, t))
    return (pa[0] + t * vx, pa[1] + t * vy, pa[2] + t * vz)


def partition_triangle_by_vertex_labels(vertices, triangle, vertex_labels):
    """Partition a single triangle into polygon pieces per label using edge intersections.

    Returns a dict: label -> list of polygons (each polygon is list of 3D points).
    """
    i, j, k = triangle
    inds = (i, j, k)
    labs = [vertex_labels[i], vertex_labels[j], vertex_labels[k]]

    # 1) simple case: all three vertices same label -> whole triangle belongs to that label
    if labs[0] == labs[1] == labs[2]:
        return {labs[0]: [[vertices[i], vertices[j], vertices[k]]]}

    # compute triangle plane for robust projection / ordering
    tri_pts = [vertices[i], vertices[j], vertices[k]]
    bf = _best_fit_plane(tri_pts)
    if bf is not None:
        origin, u, v, normal, max_res = bf
    else:
        # fallback to XY plane
        origin = (0.0, 0.0, 0.0)
        u = (1.0, 0.0, 0.0)
        v = (0.0, 1.0, 0.0)

    # 2) compute bisector points on edges where labels differ; skip unstable edges (None)
    edge_pairs = [(i, j), (j, k), (k, i)]
    edge_intersections = {}
    for a, b in edge_pairs:
        La = vertex_labels[a]
        Lb = vertex_labels[b]
        if La == Lb:
            continue
        pa = vertices[a]
        pb = vertices[b]
        sa = vertices[La]
        sb = vertices[Lb]
        p = _interp_edge_bisector(pa, pb, sa, sb)
        if p is None:
            # unstable bisector; skip adding a bisector point here
            continue
        # store one bisector per undirected edge
        edge_intersections[(min(a, b), max(a, b))] = (p, La, Lb)

    # 3) collect points for each label: triangle vertices that belong to label + bisector points adjacent to label
    pieces = defaultdict(list)
    for idx in inds:
        lab = vertex_labels[idx]
        pieces[lab].append(vertices[idx])
    for (a, b), (p, La, Lb) in edge_intersections.items():
        pieces[La].append(p)
        pieces[Lb].append(p)

    # 4) for each label, project to plane, dedupe by rounded plane coords, sort by angle, and reconstruct coplanar 3D pts
    ordered = {}
    for lab, pts in pieces.items():
        if not pts:
            continue
        proj = []
        seen = set()
        # project points to (x,y) in triangle plane and dedupe by rounded coords
        # choose rounding digits relative to scale
        scale = max(1.0, max(abs(c) for p in tri_pts for c in p))
        digits = 9
        for q in pts:
            vec = np.array(q) - np.array(origin)
            x = float(np.dot(vec, np.array(u)))
            y = float(np.dot(vec, np.array(v)))
            key = (round(x, digits), round(y, digits))
            if key in seen:
                continue
            seen.add(key)
            proj.append((x, y, q))

        if len(proj) < 3:
            # If we don't have enough distinct points to form a polygon,
            # add the triangle centroid (projected) as a fallback. This
            # ensures we can always form a minimal polygon fragment for
            # this label when appropriate.
            tv0, tv1, tv2 = tri_pts
            vec0 = np.array(tv0) - np.array(origin)
            vec1 = np.array(tv1) - np.array(origin)
            vec2 = np.array(tv2) - np.array(origin)
            tx = (
                float(np.dot(vec0, np.array(u)))
                + float(np.dot(vec1, np.array(u)))
                + float(np.dot(vec2, np.array(u)))
            ) / 3.0
            ty = (
                float(np.dot(vec0, np.array(v)))
                + float(np.dot(vec1, np.array(v)))
                + float(np.dot(vec2, np.array(v)))
            ) / 3.0
            tri_centroid_3d = (
                (tv0[0] + tv1[0] + tv2[0]) / 3.0,
                (tv0[1] + tv1[1] + tv2[1]) / 3.0,
                (tv0[2] + tv1[2] + tv2[2]) / 3.0,
            )
            tkey = (round(tx, digits), round(ty, digits))
            if tkey not in seen:
                seen.add(tkey)
                proj.append((tx, ty, tri_centroid_3d))
        if len(proj) < 3:
            # still cannot form polygon
            continue

        # centroid in plane coords
        cx = sum(p[0] for p in proj) / len(proj)
        cy = sum(p[1] for p in proj) / len(proj)

        def ang(item):
            x, y, _ = item
            return np.arctan2(y - cy, x - cx)

        proj.sort(key=ang)

        ordered_pts = []
        # Use the original 3D points (stored as the third element) to avoid
        # reconstructing from plane coordinates which can drift off the
        # original triangle surface when merged across multiple triangles.
        for x, y, q in proj:
            ordered_pts.append(q)
        ordered[lab] = ordered_pts

    return ordered


def build_voronoi_regions_by_splitting_triangles(vertices, triangles, vertex_labels):
    """Build region polygons by splitting each triangle according to vertex_labels.

    Returns dict label -> list of polygons (each polygon is list of 3D points).
    """
    regions = defaultdict(list)
    for tri in triangles:
        parts = partition_triangle_by_vertex_labels(vertices, tri, vertex_labels)
        for lab, polys in parts.items():
            for poly in polys:
                regions[lab].append(poly)
    return regions


def normalize_polygon(poly):
    """Return polygon as list of 3-tuples (x,y,z).

    Handles cases where poly may be a flat list of numbers or numpy arrays.
    """
    out = []
    # case: empty
    if not poly:
        return out
    # if poly is a flat list of numbers
    if all(isinstance(x, (int, float, np.floating, np.integer)) for x in poly):
        if len(poly) % 3 != 0:
            return []
        for i in range(0, len(poly), 3):
            out.append((float(poly[i]), float(poly[i + 1]), float(poly[i + 2])))
        return out
    # otherwise expect iterable of 3-length items
    for v in poly:
        if v is None:
            continue
        try:
            if len(v) == 3:
                out.append((float(v[0]), float(v[1]), float(v[2])))
                continue
        except Exception:
            pass
        # fallback: ignore malformed entries
    return out


def _best_fit_plane(pts):
    """Return (origin, u, v, normal, max_residual)

    origin: centroid
    u, v: orthonormal basis in plane
    normal: unit normal
    max_residual: maximum absolute distance of points from plane
    """
    if not pts:
        return None
    arr = np.array([(p[0], p[1], p[2]) for p in pts], dtype=float)
    centroid = arr.mean(axis=0)
    M = arr - centroid
    # SVD: principal components
    try:
        U, S, Vt = np.linalg.svd(M, full_matrices=False)
    except Exception:
        return None
    # Vt rows are principal directions; normal is last row
    if Vt.shape[0] < 3:
        return None
    u = Vt[0]
    v = Vt[1]
    normal = Vt[2]
    # ensure orthonormal
    # compute residuals (distance to plane along normal)
    dists = np.abs(np.dot(M, normal))
    max_res = float(np.max(dists))
    return (tuple(centroid), tuple(u), tuple(v), tuple(normal), max_res)


def _project_and_order_on_plane(pts, origin, u, v):
    """Project 3D pts to plane basis (u,v) at origin, return ordered 3D pts around centroid in plane coords."""
    coords = []
    for p in pts:
        vec = np.array(p) - np.array(origin)
        x = float(np.dot(vec, np.array(u)))
        y = float(np.dot(vec, np.array(v)))
        coords.append((x, y, p))
    # compute centroid in plane coords
    cx = sum(c[0] for c in coords) / len(coords)
    cy = sum(c[1] for c in coords) / len(coords)

    def ang(c):
        return np.arctan2(c[1] - cy, c[0] - cx)

    coords.sort(key=ang)
    # reconstruct ordered 3D points from plane coords to ensure exact coplanarity
    ordered = []
    for x, y, _ in coords:
        p3 = (
            origin[0] + x * u[0] + y * v[0],
            origin[1] + x * u[1] + y * v[1],
            origin[2] + x * u[2] + y * v[2],
        )
        ordered.append(p3)
    return ordered


def _polygon_area3d(pts):
    """Compute polygon area by projecting to its best-fit plane and using shoelace."""
    if not pts or len(pts) < 3:
        return 0.0
    bf = _best_fit_plane(pts)
    if bf is None:
        return 0.0
    origin, u, v, normal, max_res = bf
    coords = []
    for p in pts:
        vec = np.array(p) - np.array(origin)
        x = float(np.dot(vec, np.array(u)))
        y = float(np.dot(vec, np.array(v)))
        coords.append((x, y))
    area = 0.0
    for i in range(len(coords)):
        x1, y1 = coords[i]
        x2, y2 = coords[(i + 1) % len(coords)]
        area += x1 * y2 - x2 * y1
    return abs(area) * 0.5


def _polygon_centroid3d(pts):
    if not pts:
        return (0.0, 0.0, 0.0)
    arr = np.array(pts, dtype=float)
    c = tuple(arr.mean(axis=0))
    return c


def merge_small_fragments(regions, vertices, tri_area_mean, eps_rel=0.01, eps_abs=1e-6):
    """Merge very small fragments into nearest larger fragments to remove slivers.

    regions: dict label -> list of polygon (3D pts)
    tri_area_mean: mean triangle area (used to set relative threshold)
    Returns: new_regions (same structure), and logs merge counts via prints.
    """
    # flatten fragments
    frags = []
    for lab, polys in regions.items():
        for poly in polys:
            pts = normalize_polygon(poly)
            area = _polygon_area3d(pts)
            centroid = _polygon_centroid3d(pts)
            frags.append(
                {
                    "label": lab,
                    "pts": pts,
                    "area": area,
                    "centroid": centroid,
                }
            )

    if not frags:
        return regions

    # threshold
    thresh = max(eps_abs, eps_rel * max(1e-12, tri_area_mean))

    # separate small and large
    large = [f for f in frags if f["area"] >= thresh]
    small = [f for f in frags if f["area"] < thresh]

    if not small:
        return regions

    # Build KD-tree of large centroids for nearest neighbor search (simple brute-force ok for small counts)
    large_centroids = [f["centroid"] for f in large]

    merged_count = 0
    for s in small:
        # find nearest large fragment by Euclidean distance
        best_idx = None
        best_dist = None
        sx, sy, sz = s["centroid"]
        for idx, lc in enumerate(large_centroids):
            dx = lc[0] - sx
            dy = lc[1] - sy
            dz = lc[2] - sz
            d2 = dx * dx + dy * dy + dz * dz
            if best_dist is None or d2 < best_dist:
                best_dist = d2
                best_idx = idx

        if best_idx is None:
            # no large fragments? skip
            continue

        # merge small into chosen large: append points and recompute area/centroid
        target = large[best_idx]
        # append small pts to target pts (simple concatenation)
        target["pts"].extend(s["pts"])
        # clean: dedupe projected points on target (use plane of target)
        # recompute centroid and area
        target["area"] = _polygon_area3d(target["pts"])
        target["centroid"] = _polygon_centroid3d(target["pts"])
        large_centroids[best_idx] = target["centroid"]
        merged_count += 1

    # rebuild regions dict from large list
    new_regions = defaultdict(list)
    for t in large:
        # try to order/proj the merged pts on best-fit plane for neatness
        pts = t["pts"]
        if len(pts) < 3:
            continue
        bf = _best_fit_plane(pts)
        if bf is not None:
            origin, u, v, normal, max_res = bf
            # project, dedupe, sort by angle
            proj = []
            seen = set()
            digits = 9
            for q in pts:
                vec = np.array(q) - np.array(origin)
                x = float(np.dot(vec, np.array(u)))
                y = float(np.dot(vec, np.array(v)))
                key = (round(x, digits), round(y, digits))
                if key in seen:
                    continue
                seen.add(key)
                proj.append((x, y, q))
            if len(proj) < 3:
                continue
            cx = sum(p[0] for p in proj) / len(proj)
            cy = sum(p[1] for p in proj) / len(proj)
            proj.sort(key=lambda item: np.arctan2(item[1] - cy, item[0] - cx))
            ordered = []
            # Use original 3D points stored in proj entries to preserve exact
            # surface positions rather than reprojecting which can produce
            # slightly off-surface vertices when merging fragments from
            # different triangles.
            for x, y, q in proj:
                ordered.append(q)
            new_regions[t["label"]].append(ordered)
        else:
            # can't order, just keep unique pts
            uniq = []
            seen = set()
            for q in pts:
                key = (round(q[0], 9), round(q[1], 9), round(q[2], 9))
                if key in seen:
                    continue
                seen.add(key)
                uniq.append(q)
            if len(uniq) >= 3:
                new_regions[t["label"]].append(uniq)

    print(f"DEBUG: merged_small_fragments merged_count={merged_count}, thresh={thresh}")
    return new_regions


def _ear_clip_triangulation(points):
    """Simple ear-clipping on polygon vertices in given order (2D XY projection).

    Assumes polygon is simple and vertices are ordered counter-clockwise.
    """

    def area(a, b, c):
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

    pts2 = [(p[0], p[1]) for p in points]
    n = len(pts2)
    if n < 3:
        return []
    # vertex index list
    V = list(range(n))
    faces = []
    counter = 0
    while len(V) > 3 and counter < n * n:
        removed = False
        for i in range(len(V)):
            i_prev = V[(i - 1) % len(V)]
            i_curr = V[i]
            i_next = V[(i + 1) % len(V)]
            a = pts2[i_prev]
            b = pts2[i_curr]
            c = pts2[i_next]
            if area(a, b, c) <= 0:
                continue  # not an ear if non-positive area (assuming CCW)
            # check no other points inside triangle abc
            ear = True
            for vi in V:
                if vi in (i_prev, i_curr, i_next):
                    continue
                p = pts2[vi]
                # barycentric / area tests
                a1 = area(p, a, b)
                a2 = area(p, b, c)
                a3 = area(p, c, a)
                if a1 >= 0 and a2 >= 0 and a3 >= 0:
                    ear = False
                    break
            if ear:
                # clip ear
                poly = BRepBuilderAPI_MakePolygon()
                poly.Add(gp_Pnt(*points[i_prev]))
                poly.Add(gp_Pnt(*points[i_curr]))
                poly.Add(gp_Pnt(*points[i_next]))
                poly.Close()
                faces.append(BRepBuilderAPI_MakeFace(poly.Wire()).Face())
                V.pop(i)
                removed = True
                break
        if not removed:
            # cannot find ear (maybe degenerate), stop to avoid infinite loop
            break
        counter += 1
    if len(V) == 3:
        i0, i1, i2 = V[0], V[1], V[2]
        poly = BRepBuilderAPI_MakePolygon()
        poly.Add(gp_Pnt(*points[i0]))
        poly.Add(gp_Pnt(*points[i1]))
        poly.Add(gp_Pnt(*points[i2]))
        poly.Close()
        faces.append(BRepBuilderAPI_MakeFace(poly.Wire()).Face())
    return faces


def sew_faces(faces, tol=1e-6):
    sew = BRepBuilderAPI_Sewing(tol)
    for f in faces:
        sew.Add(f)
    sew.Perform()
    return sew.SewedShape()


if __name__ == "__main__":
    # fixed parameters (no argparse) per user request
    n = 50
    seed = 123
    radius = 50.0
    zvar = 10.0

    pts = generate_points(n=n, radius=radius, z_variation=zvar, seed=seed)
    wire = make_outline_wire(pts)

    # triangulate (faces for display)
    faces = make_tri_faces_by_proximity(pts)
    sewn = sew_faces(faces)

    display, start_display, add_menu, add_function_to_menu = init_display()
    # don't display the full sewn face (it can occlude region faces)
    # display.DisplayShape(sewn, update=False)

    # Attempt to build a Delaunay mesh and compute geodesic Voronoi on vertices
    points_xy = [(p[0], p[1]) for p in pts]
    vertices, triangles = build_mesh_delaunay(points_xy, points3d=pts)

    # sample a point cloud on the surface (returning triangle index)
    samples = sample_points_on_mesh(vertices, triangles, m=20, seed=seed)
    # display sampled points (as small yellow vertices)
    for sp, tidx in samples:
        vtx = BRepBuilderAPI_MakeVertex(gp_Pnt(*sp)).Vertex()
        display.DisplayShape(vtx)

    # Re-mesh approach: create a new triangulation that includes all samples as vertices.
    # 1) append sample coords to vertex list
    sample_coords = [sp for sp, _ in samples]
    orig_vert_count = len(vertices)
    all_vertices = list(vertices) + list(sample_coords)
    all_xy = [(v[0], v[1]) for v in all_vertices]

    # 2) run Delaunay on combined XY to produce a consistent mesh
    new_vertices, new_triangles = build_mesh_delaunay(all_xy, points3d=all_vertices)

    # 3) filter triangles by whether their centroid is inside the original polygon outline
    poly_xy = [(p[0], p[1]) for p in pts]

    def point_in_polygon_xy(pt, polygon):
        # ray-casting algorithm for point-in-polygon (polygon: list of (x,y))
        x, y = pt
        inside = False
        n = len(polygon)
        for i in range(n):
            xi, yi = polygon[i]
            xj, yj = polygon[(i + 1) % n]
            intersect = ((yi > y) != (yj > y)) and (
                x < (xj - xi) * (y - yi) / (yj - yi + 1e-30) + xi
            )
            if intersect:
                inside = not inside
        return inside

    def _on_segment(a, b, c):
        # check if point c is on segment ab (2D)
        (x1, y1), (x2, y2), (x3, y3) = (a, b, c)
        return (
            min(x1, x2) - 1e-12 <= x3 <= max(x1, x2) + 1e-12
            and min(y1, y2) - 1e-12 <= y3 <= max(y1, y2) + 1e-12
        )

    def _orient(a, b, c):
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

    def segments_intersect(p1, p2, q1, q2):
        # robust 2D segment intersection
        o1 = _orient(p1, p2, q1)
        o2 = _orient(p1, p2, q2)
        o3 = _orient(q1, q2, p1)
        o4 = _orient(q1, q2, p2)
        if o1 == 0 and _on_segment(p1, p2, q1):
            return True
        if o2 == 0 and _on_segment(p1, p2, q2):
            return True
        if o3 == 0 and _on_segment(q1, q2, p1):
            return True
        if o4 == 0 and _on_segment(q1, q2, p2):
            return True
        return (o1 * o2 < 0) and (o3 * o4 < 0)

    def triangle_edge_intersects_polygon(v_pts, polygon):
        # v_pts: list of three (x,y) tuples for triangle vertices
        tri_edges = [(v_pts[0], v_pts[1]), (v_pts[1], v_pts[2]), (v_pts[2], v_pts[0])]
        n = len(polygon)
        for a, b in tri_edges:
            for i in range(n):
                c = polygon[i]
                d = polygon[(i + 1) % n]
                if segments_intersect(a, b, c, d):
                    return True
        return False

    new_triangles_filtered = []
    for i, j, k in new_triangles:
        vi = new_vertices[i]
        vj = new_vertices[j]
        vk = new_vertices[k]
        cx = (vi[0] + vj[0] + vk[0]) / 3.0
        cy = (vi[1] + vj[1] + vk[1]) / 3.0
        # include triangle if centroid inside OR any vertex inside OR any edge intersects polygon
        centroid_in = point_in_polygon_xy((cx, cy), poly_xy)
        any_vertex_in = (
            point_in_polygon_xy((vi[0], vi[1]), poly_xy)
            or point_in_polygon_xy((vj[0], vj[1]), poly_xy)
            or point_in_polygon_xy((vk[0], vk[1]), poly_xy)
        )
        tri_xy = [(vi[0], vi[1]), (vj[0], vj[1]), (vk[0], vk[1])]
        edge_inter = triangle_edge_intersects_polygon(tri_xy, poly_xy)
        if centroid_in or any_vertex_in or edge_inter:
            new_triangles_filtered.append((i, j, k))

    # DEBUG: draw only the re-meshed (subdivided) triangles that include sample vertices
    # This highlights the small triangles created by inserting sample points (midpoints/samples).
    for i, j, k in new_triangles_filtered:
        # show triangle only if it contains at least one sampled vertex (appended samples)
        if not (i >= orig_vert_count or j >= orig_vert_count or k >= orig_vert_count):
            continue
        vi = new_vertices[i]
        vj = new_vertices[j]
        vk = new_vertices[k]
        pl = BRepBuilderAPI_MakePolygon()
        pl.Add(gp_Pnt(*vi))
        pl.Add(gp_Pnt(*vj))
        pl.Add(gp_Pnt(*vk))
        pl.Close()
        try:
            display.DisplayShape(pl.Wire(), update=False, color="DARKGRAY")
        except Exception:
            pass

    # DEBUG: highlight sample vertices (the appended seed vertices) in blue
    for idx in range(orig_vert_count, orig_vert_count + len(sample_coords)):
        sv = new_vertices[idx]
        vtx = BRepBuilderAPI_MakeVertex(gp_Pnt(*sv)).Vertex()
        display.DisplayShape(vtx, update=False, color="BLUE1")

    # ensure debug mesh is visible
    display.FitAll()

    # 4) seed indices are the indices of appended samples
    seed_indices = list(range(orig_vert_count, orig_vert_count + len(sample_coords)))

    if not seed_indices:
        print("No seeds found; skipping Voronoi.")
        display.DisplayShape(wire)

    labels, dists = geodesic_voronoi(new_vertices, new_triangles_filtered, seed_indices)
    # --- Visualize per-vertex Voronoi labels and boundary edges (helps verify partition) ---
    # create a color map for labels (cycle PALETTE)
    uniq = sorted(set(l for l in labels if l is not None and l >= 0))
    label_to_color = {lab: PALETTE[i % len(PALETTE)] for i, lab in enumerate(uniq)}

    # color vertices by label
    for vidx, v in enumerate(new_vertices):
        lab = labels[vidx]
        if lab is None or lab < 0:
            continue
        color = label_to_color.get(lab, PALETTE[0])
        vtx = BRepBuilderAPI_MakeVertex(gp_Pnt(*v)).Vertex()
        display.DisplayShape(vtx, update=False, color=color)

    # draw boundary edges where labels differ
    for i, j, k in new_triangles_filtered:
        for a, b in ((i, j), (j, k), (k, i)):
            if labels[a] != labels[b]:
                pa = new_vertices[a]
                pb = new_vertices[b]
                pl = BRepBuilderAPI_MakePolygon()
                pl.Add(gp_Pnt(*pa))
                pl.Add(gp_Pnt(*pb))
                pl.Close()
                display.DisplayShape(pl.Wire(), update=False, color="BLACK")

    # force redraw
    display.FitAll()

    # Diagnostic logging: seeds, mesh sizes, label stats
    uniq_labels = set([l for l in labels if l is not None and l >= 0])
    print(
        f"DEBUG: seeds={len(seed_indices)}, total_vertices={len(new_vertices)}, filtered_triangles={len(new_triangles_filtered)}"
    )
    print(
        f"DEBUG: unique_label_count={len(uniq_labels)} (showing up to 10): {list(uniq_labels)[:10]}"
    )
    # small sample of label assignments for first 20 vertices
    print("DEBUG: first 20 vertex labels:", labels[:20])

    # build region polygons by splitting triangles along bisectors so Voronoi regions
    # are not limited to whole triangles
    regions = build_voronoi_regions_by_splitting_triangles(
        new_vertices, new_triangles_filtered, labels
    )

    # Merge very small fragments to reduce slivers and missing-seed issues
    # compute mean triangle area for scale
    tri_areas = []
    for i, j, k in new_triangles_filtered:
        vi = np.array(new_vertices[i])
        vj = np.array(new_vertices[j])
        vk = np.array(new_vertices[k])
        ux = vj - vi
        vx = vk - vi
        cross = np.cross(ux, vx)
        tri_areas.append(0.5 * np.linalg.norm(cross))
    tri_area_mean = float(np.mean(tri_areas)) if tri_areas else 0.0
    regions = merge_small_fragments(
        regions, new_vertices, tri_area_mean, eps_rel=0.01, eps_abs=1e-8
    )

    # DEBUG: find seeds that ended up with no region pieces
    seeds_with_no_region = []
    for s in seed_indices:
        if s not in regions or len(regions.get(s, [])) == 0:
            seeds_with_no_region.append(s)
    print(f"DEBUG: seeds_with_no_region_count={len(seeds_with_no_region)}")
    if seeds_with_no_region:
        print(
            "DEBUG: seed indices with no region (sample up to 20):",
            seeds_with_no_region[:20],
        )
        # print coordinates
        coords = [new_vertices[s] for s in seeds_with_no_region[:50]]
        print("DEBUG: coords of first missing seeds:", coords[:10])
        # display them as red vertices for quick visual inspection
        for s in seeds_with_no_region:
            sv = new_vertices[s]
            vtx = BRepBuilderAPI_MakeVertex(gp_Pnt(*sv)).Vertex()
            display.DisplayShape(vtx, update=False, color="RED")

    # Diagnostic: region sizes
    region_sizes = {lab: len(polys) for lab, polys in regions.items()}
    total_regions = len(region_sizes)
    print(f"DEBUG: built regions count={total_regions}")
    # show a few region sizes
    items = list(region_sizes.items())
    items.sort(key=lambda x: -x[1])
    print("DEBUG: top 10 region sizes:", items[:10])

    # --- Validation checks for Voronoi correctness ---
    def _polygon_area(pts):
        # project polygon to best-fit plane and compute 2D polygon area (shoelace)
        if not pts or len(pts) < 3:
            return 0.0
        bf = _best_fit_plane(pts)
        if bf is None:
            return 0.0
        origin, u, v, normal, max_res = bf
        coords = []
        for p in pts:
            vec = np.array(p) - np.array(origin)
            x = float(np.dot(vec, np.array(u)))
            y = float(np.dot(vec, np.array(v)))
            coords.append((x, y))
        area = 0.0
        for i in range(len(coords)):
            x1, y1 = coords[i]
            x2, y2 = coords[(i + 1) % len(coords)]
            area += x1 * y2 - x2 * y1
        return abs(area) * 0.5

    # total area from original filtered triangles
    tri_area_sum = 0.0
    for i, j, k in new_triangles_filtered:
        vi = np.array(new_vertices[i])
        vj = np.array(new_vertices[j])
        vk = np.array(new_vertices[k])
        ux = vj - vi
        vx = vk - vi
        cross = np.cross(ux, vx)
        tri_area_sum += 0.5 * np.linalg.norm(cross)

    # total area from region polygons
    region_area_sum = 0.0
    per_label_area = {}
    for lab, polys in regions.items():
        a = 0.0
        for poly in polys:
            pts = normalize_polygon(poly)
            aa = _polygon_area(pts)
            a += aa
        per_label_area[lab] = a
        region_area_sum += a

    print(
        f"VALIDATION: triangle_area_sum={tri_area_sum:.6f}, region_area_sum={region_area_sum:.6f}, diff={abs(tri_area_sum-region_area_sum):.6f}"
    )

    # Check each seed is contained in at least one of its region polygons
    def _point_in_poly_2d(pt2, poly2):
        x, y = pt2
        inside = False
        n = len(poly2)
        for i in range(n):
            xi, yi = poly2[i]
            xj, yj = poly2[(i + 1) % n]
            intersect = ((yi > y) != (yj > y)) and (
                x < (xj - xi) * (y - yi) / (yj - yi + 1e-30) + xi
            )
            if intersect:
                inside = not inside
        return inside

    seeds_missing = []
    for s in seed_indices:
        found = False
        seed_pt = new_vertices[s]
        polys = regions.get(s, [])
        for poly in polys:
            pts = normalize_polygon(poly)
            if len(pts) < 3:
                continue
            bf = _best_fit_plane(pts)
            if bf is None:
                continue
            origin, u, v, normal, max_res = bf
            # project seed into plane coords
            vec = np.array(seed_pt) - np.array(origin)
            sx = float(np.dot(vec, np.array(u)))
            sy = float(np.dot(vec, np.array(v)))
            coords2 = []
            for p in pts:
                vecp = np.array(p) - np.array(origin)
                x = float(np.dot(vecp, np.array(u)))
                y = float(np.dot(vecp, np.array(v)))
                coords2.append((x, y))
            if _point_in_poly_2d((sx, sy), coords2):
                found = True
                break
        if not found:
            seeds_missing.append(s)

    print(f"VALIDATION: seeds_with_no_containing_polygon={len(seeds_missing)}")
    if seeds_missing:
        print("VALIDATION: sample missing seed indices:", seeds_missing[:20])

    # Detailed per-label diagnostics for first few labels
    print("DEBUG: per-label diagnostics (first 10 labels):")
    cnt = 0
    for lab in sorted(list(regions.keys()))[:10]:
        polys = regions.get(lab, [])
        a = per_label_area.get(lab, 0.0)
        print(f"  label={lab}, polys={len(polys)}, area={a:.6f}")
        for pi, poly in enumerate(polys[:3]):
            pts = normalize_polygon(poly)
            print(f"    poly#{pi} vertices={len(pts)}: {pts[:6]}")
        cnt += 1
    if cnt == 0:
        print("  (no regions to show)")

    # For display simplicity we only render the Voronoi wires (true Voronoi boundaries)
    # No compounds/faces are constructed here to keep rendering fast and faithful
    # to the computed Voronoi partition.
    print("INFO: displaying Voronoi wires only (no sewn triangles or filled faces)")

    # Report bisector failure diagnostics (edges where algebraic solve was unstable)
    bf_count = len(bisector_failures)
    print(f"DEBUG: bisector_failures_count= {bf_count}")
    if bf_count:
        print("DEBUG: bisector_failures (sample up to 10):", bisector_failures[:10])

    # deterministic color assignment per label
    labels_sorted = sorted(regions.keys())
    # color_map contains Quantity_Color (not always used by DisplayShape),
    # but for reliability we also make a string-based color map cycling PALETTE.
    color_map = {}
    color_map_str = {}
    nlab = max(1, len(labels_sorted))
    for i, lab in enumerate(labels_sorted):
        h = float(i) / nlab
        r, g, b = colorsys.hsv_to_rgb(h, 0.6, 0.95)
        color_map[lab] = Quantity_Color(r, g, b, Quantity_TOC_RGB)
        color_map_str[lab] = PALETTE[i % len(PALETTE)]

    # Display Voronoi boundaries as wires only (fast, minimal)
    for lab, poly_list in regions.items():
        colstr = color_map_str.get(lab, PALETTE[0])
        for poly in poly_list:
            pts = normalize_polygon(poly)
            if len(pts) < 2:
                continue
            pl = BRepBuilderAPI_MakePolygon()
            for p in pts:
                pl.Add(gp_Pnt(*p))
            pl.Close()
            display.DisplayShape(pl.Wire(), color="RED")

    # show outline and fit view (don't display the filled sewn face; it hides/occludes partition wires)
    display.DisplayShape(sewn, transparency=0.5)
    display.FitAll()
    start_display()
