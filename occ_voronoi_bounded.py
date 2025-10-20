#!/usr/bin/env python3
"""
Voronoi diagram for points restricted to [0,1]x[0,1] domain using PythonOCC-core
Handles boundary conditions properly by adding mirror points and clipping regions
"""

import numpy as np
from scipy.spatial import Voronoi
from OCC.Core.BRep import BRep_Builder
from OCC.Core.BRepBuilderAPI import (
    BRepBuilderAPI_MakeVertex,
    BRepBuilderAPI_MakeEdge,
    BRepBuilderAPI_MakeWire,
    BRepBuilderAPI_MakeFace,
    BRepBuilderAPI_MakePolygon,
)
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Ax3
from OCC.Core.TopoDS import TopoDS_Compound
from OCC.Core.TopTools import TopTools_ListOfShape
from OCC.Core.Quantity import (
    Quantity_NOC_PURPLE1,
    Quantity_NOC_WHEAT,
    Quantity_NOC_PINK,
)
from OCC.Display.SimpleGui import init_display
from OCC.Core.Geom import (
    Geom_Surface,
    Geom_SphericalSurface,
    Geom_ToroidalSurface,
    Geom_CylindricalSurface,
)
from OCC.Core.GCE2d import GCE2d_MakeSegment, GCE2d_MakeCircle
from OCC.Core.gp import gp_Pnt2d, gp_Circ2d, gp_Ax2d, gp_Dir2d
from OCC.Core.TColgp import TColgp_Array1OfPnt2d, TColgp_HArray1OfPnt2d
from OCC.Core.Geom2dAPI import Geom2dAPI_Interpolate


class VoronoiDomain:
    """
    Class to manage UV parameter domain bounds for Voronoi diagram with periodic boundary support
    """

    def __init__(
        self, umin=0.0, umax=1.0, vmin=0.0, vmax=1.0, periodic_u=False, periodic_v=False
    ):
        """
        Initialize UV domain bounds with optional periodic boundaries

        Args:
            umin, umax: U-parameter bounds
            vmin, vmax: V-parameter bounds
            periodic_u: If True, U-direction has periodic boundary (like longitude on sphere)
            periodic_v: If True, V-direction has periodic boundary (like latitude wrapping)
        """
        self.umin = float(umin)
        self.umax = float(umax)
        self.vmin = float(vmin)
        self.vmax = float(vmax)
        self.periodic_u = periodic_u
        self.periodic_v = periodic_v

        # Validate bounds
        if self.umin >= self.umax:
            raise ValueError("umin must be less than umax")
        if self.vmin >= self.vmax:
            raise ValueError("vmin must be less than vmax")

    def periodic_distance(self, p1, p2):
        """
        Calculate periodic distance between two points in the domain.
        """
        du = abs(p1[0] - p2[0])
        dv = abs(p1[1] - p2[1])
        if self.periodic_u:
            du = min(du, self.u_width - du)
        if self.periodic_v:
            dv = min(dv, self.v_height - dv)
        return np.sqrt(du ** 2 + dv ** 2)

    def pairwise_periodic_distances(points, domain):
        """
        Compute pairwise periodic distances for a set of points in the domain.
        Returns a (N, N) array of distances.
        """
        N = len(points)
        dists = np.zeros((N, N))
        for i in range(N):
            for j in range(N):
                dists[i, j] = domain.periodic_distance(points[i], points[j])
        return dists

    @property
    def u_width(self):
        """U-parameter domain width"""
        return self.umax - self.umin

    @property
    def v_height(self):
        """V-parameter domain height"""
        return self.vmax - self.vmin

    @property
    def center(self):
        """Domain center point in UV coordinates"""
        return [(self.umin + self.umax) / 2, (self.vmin + self.vmax) / 2]

    @property
    def bounds(self):
        """Return bounds as tuple (umin, umax, vmin, vmax)"""
        return (self.umin, self.umax, self.vmin, self.vmax)

    def contains_point(self, point):
        """Check if point is within UV domain"""
        u, v = point
        return self.umin <= u <= self.umax and self.vmin <= v <= self.vmax

    def wrap_point(self, point):
        """Wrap UV coordinates for periodic boundaries"""
        u, v = point
        if self.periodic_u:
            u = self.umin + (u - self.umin) % self.u_width
        if self.periodic_v:
            v = self.vmin + (v - self.vmin) % self.v_height
        return [u, v]

    def get_periodic_copies(self, point, layers=1):
        """Get periodic copies of a UV point for boundary handling"""
        u, v = point
        copies = []

        # Generate offsets based on periodic directions
        u_offsets = []
        v_offsets = []

        if self.periodic_u:
            for i in range(-layers, layers + 1):
                u_offsets.append(i * self.u_width)
        else:
            u_offsets = [0]

        if self.periodic_v:
            for i in range(-layers, layers + 1):
                v_offsets.append(i * self.v_height)
        else:
            v_offsets = [0]

        # Generate all combinations
        for du in u_offsets:
            for dv in v_offsets:
                if du != 0 or dv != 0:  # Exclude original point
                    copies.append([u + du, v + dv])

        return copies

    def __str__(self):
        periodic_str = ""
        if self.periodic_u or self.periodic_v:
            periodic_flags = []
            if self.periodic_u:
                periodic_flags.append("periodic_u")
            if self.periodic_v:
                periodic_flags.append("periodic_v")
            periodic_str = f" ({', '.join(periodic_flags)})"
        return f"VoronoiDomain([{self.umin:.2f},{self.umax:.2f}] × [{self.vmin:.2f},{self.vmax:.2f}]){periodic_str}"


def generate_points_in_domain(num_points, domain, seed=None):
    """
    Generate random points within specified UV domain bounds

    Args:
        num_points: Number of points to generate
        domain: VoronoiDomain object specifying UV bounds
        seed: Random seed for reproducibility

    Returns:
        numpy array of shape (num_points, 2) with UV coordinates within domain bounds
    """
    if seed is not None:
        np.random.seed(seed)

    # Generate points in [0,1]×[0,1] then scale and translate to UV domain bounds
    points_01 = np.random.rand(num_points, 2)

    # Scale and translate to UV domain bounds
    points = np.zeros_like(points_01)
    points[:, 0] = domain.umin + points_01[:, 0] * domain.u_width  # U coordinate
    points[:, 1] = domain.vmin + points_01[:, 1] * domain.v_height  # V coordinate

    return points


def add_boundary_points(points, domain, boundary_density=20):
    """
    Add virtual boundary points to handle edge effects in Voronoi diagram
    For periodic domains, adds periodic copies of original points instead of artificial boundary points

    Args:
        points: Original points within UV domain bounds
        domain: VoronoiDomain object specifying UV bounds
        boundary_density: Number of points per boundary edge (for non-periodic boundaries)

    Returns:
        Extended point set including boundary points or periodic copies
    """
    boundary_points = []

    if domain.periodic_u or domain.periodic_v:
        # For periodic domains, add periodic copies of original points
        for point in points:
            periodic_copies = domain.get_periodic_copies(point, layers=1)
            boundary_points.extend(periodic_copies)
    else:
        # For non-periodic domains, use artificial boundary points
        # Calculate margin based on UV domain size
        margin_u = domain.u_width * 0.5  # 50% margin
        margin_v = domain.v_height * 0.5

        # Extended bounds for virtual boundary points
        extended_umin = domain.umin - margin_u
        extended_umax = domain.umax + margin_u
        extended_vmin = domain.vmin - margin_v
        extended_vmax = domain.vmax + margin_v

        # Bottom edge
        u_vals = np.linspace(extended_umin, extended_umax, boundary_density)
        for u in u_vals:
            boundary_points.append([u, extended_vmin])

        # Top edge
        for u in u_vals:
            boundary_points.append([u, extended_vmax])

        # Left edge
        v_vals = np.linspace(extended_vmin, extended_vmax, boundary_density)
        for v in v_vals:
            boundary_points.append([extended_umin, v])

        # Right edge
        for v in v_vals:
            boundary_points.append([extended_umax, v])

        # Corner points for better boundary handling
        corners = [
            [extended_umin, extended_vmin],
            [extended_umin, extended_vmax],
            [extended_umax, extended_vmin],
            [extended_umax, extended_vmax],
        ]
        boundary_points.extend(corners)

    # Combine original points with boundary points
    all_points = np.vstack([points, np.array(boundary_points)])

    return all_points, len(points)  # Return extended points and original count


def clip_line_to_domain(p1, p2, domain):
    """
    Clip a line segment to the domain bounds using Cohen-Sutherland algorithm

    Args:
        p1, p2: Line endpoints as [u, v] coordinates
        domain: VoronoiDomain object specifying UV bounds

    Returns:
        Tuple of clipped endpoints (p1_clipped, p2_clipped) or None if line is outside
    """

    # Cohen-Sutherland outcodes
    def compute_outcode(u, v):
        code = 0
        if u < domain.umin:  # Left
            code |= 1
        elif u > domain.umax:  # Right
            code |= 2
        if v < domain.vmin:  # Bottom
            code |= 4
        elif v > domain.vmax:  # Top
            code |= 8
        return code

    u1, v1 = p1
    u2, v2 = p2

    outcode1 = compute_outcode(u1, v1)
    outcode2 = compute_outcode(u2, v2)

    while True:
        # Both points inside
        if (outcode1 | outcode2) == 0:
            return ([u1, v1], [u2, v2])

        # Both points outside same region
        if (outcode1 & outcode2) != 0:
            return None

        # At least one point outside, find intersection
        if outcode1 != 0:
            outcode_out = outcode1
        else:
            outcode_out = outcode2

        # Find intersection point
        if outcode_out & 8:  # Top
            u = u1 + (u2 - u1) * (domain.vmax - v1) / (v2 - v1)
            v = domain.vmax
        elif outcode_out & 4:  # Bottom
            u = u1 + (u2 - u1) * (domain.vmin - v1) / (v2 - v1)
            v = domain.vmin
        elif outcode_out & 2:  # Right
            v = v1 + (v2 - v1) * (domain.umax - u1) / (u2 - u1)
            u = domain.umax
        elif outcode_out & 1:  # Left
            v = v1 + (v2 - v1) * (domain.umin - u1) / (u2 - u1)
            u = domain.umin

        # Replace point outside with intersection
        if outcode_out == outcode1:
            u1, v1 = u, v
            outcode1 = compute_outcode(u1, v1)
        else:
            u2, v2 = u, v
            outcode2 = compute_outcode(u2, v2)


def clip_polygon_to_standard_domain(vertices, domain):
    """Standard Sutherland-Hodgman clipping for non-periodic domains"""

    def is_inside_periodic(point, edge):
        u, v = point
        if edge == "left":
            if domain.periodic_u:
                return True
            return u >= domain.umin
        elif edge == "right":
            if domain.periodic_u:
                return True
            return u <= domain.umax
        elif edge == "bottom":
            if domain.periodic_v:
                return True
            return v >= domain.vmin
        elif edge == "top":
            if domain.periodic_v:
                return True
            return v <= domain.vmax
        return True

    def compute_intersection_periodic(p1, p2, edge):
        u1, v1 = p1
        u2, v2 = p2
        if edge == "left" and not domain.periodic_u and u1 != u2:
            v = v1 + (v2 - v1) * (domain.umin - u1) / (u2 - u1)
            return [domain.umin, v]
        elif edge == "right" and not domain.periodic_u and u1 != u2:
            v = v1 + (v2 - v1) * (domain.umax - u1) / (u2 - u1)
            return [domain.umax, v]
        elif edge == "bottom" and not domain.periodic_v and v1 != v2:
            u = u1 + (u2 - u1) * (domain.vmin - v1) / (v2 - v1)
            return [u, domain.vmin]
        elif edge == "top" and not domain.periodic_v and v1 != v2:
            u = u1 + (u2 - u1) * (domain.vmax - v1) / (v2 - v1)
            return [u, domain.vmax]
        return None

    clipped = vertices[:]
    for edge in ["left", "right", "bottom", "top"]:
        if not clipped:
            break
        input_list = clipped[:]
        clipped = []
        if input_list:
            s = input_list[-1]
            for e in input_list:
                if is_inside_periodic(e, edge):
                    if not is_inside_periodic(s, edge):
                        intersection = compute_intersection_periodic(s, e, edge)
                        if intersection is not None:
                            clipped.append(intersection)
                    clipped.append(e)
                elif is_inside_periodic(s, edge):
                    intersection = compute_intersection_periodic(s, e, edge)
                    if intersection is not None:
                        clipped.append(intersection)
                s = e
    return clipped


def generate_bounded_voronoi_diagram(points, domain, boundary_density=20):
    """
    Generate Voronoi diagram bounded to specified domain

    Args:
        points: Points within domain bounds
        domain: VoronoiDomain object specifying bounds
        boundary_density: Density of boundary points for proper edge handling

    Returns:
        List of Voronoi regions (each region is list of vertices)
    """

    # 周期コピー点群を生成（分割計算には主周期点＋周期コピー点を使う。region assignmentは主周期点のみ）
    points_main = np.array(points)
    periodic_points = [points_main]
    if domain.periodic_u:
        u_period = domain.u_width
        periodic_points += [points_main + np.array([u_period, 0]), points_main - np.array([u_period, 0])]
    if domain.periodic_v:
        v_period = domain.v_height
        periodic_points += [points_main + np.array([0, v_period]), points_main - np.array([0, v_period])]
    all_points = np.vstack(periodic_points)

    # 境界点追加は非周期方向のみ
    if domain.periodic_u or domain.periodic_v:
        extended_points = all_points
        original_count = len(points_main)
    else:
        extended_points, original_count = add_boundary_points(all_points, domain, boundary_density)

    # Voronoi分割（scipyは距離計算を周期的にできないので、region assignmentを周期距離で再計算する）
    vor = Voronoi(extended_points)

    # region assignment: 主周期点のみ（まず）
    regions = []
    region_point_indices = []
    for i in range(original_count):
        region_idx = vor.point_region[i]
        region = vor.regions[region_idx]
        if not region or -1 in region:
            continue
        vertices = [vor.vertices[v] for v in region if v >= 0]
        if len(vertices) >= 3:
            clipped = clip_polygon_to_standard_domain(vertices, domain)
            if len(clipped) >= 3:
                regions.append(clipped)
                region_point_indices.append(i)

    # 空白領域検出と追加点選定
    # 1. Domainを細かいグリッドでサンプリング
    grid_density = 100
    u_grid = np.linspace(domain.umin, domain.umax, grid_density)
    v_grid = np.linspace(domain.vmin, domain.vmax, grid_density)
    grid_points = np.array([[u, v] for u in u_grid for v in v_grid])

    # 2. 各グリッド点がどのregionに含まれるか判定
    covered = np.zeros(len(grid_points), dtype=bool)
    from matplotlib.path import Path
    for region in regions:
        path = Path(region)
        covered |= path.contains_points(grid_points)

    # 3. 未カバー点があれば、その近傍の周期コピー点を追加
    if not np.all(covered):
        # 周期コピー点のインデックス
        periodic_indices = list(range(original_count, len(all_points)))
        # 追加候補点
        candidate_points = all_points[periodic_indices]
        # 追加点インデックス
        add_indices = set()
        for idx, pt in enumerate(candidate_points):
            # その点のregionが空白領域をカバーするか判定
            region_idx = vor.point_region[periodic_indices[idx]]
            region = vor.regions[region_idx]
            if not region or -1 in region:
                continue
            vertices = [vor.vertices[v] for v in region if v >= 0]
            if len(vertices) >= 3:
                clipped = clip_polygon_to_standard_domain(vertices, domain)
                if len(clipped) >= 3:
                    path = Path(clipped)
                    if np.any(path.contains_points(grid_points[~covered])):
                        add_indices.add(periodic_indices[idx])
        # 追加点のregionもregionsに加える
        for idx in add_indices:
            region_idx = vor.point_region[idx]
            region = vor.regions[region_idx]
            if not region or -1 in region:
                continue
            vertices = [vor.vertices[v] for v in region if v >= 0]
            if len(vertices) >= 3:
                clipped = clip_polygon_to_standard_domain(vertices, domain)
                if len(clipped) >= 3:
                    regions.append(clipped)
    return regions, vor


def create_voronoi_faces(regions):
    """
    Create OpenCASCADE faces from Voronoi regions

    Args:
        regions: List of Voronoi regions (each is list of vertices)

    Returns:
        List of TopoDS_Face objects
    """
    faces = []

    for region in regions:
        if len(region) >= 3:
            try:
                # Create polygon from vertices
                polygon_builder = BRepBuilderAPI_MakePolygon()

                for vertex in region:
                    pnt = gp_Pnt(float(vertex[0]), float(vertex[1]), 0.0)
                    polygon_builder.Add(pnt)

                polygon_builder.Close()

                if polygon_builder.IsDone():
                    wire = polygon_builder.Wire()

                    # Create face from wire
                    face_builder = BRepBuilderAPI_MakeFace(wire)
                    if face_builder.IsDone():
                        faces.append(face_builder.Face())

            except Exception as e:
                print(f"Error creating face from region: {e}")
                continue

    return faces


def create_domain_boundary(domain):
    """
    Create the boundary of specified domain

    Args:
        domain: VoronoiDomain object specifying bounds

    Returns:
        TopoDS_Face representing the domain rectangle
    """
    # Create corner points
    p1 = gp_Pnt(domain.umin, domain.vmin, 0.0)
    p2 = gp_Pnt(domain.umax, domain.vmin, 0.0)
    p3 = gp_Pnt(domain.umax, domain.vmax, 0.0)
    p4 = gp_Pnt(domain.umin, domain.vmax, 0.0)

    # Create boundary polygon
    polygon = BRepBuilderAPI_MakePolygon()
    polygon.Add(p1)
    polygon.Add(p2)
    polygon.Add(p3)
    polygon.Add(p4)
    polygon.Close()

    # Create face
    face_builder = BRepBuilderAPI_MakeFace(polygon.Wire())
    return face_builder.Face()


def get_clipped_voronoi_edges(vor, domain, original_count):
    """
    Get Voronoi edges that are properly clipped to domain bounds

    Args:
        vor: Voronoi diagram object
        domain: VoronoiDomain object specifying bounds
        original_count: Number of original points (excluding boundary points)

    Returns:
        List of edge segments within domain bounds
    """
    clipped_edges = []

    for ridge_idx, ridge in enumerate(vor.ridge_vertices):
        if -1 not in ridge:  # Only finite edges
            ridge_points = vor.ridge_points[ridge_idx]

            # Include edges that involve at least one original point
            if ridge_points[0] < original_count or ridge_points[1] < original_count:
                v1 = vor.vertices[ridge[0]]
                v2 = vor.vertices[ridge[1]]

                # Clip edge to domain bounds
                clipped_edge = clip_line_to_domain(
                    [v1[0], v1[1]], [v2[0], v2[1]], domain
                )

                if clipped_edge is not None:
                    # Only add edges that have significant length within unit square
                    p1, p2 = clipped_edge
                    edge_length = np.sqrt((p2[0] - p1[0]) ** 2 + (p2[1] - p1[1]) ** 2)
                    if edge_length > 1e-6:  # Filter out very short segments
                        clipped_edges.append(clipped_edge)

    return clipped_edges


def create_periodic_voronoi_edges(regions, domain):
    """
    Create periodic connections for Voronoi edges across boundaries

    Args:
        regions: List of Voronoi regions
        domain: VoronoiDomain with periodic flags

    Returns:
        List of additional edge segments connecting across periodic boundaries
    """
    periodic_edges = []

    if not (domain.periodic_u or domain.periodic_v):
        return periodic_edges

    # For each region, check if it touches a periodic boundary
    for region in regions:
        if len(region) < 3:
            continue

        # Check for vertices at domain boundaries that need periodic connection
        boundary_vertices = []

        for vertex in region:
            u, v = vertex
            at_boundary = False

            # Check if vertex is at a periodic boundary
            if domain.periodic_u:
                if abs(u - domain.umin) < 1e-6 or abs(u - domain.umax) < 1e-6:
                    at_boundary = True

            if domain.periodic_v:
                if abs(v - domain.vmin) < 1e-6 or abs(v - domain.vmax) < 1e-6:
                    at_boundary = True

            if at_boundary:
                boundary_vertices.append(vertex)

        # Create periodic connections between boundary vertices
        if len(boundary_vertices) >= 2:
            for i in range(len(boundary_vertices)):
                v1 = boundary_vertices[i]
                v2 = boundary_vertices[(i + 1) % len(boundary_vertices)]

                # Create wrapped versions of vertices for connection
                v1_wrapped = domain.wrap_point(v1)
                v2_wrapped = domain.wrap_point(v2)

                # Add edge if it crosses a periodic boundary
                if (
                    abs(v1[0] - v1_wrapped[0]) > 1e-6
                    or abs(v1[1] - v1_wrapped[1]) > 1e-6
                    or abs(v2[0] - v2_wrapped[0]) > 1e-6
                    or abs(v2[1] - v2_wrapped[1]) > 1e-6
                ):
                    periodic_edges.append([v1_wrapped, v2_wrapped])

    return periodic_edges


def create_voronoi_on_parametric_surface(
    uv_points, surface, domain, boundary_density=20
):
    """
    Create Voronoi diagram on a parametric surface using UV coordinates

    Args:
        uv_points: List of [u, v] coordinates in parameter space
        surface: Geom_Surface object (Sphere, Torus, etc.)
        domain: VoronoiDomain defining UV parameter bounds with periodic flags
        boundary_density: Density for boundary point generation

    Returns:
        List of TopoDS_Wire objects representing Voronoi cells on surface
    """
    # Generate Voronoi in UV parameter space
    regions, vor = generate_bounded_voronoi_diagram(uv_points, domain, boundary_density)

    # Convert each region to a wire on the parametric surface
    surface_wires = []
    seen_centroids = set()

    for region in regions:
        if len(region) < 3:
            continue

        # 領域の重複除去: セントロイドが主周期内にあるものだけ
        centroid = np.mean(region, axis=0)
        u, v = centroid
        # wrapして主周期内に戻す
        u_wrapped = domain.umin + (u - domain.umin) % domain.u_width
        v_wrapped = domain.vmin + (v - domain.vmin) % domain.v_height
        key = (round(u_wrapped, 8), round(v_wrapped, 8))
        if key in seen_centroids:
            continue
        seen_centroids.add(key)
        # 主周期内のみ
        if not (domain.umin <= u_wrapped <= domain.umax and domain.vmin <= v_wrapped <= domain.vmax):
            continue

        try:
            # Create edges in parameter space
            edges = []
            n = len(region)

            for i in range(n):
                p1_uv = region[i]
                p2_uv = region[(i + 1) % n]

                # Handle periodic wrapping for UV coordinates
                if domain.periodic_u or domain.periodic_v:
                    p1_uv = domain.wrap_point(p1_uv)
                    p2_uv = domain.wrap_point(p2_uv)

                # Create 2D segment in UV space
                p1_2d = gp_Pnt2d(p1_uv[0], p1_uv[1])
                p2_2d = gp_Pnt2d(p2_uv[0], p2_uv[1])

                # Handle periodic boundary crossing
                if domain.periodic_u:
                    du = p2_uv[0] - p1_uv[0]
                    if abs(du) > domain.u_width / 2:
                        # Edge crosses periodic boundary, need special handling
                        if du > 0:
                            p2_2d = gp_Pnt2d(p2_uv[0] - domain.u_width, p2_uv[1])
                        else:
                            p2_2d = gp_Pnt2d(p2_uv[0] + domain.u_width, p2_uv[1])

                if domain.periodic_v:
                    dv = p2_uv[1] - p1_uv[1]
                    if abs(dv) > domain.v_height / 2:
                        # Edge crosses periodic boundary
                        if dv > 0:
                            p2_2d = gp_Pnt2d(p2_uv[0], p2_uv[1] - domain.v_height)
                        else:
                            p2_2d = gp_Pnt2d(p2_uv[0], p2_uv[1] + domain.v_height)

                # Create 2D segment
                seg2d = GCE2d_MakeSegment(p1_2d, p2_2d).Value()

                # Map to 3D surface
                edge = BRepBuilderAPI_MakeEdge(seg2d, surface).Edge()
                edges.append(edge)

            # Create wire from edges
            wire_maker = BRepBuilderAPI_MakeWire()
            for edge in edges:
                wire_maker.Add(edge)

            if wire_maker.IsDone():
                surface_wires.append(wire_maker.Wire())

        except Exception as e:
            print(f"Error creating surface wire for region: {e}")
            continue

    return surface_wires


def get_boundary_edges(domain):
    """
    Get the boundary edges of the specified domain

    Args:
        domain: VoronoiDomain object specifying bounds

    Returns:
        List of boundary edge segments
    """
    boundary_edges = [
        ([domain.umin, domain.vmin], [domain.umax, domain.vmin]),  # Bottom edge
        ([domain.umax, domain.vmin], [domain.umax, domain.vmax]),  # Right edge
        ([domain.umax, domain.vmax], [domain.umin, domain.vmax]),  # Top edge
        ([domain.umin, domain.vmax], [domain.umin, domain.vmin]),  # Left edge
    ]
    return boundary_edges


def visualize_bounded_voronoi(
    points, regions, domain, vor=None, surface=None, display_options=None
):
    """
    Visualize the bounded Voronoi diagram using OpenCASCADE

    Args:
        points: Original points within domain bounds
        regions: Voronoi regions to display
        domain: VoronoiDomain object specifying bounds
        vor: Voronoi object (optional, for displaying edges)
        surface: Geom_Surface object for parametric surface display (optional)
        display_options: Dictionary of display options
    """
    if display_options is None:
        display_options = {
            "show_points": True,
            "show_regions": True,
            "show_edges": True,
            "show_boundary": True,
            "region_transparency": 0.7,
        }

    # Initialize display
    display, start_display, add_menu, add_function_to_menu = init_display()

    # Display domain boundary
    if display_options.get("show_boundary", True):
        boundary_face = create_domain_boundary(domain)
        display.DisplayShape(
            boundary_face, color="black", transparency=0.9, update=False
        )

    # Display Voronoi regions as faces
    if display_options.get("show_regions", True):
        # 主周期内のセル
        faces = create_voronoi_faces(regions)
        colors = ["red", "green", "blue1", "yellow", "cyan1", "orange"]
        for i, face in enumerate(faces):
            color = colors[i % len(colors)]
            transparency = display_options.get("region_transparency", 0.7)
            display.DisplayShape(
                face, color=color, transparency=transparency, update=False
            )

        # 周期方向に分割を延長
        periodic_shifts = []
        if domain.periodic_u:
            periodic_shifts += [[domain.u_width, 0], [-domain.u_width, 0]]
        if domain.periodic_v:
            periodic_shifts += [[0, domain.v_height], [0, -domain.v_height]]
        # 複数周期方向の場合は全組み合わせ
        if domain.periodic_u and domain.periodic_v:
            periodic_shifts += [
                [domain.u_width, domain.v_height],
                [domain.u_width, -domain.v_height],
                [-domain.u_width, domain.v_height],
                [-domain.u_width, -domain.v_height],
            ]

        for shift in periodic_shifts:
            shifted_regions = []
            for region in regions:
                shifted = [[pt[0] + shift[0], pt[1] + shift[1]] for pt in region]
                shifted_regions.append(shifted)
            shifted_faces = create_voronoi_faces(shifted_regions)
            for i, face in enumerate(shifted_faces):
                color = colors[i % len(colors)]
                transparency = display_options.get("region_transparency", 0.7)
                display.DisplayShape(
                    face, color=color, transparency=transparency, update=False
                )

    # Display parametric surface if provided
    if surface is not None:
        # Create surface face for display
        if isinstance(surface, Geom_SphericalSurface):
            surface_face = BRepBuilderAPI_MakeFace(
                surface.Sphere(), 0, 2 * np.pi, -np.pi / 2, np.pi / 2
            ).Face()
        elif isinstance(surface, Geom_ToroidalSurface):
            surface_face = BRepBuilderAPI_MakeFace(
                surface.Torus(), 0, 2 * np.pi, 0, 2 * np.pi
            ).Face()
        elif isinstance(surface, Geom_CylindricalSurface):
            surface_face = BRepBuilderAPI_MakeFace(
                surface.Cylinder(), 0, 2 * np.pi, domain.ymin, domain.ymax
            ).Face()
        else:
            surface_face = None

        if surface_face is not None:
            display.DisplayShape(
                surface_face, transparency=0.8, update=False
            )

        # Create and display Voronoi wires on surface
        if display_options.get("show_regions", True):
            surface_wires = create_voronoi_on_parametric_surface(
                points, surface, domain
            )
            print(f"Displaying {len(surface_wires)} Voronoi wires on surface")

            colors = ["red", "green", "blue1", "yellow", "cyan1", "orange"]
            for i, wire in enumerate(surface_wires):
                color = colors[i % len(colors)]
                display.DisplayShape(wire, color=color, update=False)

    # Display Voronoi edges (properly clipped to domain bounds)
    if display_options.get("show_edges", True) and vor is not None and surface is None:
        clipped_edges = get_clipped_voronoi_edges(vor, domain, len(points))

        print(f"Displaying {len(clipped_edges)} clipped Voronoi edges")
        for edge_segment in clipped_edges:
            p1_clipped, p2_clipped = edge_segment
            pnt1 = gp_Pnt(float(p1_clipped[0]), float(p1_clipped[1]), 0.0)
            pnt2 = gp_Pnt(float(p2_clipped[0]), float(p2_clipped[1]), 0.0)
            edge = BRepBuilderAPI_MakeEdge(pnt1, pnt2).Edge()
            display.DisplayShape(edge, color="blue1", update=False)

        # Display periodic connections if domain has periodic boundaries
        if domain.periodic_u or domain.periodic_v:
            periodic_edges = create_periodic_voronoi_edges(regions, domain)
            print(f"Displaying {len(periodic_edges)} periodic connection edges")

            for edge_segment in periodic_edges:
                p1, p2 = edge_segment
                pnt1 = gp_Pnt(float(p1[0]), float(p1[1]), 0.01)
                pnt2 = gp_Pnt(float(p2[0]), float(p2[1]), 0.01)
                edge = BRepBuilderAPI_MakeEdge(pnt1, pnt2).Edge()
                display.DisplayShape(edge, color="purple", update=False)

    # Display boundary edges more prominently
    if display_options.get("show_boundary", True):
        boundary_edges = get_boundary_edges(domain)

        for edge_segment in boundary_edges:
            p1, p2 = edge_segment
            pnt1 = gp_Pnt(float(p1[0]), float(p1[1]), 0.01)  # Slightly elevated
            pnt2 = gp_Pnt(float(p2[0]), float(p2[1]), 0.01)
            edge = BRepBuilderAPI_MakeEdge(pnt1, pnt2).Edge()
            display.DisplayShape(edge, color="black", update=False)

    # Display original points and their periodic copies
    if display_options.get("show_points", True):
        # 主周期内
        for i, point in enumerate(points):
            pnt = gp_Pnt(float(point[0]), float(point[1]), 0.02)
            vertex_shape = BRepBuilderAPI_MakeVertex(pnt).Shape()
            display.DisplayShape(vertex_shape, color="red", update=False)

        # 周期コピー
        periodic_shifts = []
        if domain.periodic_u:
            periodic_shifts += [[domain.u_width, 0], [-domain.u_width, 0]]
        if domain.periodic_v:
            periodic_shifts += [[0, domain.v_height], [0, -domain.v_height]]
        if domain.periodic_u and domain.periodic_v:
            periodic_shifts += [
                [domain.u_width, domain.v_height],
                [domain.u_width, -domain.v_height],
                [-domain.u_width, domain.v_height],
                [-domain.u_width, -domain.v_height],
            ]
        for shift in periodic_shifts:
            for point in points:
                pnt = gp_Pnt(
                    float(point[0] + shift[0]), float(point[1] + shift[1]), 0.02
                )
                vertex_shape = BRepBuilderAPI_MakeVertex(pnt).Shape()
                display.DisplayShape(vertex_shape, color="red", update=False)

    # Fit view and start display
    display.FitAll()
    print("Starting Voronoi diagram visualization...")
    print("Controls:")
    print("- Mouse wheel: Zoom in/out")
    print("- Left mouse button: Rotate")
    print("- Middle mouse button: Pan")
    print("- Press 'q' or close window to exit")

    start_display()


if __name__ == "__main__":
    """
    Main function demonstrating bounded Voronoi diagram generation with periodic boundaries and surfaces
    """
    print("Bounded Voronoi Diagram with Periodic Boundaries and Parametric Surfaces")
    print("=" * 80)

    num_points = 30
    seed = 41

    # Choose example to run
    example_choice = 5  # Change this to run different examples
    boundary_density = 25

    # Example selection and domain setup
    if example_choice == 1:
        print("\n1. Standard non-periodic domain [-2,3]x[-1,2]:")
        domain = VoronoiDomain(-2.0, 3.0, -1.0, 2.0, periodic_u=False, periodic_v=False)
    elif example_choice == 2:
        print("\n2. Periodic domain [0, 2π]x[-1, 1] with U-periodic:")
        domain = VoronoiDomain(
            0.0, 2 * np.pi, -1.0, 1.0, periodic_u=True, periodic_v=False
        )
    elif example_choice == 3:
        print("\n3. Periodic domain [0, 2π]x[-1, 1] with V-periodic:")
        domain = VoronoiDomain(
            0.0, 2 * np.pi, -1.0, 1.0, periodic_u=False, periodic_v=True
        )
    elif example_choice == 4:
        print("\n4. Spherical surface [0, 2π]x[0, π] with periodic longitude:")
        domain = VoronoiDomain(
            0.0, 2 * np.pi, 0.0, np.pi, periodic_u=True, periodic_v=False
        )
    elif example_choice == 5:
        print("\n5. Toroidal surface [0, 2π]x[0, 2π] with both directions periodic:")
        domain = VoronoiDomain(
            0.0, 2 * np.pi, 0.0, 2 * np.pi, periodic_u=True, periodic_v=True
        )
    else:
        example_choice = 1
        domain = VoronoiDomain(-2.0, 3.0, -1.0, 2.0, periodic_u=False, periodic_v=False)

    print(f"   Domain: {domain}")

    # Generate points in domain
    points = generate_points_in_domain(num_points, domain, seed=seed)

    # Voronoi calculation（主周期内の点群のみ渡す）
    regions, vor = generate_bounded_voronoi_diagram(
        points, domain, boundary_density=boundary_density
    )

    # Surface setup
    surface = None
    if example_choice == 4:
        sphere_surface = Geom_SphericalSurface(gp_Ax3(), 10.0)
        surface = sphere_surface
    elif example_choice == 6:
        torus_surface = Geom_ToroidalSurface(gp_Ax3(), 15.0, 5.0)
        surface = torus_surface

    print(f"   Generated {num_points} random points")
    print(f"   Sample points: {points[:3]}")
    print(f"   Created {len(regions)} bounded Voronoi regions")

    # Display options
    display_options = {
        "show_points": True,
        "show_regions": True,
        "show_edges": True,
        "show_boundary": True,
        "region_transparency": 0.6,
    }

    # Visualize result
    print(f"\nLaunching visualization for {domain}")
    if surface is not None:
        print("Note: Voronoi diagram projected onto parametric surface")
    if domain.periodic_u or domain.periodic_v:
        print("Note: Periodic boundaries enabled - diagram wraps across boundaries")

    visualize_bounded_voronoi(points, regions, domain, vor, surface, display_options)
