import numpy as np
from scipy.spatial import Voronoi, ConvexHull
from OCC.Core.BRep import BRep_Builder, BRep_Tool
from OCC.Core.BRepBuilderAPI import (
    BRepBuilderAPI_MakeVertex,
    BRepBuilderAPI_MakeEdge,
    BRepBuilderAPI_MakeWire,
    BRepBuilderAPI_MakeFace,
)
from OCC.Core.gp import gp_Pnt
from OCC.Core.TopoDS import TopoDS_Compound, TopoDS_Edge, TopoDS_Face
from OCC.Core.TopExp import TopExp_Explorer, topexp
from OCC.Core.TopAbs import TopAbs_EDGE, TopAbs_FACE
from OCC.Core.TopTools import TopTools_IndexedDataMapOfShapeListOfShape
from OCC.Display.SimpleGui import init_display


def generate_points_on_plane(num_points):
    """
    Function to generate random points on a plane

    Args:
        num_points: Number of points to generate. Example: 30

    Returns:
        Generated points. Shape is (num_points, 3).
    """
    points = np.random.rand(num_points, 2) * 100
    # Place points within a 100x100 range
    return points


def clip_voronoi_to_boundary(vor, boundary_box):
    """
    Clip Voronoi diagram to a rectangular boundary

    Args:
        vor: Voronoi diagram object
        boundary_box: [xmin, xmax, ymin, ymax]

    Returns:
        List of clipped regions (each region is a list of vertices)
    """
    xmin, xmax, ymin, ymax = boundary_box

    clipped_regions = []

    for point_idx, region_idx in enumerate(vor.point_region):
        region = vor.regions[region_idx]

        if not region or -1 in region:
            # Handle infinite regions by clipping to boundary
            vertices = clip_infinite_region(vor, point_idx, boundary_box)
        else:
            # Finite region - just get vertices
            vertices = [vor.vertices[i] for i in region]

        if len(vertices) >= 3:  # Valid polygon
            # Ensure vertices are within boundary
            clipped_vertices = clip_polygon_to_box(vertices, boundary_box)
            if len(clipped_vertices) >= 3:
                clipped_regions.append(clipped_vertices)

    return clipped_regions


def clip_infinite_region(vor, point_idx, boundary_box):
    """
    Clip an infinite Voronoi region to boundary box
    """
    xmin, xmax, ymin, ymax = boundary_box
    center = vor.points[point_idx]
    region_idx = vor.point_region[point_idx]
    region = vor.regions[region_idx]

    if not region:
        return []

    # Get finite vertices
    finite_vertices = []
    for vertex_idx in region:
        if vertex_idx >= 0:
            finite_vertices.append(vor.vertices[vertex_idx])

    if len(finite_vertices) < 2:
        # Create a simple boundary region around the point
        margin = 5.0
        return [
            [max(xmin, center[0] - margin), max(ymin, center[1] - margin)],
            [min(xmax, center[0] + margin), max(ymin, center[1] - margin)],
            [min(xmax, center[0] + margin), min(ymax, center[1] + margin)],
            [max(xmin, center[0] - margin), min(ymax, center[1] + margin)],
        ]

    # For infinite regions, we need to extend rays to boundary
    vertices = []

    # Add finite vertices
    for vertex in finite_vertices:
        vertices.append(vertex)

    # Add boundary intersections for infinite rays
    boundary_points = []

    # Find intersections with boundary
    for i in range(len(finite_vertices)):
        v1 = finite_vertices[i]
        v2 = finite_vertices[(i + 1) % len(finite_vertices)]

        # Extend ray from center through vertices to boundary
        direction = np.array(v1) - np.array(center)
        if np.linalg.norm(direction) > 0:
            direction = direction / np.linalg.norm(direction)

            # Find boundary intersection
            t_values = []
            if direction[0] != 0:
                t_x1 = (xmin - center[0]) / direction[0]
                t_x2 = (xmax - center[0]) / direction[0]
                t_values.extend([t_x1, t_x2])
            if direction[1] != 0:
                t_y1 = (ymin - center[1]) / direction[1]
                t_y2 = (ymax - center[1]) / direction[1]
                t_values.extend([t_y1, t_y2])

            # Find valid intersections
            for t in t_values:
                if t > 0:
                    intersect = center + t * direction
                    if xmin <= intersect[0] <= xmax and ymin <= intersect[1] <= ymax:
                        boundary_points.append(intersect)
                        break

    # Combine vertices and boundary points
    all_vertices = vertices + boundary_points

    if len(all_vertices) >= 3:
        # Sort vertices by angle from center
        center_np = np.array(center)
        angles = []
        for v in all_vertices:
            diff = np.array(v) - center_np
            angle = np.arctan2(diff[1], diff[0])
            angles.append(angle)

        sorted_indices = np.argsort(angles)
        sorted_vertices = [all_vertices[i] for i in sorted_indices]
        return sorted_vertices

    return vertices


def clip_polygon_to_box(vertices, boundary_box):
    """
    Clip polygon to bounding box using Sutherland-Hodgman algorithm
    """
    xmin, xmax, ymin, ymax = boundary_box

    def is_inside(point, edge):
        x, y = point
        if edge == "left":
            return x >= xmin
        elif edge == "right":
            return x <= xmax
        elif edge == "bottom":
            return y >= ymin
        elif edge == "top":
            return y <= ymax

    def compute_intersection(p1, p2, edge):
        x1, y1 = p1
        x2, y2 = p2

        if edge == "left":
            if x2 != x1:
                y = y1 + (y2 - y1) * (xmin - x1) / (x2 - x1)
                return [xmin, y]
        elif edge == "right":
            if x2 != x1:
                y = y1 + (y2 - y1) * (xmax - x1) / (x2 - x1)
                return [xmax, y]
        elif edge == "bottom":
            if y2 != y1:
                x = x1 + (x2 - x1) * (ymin - y1) / (y2 - y1)
                return [x, ymin]
        elif edge == "top":
            if y2 != y1:
                x = x1 + (x2 - x1) * (ymax - y1) / (y2 - y1)
                return [x, ymax]
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
                if is_inside(e, edge):
                    if not is_inside(s, edge):
                        intersection = compute_intersection(s, e, edge)
                        if intersection:
                            clipped.append(intersection)
                    clipped.append(e)
                elif is_inside(s, edge):
                    intersection = compute_intersection(s, e, edge)
                    if intersection:
                        clipped.append(intersection)
                s = e

    return clipped


def clip_voronoi_to_convex_hull(vor, points, hull):
    """
    Clip Voronoi diagram to convex hull of points

    Args:
        vor: Voronoi diagram object
        points: Original point cloud
        hull: ConvexHull object

    Returns:
        List of clipped regions
    """
    hull_vertices = points[hull.vertices]

    clipped_regions = []

    for point_idx, region_idx in enumerate(vor.point_region):
        region = vor.regions[region_idx]

        if not region or -1 in region:
            # Handle infinite regions by clipping to convex hull
            vertices = clip_infinite_region_to_hull(vor, point_idx, hull_vertices)
        else:
            # Finite region - just get vertices
            vertices = [vor.vertices[i] for i in region]

        if len(vertices) >= 3:  # Valid polygon
            # Clip polygon to convex hull
            clipped_vertices = clip_polygon_to_convex_hull(vertices, hull_vertices)
            if len(clipped_vertices) >= 3:
                clipped_regions.append(clipped_vertices)

    return clipped_regions


def clip_infinite_region_to_hull(vor, point_idx, hull_vertices):
    """
    Clip an infinite Voronoi region to convex hull
    """
    center = vor.points[point_idx]
    region_idx = vor.point_region[point_idx]
    region = vor.regions[region_idx]

    if not region:
        return []

    # Get finite vertices
    finite_vertices = []
    for vertex_idx in region:
        if vertex_idx >= 0:
            finite_vertices.append(vor.vertices[vertex_idx])

    if len(finite_vertices) < 2:
        # Point is likely outside or on the boundary
        return []

    # Extend infinite rays to hull boundary
    vertices = finite_vertices[:]

    # Find intersections with hull edges
    for i in range(len(finite_vertices)):
        v1 = finite_vertices[i]
        direction = np.array(v1) - np.array(center)
        if np.linalg.norm(direction) > 0:
            direction = direction / np.linalg.norm(direction)

            # Find intersection with hull edges
            intersection = find_ray_hull_intersection(center, direction, hull_vertices)
            if intersection is not None:
                vertices.append(intersection)

    if len(vertices) >= 3:
        # Sort vertices by angle from center
        center_np = np.array(center)
        angles = []
        for v in vertices:
            diff = np.array(v) - center_np
            angle = np.arctan2(diff[1], diff[0])
            angles.append(angle)

        sorted_indices = np.argsort(angles)
        sorted_vertices = [vertices[i] for i in sorted_indices]
        return sorted_vertices

    return vertices


def find_ray_hull_intersection(start, direction, hull_vertices):
    """
    Find intersection of a ray with convex hull edges
    """
    min_t = float("inf")
    best_intersection = None

    for i in range(len(hull_vertices)):
        p1 = hull_vertices[i]
        p2 = hull_vertices[(i + 1) % len(hull_vertices)]

        # Line equation: p1 + s * (p2 - p1)
        # Ray equation: start + t * direction
        # Solve: p1 + s * (p2 - p1) = start + t * direction

        edge_vec = p2 - p1
        start_to_p1 = p1 - start

        # 2D cross product
        denom = edge_vec[0] * direction[1] - edge_vec[1] * direction[0]

        if abs(denom) > 1e-10:  # Lines are not parallel
            s = (direction[0] * start_to_p1[1] - direction[1] * start_to_p1[0]) / denom
            t = (edge_vec[0] * start_to_p1[1] - edge_vec[1] * start_to_p1[0]) / denom

            if 0 <= s <= 1 and t > 0 and t < min_t:  # Valid intersection
                min_t = t
                best_intersection = start + t * direction

    return best_intersection


def clip_polygon_to_convex_hull(vertices, hull_vertices):
    """
    Clip polygon to convex hull using Sutherland-Hodgman algorithm
    """
    clipped = vertices[:]

    # For each edge of the convex hull
    for i in range(len(hull_vertices)):
        if not clipped:
            break

        p1 = hull_vertices[i]
        p2 = hull_vertices[(i + 1) % len(hull_vertices)]

        input_list = clipped[:]
        clipped = []

        if input_list:
            s = input_list[-1]

            for e in input_list:
                if is_inside_edge(e, p1, p2):
                    if not is_inside_edge(s, p1, p2):
                        intersection = line_intersection(s, e, p1, p2)
                        if intersection is not None:
                            clipped.append(intersection)
                    clipped.append(e)
                elif is_inside_edge(s, p1, p2):
                    intersection = line_intersection(s, e, p1, p2)
                    if intersection is not None:
                        clipped.append(intersection)
                s = e

    return clipped


def is_inside_edge(point, edge_p1, edge_p2):
    """
    Check if point is on the inside (left) side of the directed edge
    """
    return (
        (edge_p2[0] - edge_p1[0]) * (point[1] - edge_p1[1])
        - (edge_p2[1] - edge_p1[1]) * (point[0] - edge_p1[0])
    ) >= 0


def line_intersection(p1, p2, p3, p4):
    """
    Find intersection of two line segments
    """
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3
    x4, y4 = p4

    denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)

    if abs(denom) < 1e-10:
        return None

    t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom

    x = x1 + t * (x2 - x1)
    y = y1 + t * (y2 - y1)

    return [x, y]


def clip_voronoi_with_proper_closure(
    vor, points, boundary_margin, use_convex_hull=True
):
    """
    Advanced Voronoi clipping with proper closure of infinite cells
    """
    # Calculate boundary
    min_x, min_y = np.min(points, axis=0)
    max_x, max_y = np.max(points, axis=0)

    if use_convex_hull and len(points) >= 3:
        try:
            hull = ConvexHull(points)
            boundary_vertices = points[hull.vertices]
            boundary_type = "convex_hull"
        except:
            boundary_vertices = None
            boundary_type = "box"
    else:
        boundary_vertices = None
        boundary_type = "box"

    # Expand boundary for infinite cells
    if boundary_type == "box":
        boundary_box = [
            min_x - boundary_margin,
            max_x + boundary_margin,
            min_y - boundary_margin,
            max_y + boundary_margin,
        ]

    clipped_regions = []

    for point_idx, region_idx in enumerate(vor.point_region):
        region = vor.regions[region_idx]
        center = vor.points[point_idx]

        if not region:
            continue

        if -1 in region:
            # Handle infinite regions with proper closure
            vertices = handle_infinite_region_advanced(
                vor,
                point_idx,
                boundary_vertices,
                boundary_box if boundary_type == "box" else None,
            )
        else:
            # Finite region
            vertices = [vor.vertices[i] for i in region if i >= 0]

        if len(vertices) >= 3:
            # Ensure proper ordering and closure
            vertices = order_vertices_ccw(vertices, center)
            if len(vertices) >= 3:
                clipped_regions.append(vertices)

    return clipped_regions


def handle_infinite_region_advanced(vor, point_idx, boundary_vertices, boundary_box):
    """
    Advanced handling of infinite Voronoi regions with proper edge closure
    """
    center = vor.points[point_idx]
    region_idx = vor.point_region[point_idx]
    region = vor.regions[region_idx]

    # Get finite vertices
    finite_vertices = []
    for vertex_idx in region:
        if vertex_idx >= 0:
            finite_vertices.append(vor.vertices[vertex_idx])

    if len(finite_vertices) == 0:
        return []

    # Find ridge points (neighboring Voronoi sites)
    ridge_points = []
    for ridge in vor.ridge_points:
        if point_idx in ridge:
            other_idx = ridge[1] if ridge[0] == point_idx else ridge[0]
            ridge_points.append(vor.points[other_idx])

    # Generate rays from finite vertices
    rays = []
    for vertex in finite_vertices:
        # Direction perpendicular to the line connecting center to vertex
        to_vertex = np.array(vertex) - np.array(center)
        # Find the ridge point that corresponds to this vertex
        best_ridge = None
        min_dist = float("inf")

        for ridge_point in ridge_points:
            # Check if this ridge corresponds to this vertex direction
            to_ridge = np.array(ridge_point) - np.array(center)
            # Perpendicular to ridge direction should align with vertex direction
            perp = np.array([-to_ridge[1], to_ridge[0]])
            perp = perp / np.linalg.norm(perp)

            dot = np.dot(to_vertex / np.linalg.norm(to_vertex), perp)
            if abs(dot) < 0.5:  # Close to perpendicular
                dist = np.linalg.norm(
                    np.array(vertex)
                    - (np.array(center) + np.dot(to_vertex, perp) * perp)
                )
                if dist < min_dist:
                    min_dist = dist
                    best_ridge = ridge_point

        if best_ridge is not None:
            # Ray direction is perpendicular to line from center to ridge point
            to_ridge = np.array(best_ridge) - np.array(center)
            ray_dir = np.array([-to_ridge[1], to_ridge[0]])
            ray_dir = ray_dir / np.linalg.norm(ray_dir)

            # Ensure ray points outward from center
            if np.dot(to_vertex, ray_dir) < 0:
                ray_dir = -ray_dir

            rays.append((vertex, ray_dir))

    # Find intersections with boundary
    boundary_intersections = []

    if boundary_box:
        xmin, xmax, ymin, ymax = boundary_box
        boundary_edges = [
            ([xmin, ymin], [xmax, ymin]),  # bottom
            ([xmax, ymin], [xmax, ymax]),  # right
            ([xmax, ymax], [xmin, ymax]),  # top
            ([xmin, ymax], [xmin, ymin]),  # left
        ]
    else:
        # Use convex hull edges
        boundary_edges = []
        for i in range(len(boundary_vertices)):
            p1 = boundary_vertices[i]
            p2 = boundary_vertices[(i + 1) % len(boundary_vertices)]
            boundary_edges.append((p1, p2))

    # Find ray-boundary intersections
    for vertex, ray_dir in rays:
        best_intersection = None
        min_t = float("inf")

        for edge_start, edge_end in boundary_edges:
            intersection = ray_line_intersection(vertex, ray_dir, edge_start, edge_end)
            if intersection is not None:
                t = np.linalg.norm(np.array(intersection) - np.array(vertex))
                if t > 1e-6 and t < min_t:  # Valid intersection, not too close
                    min_t = t
                    best_intersection = intersection

        if best_intersection is not None:
            boundary_intersections.append(best_intersection)

    # Combine finite vertices and boundary intersections
    all_vertices = finite_vertices + boundary_intersections

    # Add boundary vertices that are inside this Voronoi cell
    if boundary_vertices is not None:
        for bv in boundary_vertices:
            if is_point_in_voronoi_cell(bv, center, ridge_points):
                all_vertices.append(bv)

    return all_vertices


def ray_line_intersection(ray_start, ray_dir, line_start, line_end):
    """
    Find intersection between ray and line segment
    """
    ray_start = np.array(ray_start)
    ray_dir = np.array(ray_dir)
    line_start = np.array(line_start)
    line_end = np.array(line_end)

    line_dir = line_end - line_start

    # Solve: ray_start + t * ray_dir = line_start + s * line_dir
    # t * ray_dir - s * line_dir = line_start - ray_start

    det = ray_dir[0] * (-line_dir[1]) - ray_dir[1] * (-line_dir[0])

    if abs(det) < 1e-10:
        return None  # Parallel

    diff = line_start - ray_start
    t = (diff[0] * (-line_dir[1]) - diff[1] * (-line_dir[0])) / det
    s = (ray_dir[0] * diff[1] - ray_dir[1] * diff[0]) / det

    if t > 1e-6 and 0 <= s <= 1:  # Ray forward, line segment
        return ray_start + t * ray_dir

    return None


def is_point_in_voronoi_cell(point, center, ridge_points):
    """
    Check if a point is inside a Voronoi cell
    """
    point = np.array(point)
    center = np.array(center)

    # Point is in Voronoi cell if it's closer to center than to any ridge point
    dist_to_center = np.linalg.norm(point - center)

    for ridge_point in ridge_points:
        ridge_point = np.array(ridge_point)
        dist_to_ridge = np.linalg.norm(point - ridge_point)
        if dist_to_ridge <= dist_to_center:
            return False

    return True


def order_vertices_ccw(vertices, center):
    """
    Order vertices counter-clockwise around center
    """
    if len(vertices) < 3:
        return vertices

    center = np.array(center)

    # Calculate angles from center to each vertex
    angle_vertex_pairs = []
    for vertex in vertices:
        vertex_arr = np.array(vertex)
        diff = vertex_arr - center
        angle = np.arctan2(diff[1], diff[0])
        angle_vertex_pairs.append((angle, vertex_arr))

    # Sort by angle
    angle_vertex_pairs.sort(key=lambda x: x[0])
    sorted_vertices = [vertex for angle, vertex in angle_vertex_pairs]

    return sorted_vertices


def create_shape_from_points(points, boundary_margin=20, use_convex_hull=True):
    """
    Function to create a Voronoi diagram from a point cloud and generate a TopoDS_Shape

    Args:
        points: Point cloud. Shape is (num_points, 2).
        boundary_margin: Margin around points for boundary box
        use_convex_hull: If True, use convex hull as boundary, otherwise use bounding box

    Returns:
        Generated shape.
    """
    # Compute Voronoi diagram
    vor = Voronoi(points)

    # Use improved clipping algorithm
    clipped_regions = clip_voronoi_with_proper_closure(
        vor, points, boundary_margin, use_convex_hull
    )

    builder = BRep_Builder()
    compound = TopoDS_Compound()
    builder.MakeCompound(compound)

    # Create shapes from clipped regions
    for region_vertices in clipped_regions:
        if len(region_vertices) < 3:
            continue  # Skip invalid polygons

        try:
            # Create wire from vertices
            wire_builder = BRepBuilderAPI_MakeWire()

            for j in range(len(region_vertices)):
                start_vertex = region_vertices[j]
                end_vertex = region_vertices[(j + 1) % len(region_vertices)]

                # Set z-coordinate to 0
                start_pnt = gp_Pnt(float(start_vertex[0]), float(start_vertex[1]), 0.0)
                end_pnt = gp_Pnt(float(end_vertex[0]), float(end_vertex[1]), 0.0)

                # Check for duplicate points
                if (
                    abs(start_pnt.X() - end_pnt.X()) > 1e-6
                    or abs(start_pnt.Y() - end_pnt.Y()) > 1e-6
                ):
                    edge_shape = BRepBuilderAPI_MakeEdge(start_pnt, end_pnt).Edge()
                    wire_builder.Add(edge_shape)

            if wire_builder.IsDone():
                wire = wire_builder.Wire()
                face_maker = BRepBuilderAPI_MakeFace(wire)
                if face_maker.IsDone():
                    face_shape = face_maker.Face()
                    builder.Add(compound, face_shape)

        except Exception as e:
            print(f"Error creating face for region: {e}")
            continue

    return compound


def extract_boundary_edges_from_shape(shape):
    """
    Extract boundary edges from a TopoDS_Shape (compound of faces)

    Args:
        shape: TopoDS_Shape containing faces

    Returns:
        List of boundary edges (edges that belong to only one face)
    """
    # Create a map of edges to faces
    edge_face_map = TopTools_IndexedDataMapOfShapeListOfShape()
    topexp.MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, edge_face_map)

    boundary_edges = []

    # Iterate through all edges
    for i in range(1, edge_face_map.Size() + 1):
        edge = edge_face_map.FindKey(i)
        face_list = edge_face_map.FindFromIndex(i)

        # If edge belongs to only one face, it's a boundary edge
        if face_list.Size() == 1:
            boundary_edges.append(edge)

    return boundary_edges


def display_boundary_edges(shape, display, color="green"):
    """
    Display boundary edges extracted from shape

    Args:
        shape: TopoDS_Shape to extract boundaries from
        display: OpenCASCADE display object
        color: Color for boundary edges
    """
    boundary_edges = extract_boundary_edges_from_shape(shape)

    print(f"Found {len(boundary_edges)} boundary edges")

    for edge in boundary_edges:
        display.DisplayShape(edge, color=color, update=False)

    return boundary_edges


def create_boundary_wire_from_shape(shape):
    """
    Create boundary wire(s) from shape by connecting boundary edges

    Args:
        shape: TopoDS_Shape to extract boundaries from

    Returns:
        List of connected boundary wires
    """
    from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeWire

    boundary_edges = extract_boundary_edges_from_shape(shape)

    if not boundary_edges:
        return []

    # Simple approach: try to create one wire from all boundary edges
    wires = []

    # Create individual wires from each boundary edge
    for edge in boundary_edges:
        try:
            wire_builder = BRepBuilderAPI_MakeWire(edge)
            if wire_builder.IsDone():
                wires.append(wire_builder.Wire())
        except Exception as e:
            print(f"Error creating wire from edge: {e}")
            continue

    print(f"Created {len(wires)} individual boundary wires")
    return wires


def get_voronoi_outline_direct(
    vor, points, boundary_margin=20, display_outline=True, display=None
):
    """
    Extract the true outline of Voronoi diagram directly from Voronoi data

    Args:
        vor: Voronoi diagram object
        points: Input points
        boundary_margin: Margin for boundary calculation
        display_outline: Whether to display the outline
        display: OpenCASCADE display object (if displaying)

    Returns:
        List of outline edges (as coordinate pairs)
    """
    # Calculate boundary
    min_x, min_y = np.min(points, axis=0)
    max_x, max_y = np.max(points, axis=0)
    boundary_box = [
        min_x - boundary_margin,
        max_x + boundary_margin,
        min_y - boundary_margin,
        max_y + boundary_margin,
    ]

    # Find boundary edges by checking which Voronoi edges intersect the boundary
    outline_edges = []

    # Get clipped regions to find boundary
    clipped_regions = clip_voronoi_with_proper_closure(
        vor, points, boundary_margin, use_convex_hull=True
    )

    # Collect all edges from clipped regions
    all_edges = []
    edge_count = {}

    for region in clipped_regions:
        for i in range(len(region)):
            p1 = region[i]
            p2 = region[(i + 1) % len(region)]

            # Skip invalid points
            p1_arr = np.array(p1)
            p2_arr = np.array(p2)
            if not (np.all(np.isfinite(p1_arr)) and np.all(np.isfinite(p2_arr))):
                continue

            # Skip degenerate edges
            dist = np.linalg.norm(p2_arr - p1_arr)
            if dist < 1e-6:
                continue

            # Create edge key (order independent)
            edge_key = tuple(sorted([tuple(p1), tuple(p2)]))

            if edge_key in edge_count:
                edge_count[edge_key] += 1
            else:
                edge_count[edge_key] = 1
                all_edges.append((p1, p2))

    # Boundary edges appear only once (not shared between regions)
    boundary_edges = []
    for p1, p2 in all_edges:
        edge_key = tuple(sorted([tuple(p1), tuple(p2)]))
        if edge_count[edge_key] == 1:
            # Check if this is a valid edge (points are different)
            p1_arr = np.array(p1)
            p2_arr = np.array(p2)
            dist = np.linalg.norm(p2_arr - p1_arr)

            # Only add edges with sufficient length and finite coordinates
            if (
                dist > 1e-6
                and np.all(np.isfinite(p1_arr))
                and np.all(np.isfinite(p2_arr))
            ):
                boundary_edges.append((p1, p2))
            else:
                print(
                    f"Skipping invalid boundary edge: dist={dist:.2e}, p1={p1}, p2={p2}"
                )

    if display_outline and display is not None:
        print(f"Displaying {len(boundary_edges)} outline edges...")
        for p1, p2 in boundary_edges:
            try:
                pnt1 = gp_Pnt(float(p1[0]), float(p1[1]), 0.0)
                pnt2 = gp_Pnt(float(p2[0]), float(p2[1]), 0.0)

                # Check if points are different enough
                dist = ((pnt1.X() - pnt2.X()) ** 2 + (pnt1.Y() - pnt2.Y()) ** 2) ** 0.5
                if dist > 1e-6:  # Only create edge if points are sufficiently different
                    edge_maker = BRepBuilderAPI_MakeEdge(pnt1, pnt2)
                    if edge_maker.IsDone():
                        edge = edge_maker.Edge()
                        display.DisplayShape(edge, color="red", update=False)
                    else:
                        print(
                            f"Failed to create edge from ({p1[0]:.3f}, {p1[1]:.3f}) to ({p2[0]:.3f}, {p2[1]:.3f})"
                        )
                else:
                    print(f"Skipping degenerate edge: points too close ({dist:.2e})")
            except Exception as e:
                print(
                    f"Error creating edge from ({p1[0]:.3f}, {p1[1]:.3f}) to ({p2[0]:.3f}, {p2[1]:.3f}): {e}"
                )
                continue

    return boundary_edges


def create_outline_wire(boundary_edges):
    """
    Create a connected wire from boundary edges

    Args:
        boundary_edges: List of (p1, p2) coordinate pairs

    Returns:
        Connected outline as OpenCASCADE wire or None if connection fails
    """
    if not boundary_edges:
        return None

    try:
        # Convert to OpenCASCADE edges
        occ_edges = []
        for p1, p2 in boundary_edges:
            pnt1 = gp_Pnt(float(p1[0]), float(p1[1]), 0.0)
            pnt2 = gp_Pnt(float(p2[0]), float(p2[1]), 0.0)

            # Check if points are different enough
            dist = ((pnt1.X() - pnt2.X()) ** 2 + (pnt1.Y() - pnt2.Y()) ** 2) ** 0.5
            if dist > 1e-6:  # Only create edge if points are sufficiently different
                edge_maker = BRepBuilderAPI_MakeEdge(pnt1, pnt2)
                if edge_maker.IsDone():
                    edge = edge_maker.Edge()
                    occ_edges.append(edge)
                else:
                    print(
                        f"Failed to create edge from ({p1[0]:.3f}, {p1[1]:.3f}) to ({p2[0]:.3f}, {p2[1]:.3f})"
                    )
            else:
                print(
                    f"Skipping degenerate edge in wire creation: points too close ({dist:.2e})"
                )

        # Try to connect edges into a wire
        from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeWire

        wire_builder = BRepBuilderAPI_MakeWire()

        for edge in occ_edges:
            wire_builder.Add(edge)

        if wire_builder.IsDone():
            return wire_builder.Wire()
        else:
            print("Failed to create connected outline wire")
            return None

    except Exception as e:
        print(f"Error creating outline wire: {e}")
        return None


def get_voronoi_outline(points, display_outline=True, display=None):
    """
    Extract the outline (boundary) of Voronoi diagram

    Args:
        points: Input points for Voronoi diagram
        display_outline: Whether to display the outline
        display: OpenCASCADE display object (if displaying)

    Returns:
        List of boundary edges forming the outline
    """
    # Create 3D shape from points
    shape = create_shape_from_points(points)

    # Extract boundary edges
    boundary_edges = extract_boundary_edges_from_shape(shape)

    if display_outline and display is not None:
        print(f"Displaying {len(boundary_edges)} outline edges...")
        for edge in boundary_edges:
            display.DisplayShape(edge, color="red", update=False)

    return boundary_edges


if __name__ == "__main__":
    # Initialize display
    display, start_display, add_menu, add_function_to_menu = init_display()

    # Generate points (restricted to a plane)
    num_points = 30
    points = generate_points_on_plane(num_points)

    print(f"Generated {len(points)} points")
    print(
        f"Points range: X=[{np.min(points[:, 0]):.2f}, {np.max(points[:, 0]):.2f}], "
        f"Y=[{np.min(points[:, 1]):.2f}, {np.max(points[:, 1]):.2f}]"
    )

    # Create Voronoi diagram for wireframe display
    print("Computing Voronoi diagram...")
    vor = Voronoi(points)

    # Calculate boundary box
    min_x, min_y = np.min(points, axis=0)
    max_x, max_y = np.max(points, axis=0)
    boundary_margin = 20
    boundary_box = [
        min_x - boundary_margin,
        max_x + boundary_margin,
        min_y - boundary_margin,
        max_y + boundary_margin,
    ]

    # Create 3D shape from points
    print("Creating 3D Voronoi shape...")
    shape = create_shape_from_points(points)

    # Display all Voronoi edges (blue lines) - both internal and boundary
    print("Displaying all Voronoi edges...")
    # Get boundary box for filtering
    min_x, min_y = np.min(points, axis=0)
    max_x, max_y = np.max(points, axis=0)
    margin = boundary_margin

    for ridge in vor.ridge_vertices:
        if -1 not in ridge:  # Only finite edges
            v1 = vor.vertices[ridge[0]]
            v2 = vor.vertices[ridge[1]]

            # Filter edges that are too far outside the point cloud
            if (
                min_x - margin <= v1[0] <= max_x + margin
                and min_y - margin <= v1[1] <= max_y + margin
                and min_x - margin <= v2[0] <= max_x + margin
                and min_y - margin <= v2[1] <= max_y + margin
            ):

                pnt1 = gp_Pnt(float(v1[0]), float(v1[1]), 0.0)
                pnt2 = gp_Pnt(float(v2[0]), float(v2[1]), 0.0)
                edge = BRepBuilderAPI_MakeEdge(pnt1, pnt2).Edge()
                display.DisplayShape(edge, color="blue1", update=False)

    # Draw points as vertices
    print("Displaying points...")
    for i, point in enumerate(points):
        pnt = gp_Pnt(float(point[0]), float(point[1]), 0.0)
        vertex_shape = BRepBuilderAPI_MakeVertex(pnt).Shape()
        display.DisplayShape(vertex_shape, color="red", update=False)

    # Fit view to show all geometry
    display.FitAll()

    print("Starting display...")
    start_display()

# https://github.com/tpaviot/pythonocc-core/issues/1410
# Handling Infinite Voronoi Cells in Point Cloud

# I tried to generate a Voronoi diagram from a point cloud and convert it into a TopoDS_Shape.
# The code works well for finite cells, but I'm struggling with handling infinite Voronoi cells.
# Specifically, I want to close these infinite cells to make them finite.
#
# Here's a snippet of my code:

# How can I modify this code to handle infinite Voronoi cells and close them to make them finite?
# By the way, most of this code was generated using Copilot.
