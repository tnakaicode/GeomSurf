from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakePolygon, BRepBuilderAPI_MakeEdge
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.Core.gp import gp_Pnt2d, gp_Pnt
from OCCUtils.Topology import topods
from OCC.Core.BRepTopAdaptor import BRepTopAdaptor_FClass2d
from OCC.Core.TopAbs import TopAbs_IN
from OCC.Display.SimpleGui import init_display
import numpy as np

display, start_display, add_menu, add_function_to_menu = init_display()


def is_point_inside_polygon(polygon_points, test_point):
    # Create a 2D polygon
    polygon = BRepBuilderAPI_MakePolygon()
    for point in polygon_points:
        polygon.Add(gp_Pnt(point[0], point[1], 0))
    polygon.Close()
    polygon_shape = polygon.Shape()

    # Create a face from the polygon
    face = BRepBuilderAPI_MakeFace(polygon_shape, True)
    face_shape = face.Shape()

    # Convert test point to gp_Pnt2d
    test_pnt = gp_Pnt2d(test_point[0], test_point[1])

    # Check if the point is inside the polygon
    classifier = BRepTopAdaptor_FClass2d(topods.Face(face_shape), 1e-6)
    state = classifier.Perform(test_pnt)

    return state == TopAbs_IN


def draw_polygon_and_point(polygon_points, test_point):

    # Create and display the polygon
    polygon = BRepBuilderAPI_MakePolygon()
    for point in polygon_points:
        polygon.Add(gp_Pnt(point[0], point[1], 0))
    polygon.Close()

    polygon_shape = polygon.Shape()
    display.DisplayShape(polygon_shape, update=True)

    # Check if the point is inside the polygon
    is_inside = is_point_inside_polygon(polygon_points, test_point)

    # Create and display the point with respective color
    # point_shape = BRepBuilderAPI_MakeEdge(gp_Pnt(test_point[0], test_point[1], 0)).Edge()
    point_shape = gp_Pnt(test_point[0], test_point[1], 0)
    color = "red" if is_inside else "blue1"
    display.DisplayShape(point_shape, color=color, update=True)


# Example usage with concave decagon points
polygon_points = [
    [0, 0],
    [2, 1],
    [4, 0],
    [5, 2],
    [3, 4],
    [3, 6],
    [5, 8],
    [3, 10],
    [2, 9],
    [0, 10],
]
test_points = [
    [2, 5],
    [4, 5],
    [2, 8],
    [2, 1],
    [4, 3],
    [1, 1],
    [3, 9],
    [3, 2],
    [5, 9],
    [2, 2],
]
for test_point in test_points:
    result = is_point_inside_polygon(polygon_points, test_point)
    print("Is the point inside the polygon?", result)

    draw_polygon_and_point(polygon_points, test_point)

start_display()
