from OCC.Core.gp import gp_Pnt2d, gp_Vec2d, gp_Dir2d
from OCC.Core.GCE2d import GCE2d_MakeSegment, GCE2d_MakeArcOfCircle
from OCC.Core.Geom2d import Geom2d_Line
from OCC.Core.Geom2dAPI import Geom2dAPI_PointsToBSpline, Geom2dAPI_InterCurveCurve
from OCC.Core.TColgp import TColgp_Array1OfPnt2d
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge2d
from OCC.Core.Geom2d import Geom2d_TrimmedCurve
from OCC.Display.SimpleGui import init_display
from OCC.Extend.ShapeFactory import make_wire
import math


def create_fillet(spline, segment, fillet_radius):

    # 曲線上の点P0を指定
    curve_param = spline.LastParameter()
    curve_point = spline.Value(curve_param)

    # 接線の方向を計算
    tangent_vec = gp_Vec2d(spline.Value(curve_param - 0.01), curve_point)
    tangent_dir = gp_Dir2d(tangent_vec)
    tangent_line_geom = Geom2d_Line(curve_point, tangent_dir)

    # 接線を延長して直線と交差する点を見つける
    intersector = Geom2dAPI_InterCurveCurve(tangent_line_geom, segment, 0.001)
    if intersector.NbPoints() == 0:
        raise RuntimeError("交差点が見つかりませんでした")

    p1 = intersector.Point(1)

    # P0-P1の線分を作成
    seg_p0_p1 = GCE2d_MakeSegment(curve_point, p1)

    # P0-P1の線分と直線のフィレットを作成
    distance = (
        (curve_point.X() - p1.X()) ** 2 + (curve_point.Y() - p1.Y()) ** 2
    ) ** 0.5
    if distance >= fillet_radius:
        # 円の中心点を計算
        h = math.sqrt(fillet_radius**2 - (distance / 2) ** 2)
        if h >= 0:
            mid_point_x = (curve_point.X() + p1.X()) / 2
            mid_point_y = (curve_point.Y() + p1.Y()) / 2
            dx = (p1.Y() - curve_point.Y()) / distance
            dy = (curve_point.X() - p1.X()) / distance
            center_point1 = gp_Pnt2d(mid_point_x + h * dx, mid_point_y - h * dy)
            center_point2 = gp_Pnt2d(mid_point_x - h * dx, mid_point_y + h * dy)

            try:
                fillet_arc1 = GCE2d_MakeArcOfCircle(
                    curve_point, center_point1, p1
                ).Value()
                fillet_arc2 = GCE2d_MakeArcOfCircle(
                    curve_point, center_point2, p1
                ).Value()
                fillet_arc = (
                    fillet_arc1
                    if (fillet_arc1.FirstParameter() == curve_point)
                    else fillet_arc2
                )
            except:
                raise RuntimeError("フィレットの作成に失敗しました")
        else:
            raise RuntimeError("フィレット半径が大きすぎます")

    # スプライン、直線、フィレットをエッジとして作成
    edge1 = BRepBuilderAPI_MakeEdge2d(spline).Edge()
    edge2 = BRepBuilderAPI_MakeEdge2d(segment).Edge()
    # fillet_edge = BRepBuilderAPI_MakeEdge2d(Geom2d_TrimmedCurve(fillet_arc)).Edge()

    # ワイヤーを作成して表示
    wire = make_wire([edge1, edge2])

    return wire, edge1, edge2, None


if __name__ == "__main__":
    # 表示を初期化
    display, start_display, add_menu, add_function_to_menu = init_display()

    # 関数を呼び出してフィレットを作成
    curve_points = [(0, 0), (1, 0.1), (2, 0)]
    line_start = (2, 0)
    line_end = (3, -1)
    fillet_radius = 0.5

    # スプライン曲線を作成
    points = TColgp_Array1OfPnt2d(1, len(curve_points))
    for i, pt in enumerate(curve_points):
        points.SetValue(i + 1, gp_Pnt2d(pt[0], pt[1]))
    spline = Geom2dAPI_PointsToBSpline(points).Curve()

    # 直線を作成
    line_start = gp_Pnt2d(line_start[0], line_start[1])
    line_end = gp_Pnt2d(line_end[0], line_end[1])
    segment = GCE2d_MakeSegment(line_start, line_end).Value()

    wire, edge1, edge2, fillet_edge = create_fillet(spline, segment, fillet_radius)

    # ワイヤーを表示
    display.DisplayShape(wire, color="GREEN")

    # 元の線分とフィレットエッジを色分けして表示
    display.DisplayShape(edge1, color="BLUE1")
    display.DisplayShape(edge2, color="BLUE1")
    # display.DisplayShape(fillet_edge, color="YELLOW")

    # 表示を調整
    display.FitAll()
    start_display()
