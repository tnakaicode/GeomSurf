import math
from OCC.Core.TopoDS import TopoDS_Face
from OCC.Core.gp import gp_Pnt
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakePolygon, BRepBuilderAPI_MakeFace
from OCC.Display.SimpleGui import init_display

# 初期化
display, start_display, add_menu, add_function_to_menu = init_display()

# 10角形の頂点を計算
center_x, center_y, radius = 0, 0, 50
points = []

for i in range(10):
    angle = 2 * math.pi * i / 10
    # 凹んだ部分の頂点を計算
    if i == 3:  # 一部を凹にする、例えば第4頂点を凹に
        x = center_x + (radius / 2) * math.cos(angle)
        y = center_y + (radius / 2) * math.sin(angle)
    else:
        x = center_x + radius * math.cos(angle)
        y = center_y + radius * math.sin(angle)
    points.append(gp_Pnt(x, y, 0))

# ポリゴンの作成
polygon = BRepBuilderAPI_MakePolygon()
for p in points:
    polygon.Add(p)
polygon.Close()

# 面の作成
face = BRepBuilderAPI_MakeFace(polygon.Wire())

# 表示
display.DisplayShape(face.Shape(), update=True)
start_display()
