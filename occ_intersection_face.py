from OCC.Core.gp import gp_Pnt, gp_Pnt2d, gp_Ax2, gp_Dir
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut
from OCC.Core.BRepFilletAPI import BRepFilletAPI_MakeFillet
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_EDGE, TopAbs_FACE
from OCC.Core.TopoDS import topods
from OCC.Core.TColgp import TColgp_Array1OfPnt2d
from OCC.Core.BRep import BRep_Tool
from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRepBndLib import brepbndlib
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeShell, BRepBuilderAPI_Sewing
from OCC.Core.BRepExtrema import BRepExtrema_DistShapeShape
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Common
from OCC.Core.TopoDS import TopoDS_Shape
from OCC.Core.TopAbs import TopAbs_COMPOUND
from OCC.Core.BRepGProp import brepgprop
from OCC.Core.GProp import GProp_GProps
from OCC.Core.ShapeAnalysis import ShapeAnalysis_Surface
from OCC.Display.SimpleGui import init_display

display, start_display, add_menu, add_function_to_menu = init_display()

# パラメータ
D = 20.0
L = 40.0
L1 = 10.0
L2 = 10.0
H = 10.0
D2 = 5.0
FSH = 0  # フィレット形状（例: 0=CONVEX, 1=CONCAVE）

# Part 1: 複雑な直方体（例としてBox）
box1 = BRepPrimAPI_MakeBox(
    gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0.1, 1), gp_Dir(1, 0, 0)), D, L, H
).Shape()

# Part 2: 下部のヘルパーBox
# box2の位置をrectangle_L2_H_shapeD32の底面に接するように調整
box2 = BRepPrimAPI_MakeBox(
    gp_Pnt(-D2, -H, -H),  # Z座標を-Hに変更して接触させる
    4 * D,
    L + L1 + 2 * H,
    H * 1.1,  # 高さをHに変更
).Shape()

# 融合
fused_shape = BRepAlgoAPI_Fuse(box1, box2).Shape()

# 共通部分を表示（デバッグ用）
common_shape = BRepAlgoAPI_Common(box1, box2).Shape()


# box1/box2の全faceをリスト化
def get_faces(shape):
    faces = []
    exp = TopExp_Explorer(shape, TopAbs_FACE)
    while exp.More():
        faces.append(topods.Face(exp.Current()))
        exp.Next()
    return faces


faces_box1 = get_faces(box1)
faces_box2 = get_faces(box2)


def is_same_plane(face1, face2, tol=1e-6):
    # 2つの面が同一平面上かどうか判定
    surf1 = BRep_Tool.Surface(face1)
    surf2 = BRep_Tool.Surface(face2)
    # 平面以外は除外
    if not surf1.IsKind("Geom_Plane") or not surf2.IsKind("Geom_Plane"):
        return False
    pln1 = surf1.Plane()
    pln2 = surf2.GetObject().Plane()
    # 法線ベクトルと位置が一致するか
    n1 = pln1.Axis().Direction()
    n2 = pln2.Axis().Direction()
    p1 = pln1.Location()
    p2 = pln2.Location()
    same_normal = n1.IsParallel(n2, tol)
    same_pos = abs(pln1.Distance(p2)) < tol
    return same_normal and same_pos


def face_center(face):
    props = GProp_GProps()
    brepgprop.SurfaceProperties(face, props)
    return props.CentreOfMass()


def is_on_surface(face_target, face_ref, tol=1e-5):
    # face_targetの重心がface_ref上にあるか判定
    center = face_center(face_target)
    surf = BRep_Tool.Surface(face_ref)
    sas = ShapeAnalysis_Surface(surf)
    uv2d = sas.ValueOfUV(center, 1e-7)
    p_proj = sas.Value(uv2d.X(), uv2d.Y())
    dist = center.Distance(p_proj)
    # UVが有効範囲内かつ距離が十分小さい
    # umin, umax, vmin, vmax = BRep_Tool().UVBounds(face_ref)
    # on_uv = (umin - tol) <= u <= (umax + tol) and (vmin - tol) <= v <= (vmax + tol)
    # return dist < tol and on_uv


# 共通部分のfaceのうち、box1/box2のface上にあるものを抽出
# intersection_faces = []
# exp = TopExp_Explorer(common_shape, TopAbs_FACE)
# while exp.More():
#     f = topods.Face(exp.Current())
#     for f1 in faces_box1 + faces_box2:
#         if is_on_surface(f, f1):
#             intersection_faces.append(f)
#             break
#     exp.Next()

# 表示
# for f in intersection_faces:
#     display.DisplayShape(f, color="BLUE", transparency=0.0, update=True)

# 差分形状を作成
cut1 = BRepAlgoAPI_Cut(box1, common_shape).Shape()
cut2 = BRepAlgoAPI_Cut(box2, common_shape).Shape()

# 各差分形状の全faceを取得
faces_cut1 = get_faces(cut1)
faces_cut2 = get_faces(cut2)

# TopoDS_Face同士が一致するものを抽出
common_faces = []
for f1 in faces_cut1:
    for f2 in faces_cut2:
        if f1.IsSame(f2):
            common_faces.append(f1)
            break

print(f"共通面の数: {len(common_faces)}")
# 表示
for f in common_faces:
    display.DisplayShape(f, color="GREEN", transparency=0.0, update=True)

# ...既存の表示...
display.DisplayShape(common_shape, color="RED", transparency=0.5, update=True)
display.DisplayShape(fused_shape, update=True, transparency=0.7)
start_display()
