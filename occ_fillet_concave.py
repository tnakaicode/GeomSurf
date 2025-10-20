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
box1 = BRepPrimAPI_MakeBox(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(
    0, 0.1, 1), gp_Dir(1, 0, 0)), D, L, H).Shape()

# Part 2: 下部のヘルパーBox
# box2の位置をrectangle_L2_H_shapeD32の底面に接するように調整
box2 = BRepPrimAPI_MakeBox(
    gp_Pnt(-D2, -H, -H),  # Z座標を-Hに変更して接触させる
    4 * D, L + L1 + 2 * H, H * 1.1  # 高さをHに変更
).Shape()

# 融合
fused_shape = BRepAlgoAPI_Fuse(box1, box2).Shape()

# フィレット設定
filletD32 = BRepFilletAPI_MakeFillet(fused_shape)
filletD32.SetFilletShape(FSH)

# フィレット半径をパラメトリックに設定
parAndRad2 = TColgp_Array1OfPnt2d(1, 3)
parAndRad2.SetValue(1, gp_Pnt2d(0, 0.001))
parAndRad2.SetValue(2, gp_Pnt2d(0.5, D / 2 * (1 / 2 * (L2 / (L2 + L))) * 2))
parAndRad2.SetValue(3, gp_Pnt2d(1, D / 2 * (L2 / (L2 + L)) ** 2))

# フィレットをかけるエッジを探索
edgesD32 = []
edge_explorerD32 = TopExp_Explorer(fused_shape, TopAbs_EDGE)
while edge_explorerD32.More():
    edgeD32 = topods.Edge(edge_explorerD32.Current())
    edgesD32.append(edgeD32)
    edge_explorerD32.Next()

edges_to_filletD32 = [3, 2]  # 適宜調整

for edge_indexD32 in edges_to_filletD32:
    if edge_indexD32 < len(edgesD32):
        filletD32.Add(parAndRad2, edgesD32[edge_indexD32])

# resultD32 = filletD32.Shape()

# Part 2をCutで除去（このままだとフィレットも消える場合あり）
# cut1 = BRepAlgoAPI_Cut(resultD32, box2).Shape()


# --- フィレット面だけを抽出する例 ---
def is_fillet_face(face):
    # 例: 面積が小さい、またはZ座標が下部にある面をフィレット面とみなす
    from OCC.Core.BRepGProp import brepgprop_SurfaceProperties
    from OCC.Core.GProp import GProp_GProps
    from OCC.Core.BRep import BRep_Tool
    import numpy as np

    props = GProp_GProps()
    brepgprop_SurfaceProperties(face, props)
    area = props.Mass()
    umin, umax, vmin, vmax = BRep_Tool().UVBounds(face)
    # 面の重心
    center = props.CentreOfMass()
    # Z座標が下部、かつ面積が小さいものをフィレット面と仮定
    return area < 2.0 and center.Z() < 1.0


# フィレット面の抽出
# fillet_faces = []
# face_exp = TopExp_Explorer(resultD32, TopAbs_FACE)
# while face_exp.More():
#    face = topods_Face(face_exp.Current())
#    if is_fillet_face(face):
#        fillet_faces.append(face)
#    face_exp.Next()

# Sewingでフィレット面をまとめる
# sewing = BRepBuilderAPI_Sewing()
# for f in fillet_faces:
#    sewing.Add(f)
# sewing.Perform()
# fillet_shell = sewing.SewedShape()

def find_contact_edges(part1, part2):
    # 接触面のエッジを特定する関数
    contact_edges = []
    face_explorer1 = TopExp_Explorer(part1, TopAbs_FACE)
    while face_explorer1.More():
        face1 = topods.Face(face_explorer1.Current())
        face_explorer2 = TopExp_Explorer(part2, TopAbs_FACE)
        while face_explorer2.More():
            face2 = topods.Face(face_explorer2.Current())
            # 面同士が接触しているかを確認
            if are_faces_touching(face1, face2):
                # 接触している面のエッジを取得
                edge_explorer = TopExp_Explorer(face1, TopAbs_EDGE)
                while edge_explorer.More():
                    edge = topods.Edge(edge_explorer.Current())
                    if edge not in contact_edges:  # 重複を防ぐ
                        contact_edges.append(edge)
                    edge_explorer.Next()
            face_explorer2.Next()
        face_explorer1.Next()
    return contact_edges


def are_faces_touching(face1, face2):
    # 面同士が接触しているかを確認する関数（修正版）
    tol = 1e-6  # 許容誤差

    # face1のバウンディングボックスを取得
    bbox1 = Bnd_Box()
    brepbndlib.Add(face1, bbox1)

    # face2のバウンディングボックスを取得
    bbox2 = Bnd_Box()
    brepbndlib.Add(face2, bbox2)

    # バウンディングボックスが重なっているかを確認
    return not bbox1.IsOut(bbox2)


def are_faces_touching2(face1, face2):
    # 面同士が接触しているかを確認する関数（厳密版）
    tol = 1e-6  # 許容誤差

    # 面同士の最小距離を計算
    dist_calc = BRepExtrema_DistShapeShape(face1, face2)
    dist_calc.Perform()

    # 計算が成功し、最小距離が許容誤差以下であれば接触しているとみなす
    if dist_calc.IsDone() and dist_calc.Value() <= tol:
        return True
    return False


def find_contact_faces(part1, part2):
    # 接触面を特定する関数
    contact_faces = []
    face_explorer1 = TopExp_Explorer(part1, TopAbs_FACE)
    while face_explorer1.More():
        face1 = face_explorer1.Current()
        face_explorer2 = TopExp_Explorer(part2, TopAbs_FACE)
        while face_explorer2.More():
            face2 = face_explorer2.Current()
            # 面同士が接触しているかを確認
            if are_faces_touching2(face1, face2):
                if face1 not in contact_faces:  # 重複を防ぐ
                    contact_faces.append(face1)
                if face2 not in contact_faces:  # 重複を防ぐ
                    contact_faces.append(face2)
            face_explorer2.Next()
        face_explorer1.Next()
    return contact_faces


def get_edges_from_faces(faces):
    # 接触面のエッジを取得する関数
    edges = []
    for face in faces:
        edge_explorer = TopExp_Explorer(face, TopAbs_EDGE)
        while edge_explorer.More():
            edge = topods.Edge(edge_explorer.Current())
            if edge not in edges:  # 重複を防ぐ
                edges.append(edge)
            edge_explorer.Next()
    return edges

# 面同士が接触しているかを確認する関数（Boolean演算版）


def are_faces_touching_with_boolean(face1, face2):
    # Boolean演算で共通部分を計算
    common = BRepAlgoAPI_Common(face1, face2).Shape()

    # 共通部分が空でないかを確認
    if not common.IsNull():
        # 共通部分がCOMPOUNDでない場合も考慮
        if common.ShapeType() != TopAbs_COMPOUND:
            return True
    return False


def find_contact_face_pairs(part1, part2):
    """
    接触している面同士を2つ組（タプル）として返す関数
    """
    contact_face_pairs = []
    face_explorer1 = TopExp_Explorer(part1, TopAbs_FACE)
    while face_explorer1.More():
        face1 = topods.Face(face_explorer1.Current())
        face_explorer2 = TopExp_Explorer(part2, TopAbs_FACE)
        while face_explorer2.More():
            face2 = topods.Face(face_explorer2.Current())
            # 面同士が接触しているかを確認
            if are_faces_touching_with_boolean(face1, face2):  # 厳密な接触判定
                contact_face_pairs.append((face1, face2))  # ペアとして追加
            face_explorer2.Next()
        face_explorer1.Next()
    return contact_face_pairs


# 接触エッジを取得
contact_edges = find_contact_edges(box1, box2)

# 接触している面のペアを取得
contact_face_pairs = find_contact_face_pairs(box2, box1)
print("Contact face pairs count:", len(contact_face_pairs))
for i, (face1, face2) in enumerate(contact_face_pairs):
    print(f"Pair {i + 1}: Face1={face1}, Face2={face2}")

# 交差部分のエッジを取得する関数


def get_intersection_edges(shape1, shape2):
    # 共通部分を計算
    common_shape = BRepAlgoAPI_Common(shape1, shape2).Shape()

    # 共通部分のエッジを探索
    intersection_edges = []
    edge_explorer = TopExp_Explorer(common_shape, TopAbs_EDGE)
    while edge_explorer.More():
        edge = topods.Edge(edge_explorer.Current())
        intersection_edges.append(edge)
        edge_explorer.Next()
    return intersection_edges


# 交差部分のエッジを取得
intersection_edges = get_intersection_edges(box1, box2)

# デバッグ: 取得したエッジの数を確認
print(f"Number of intersection edges: {len(intersection_edges)}")
for i, edge in enumerate(intersection_edges):
    print(f"Edge {i + 1}: {edge}")

# 共通部分を表示（デバッグ用）
common_shape = BRepAlgoAPI_Common(box1, box2).Shape()
display.DisplayShape(common_shape, color='RED', transparency=0.5, update=True)

# フィレットを作成
fillet = BRepFilletAPI_MakeFillet(fused_shape)
fillet.SetFilletShape(1)
for edge in intersection_edges:
    fillet.Add(5.0, edge)  # 半径5.0のフィレットを適用

# フィレット結果を取得
# result_shape = fillet.Shape()

# --- 表示 ---
display.DisplayShape(fused_shape, update=True, transparency=0.7)
# display.DisplayShape(box1, color='GREEN', update=True, transparency=0.9)
# display.DisplayShape(box2, color='YELLOW', update=True, transparency=0.9)
display.DisplayShape(intersection_edges, color='BLUE1', update=True)
# display.DisplayShape(fillet_shell, color='RED', transparency=0.7, update=True)
start_display()
