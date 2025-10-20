import numpy as np
from scipy.spatial import Delaunay
from matplotlib import colormaps
from OCC.Core.MeshDS import MeshDS_DataSource
from OCC.Core.MeshVS import (
    MeshVS_DMF_OCCMask,
    MeshVS_Mesh,
    MeshVS_ElementalColorPrsBuilder,
    MeshVS_DA_ShowEdges,
    MeshVS_DMF_ElementalColorDataPrs,
    MeshVS_DMF_NodalColorDataPrs,
    MeshVS_NodalColorPrsBuilder,
)
from OCC.Display.SimpleGui import init_display
from OCC.Core.Aspect import Aspect_SequenceOfColor
from OCC.Core.Quantity import (
    Quantity_Color,
    Quantity_TOC_RGB,
    Quantity_NOC_PURPLE,
    Quantity_NOC_ORANGE,
    Quantity_NOC_GREEN,
    Quantity_NOC_BLUE1,
    Quantity_NOC_BLACK,
)
from OCC.Core.TColStd import TColStd_DataMapOfIntegerReal


def create_mesh_data(xx, yy, zz, data):
    # メッシュデータの生成
    xyz = np.column_stack((xx.flatten(), yy.flatten(), zz.flatten()))
    tri = Delaunay(xyz[:, :2])

    vertices = xyz
    faces = tri.simplices

    # 面の色の計算
    # facesに格納されたインデックスに対応するdataの平均を計算
    face_values = np.mean(data.flatten()[faces], axis=-1)

    # メッシュデータソースの作成
    mesh_ds = MeshDS_DataSource(vertices, faces)

    return mesh_ds, face_values, vertices, faces


if __name__ == "__main__":
    # 表示の初期化
    display, start_display, add_menu, add_function_to_menu = init_display()

    # データセット1の生成
    X1, Y1 = 50, 50
    x1 = np.linspace(-5, 0, X1)
    y1 = np.linspace(-5, 0, Y1)
    xx1, yy1 = np.meshgrid(x1, y1, sparse=False)
    zz1 = np.sin(xx1 ** 2 + yy1 ** 2) / (xx1 ** 2 + yy1 ** 2)
    data1 = np.sin(xx1 ** 2 + 2 * yy1 ** 2) / (xx1 ** 2 + yy1 ** 2)

    # データセット2の生成
    X2, Y2 = 50, 50
    x2 = np.linspace(0.5, 5.5, X2)
    y2 = np.linspace(0.5, 5.5, Y2)
    xx2, yy2 = np.meshgrid(x2, y2, sparse=False)
    zz2 = np.sin(xx2 ** 2 + yy2 ** 2) / (xx2 ** 2 + yy2 ** 2)
    data2 = np.sin(xx2 ** 2 + 2 * yy2 ** 2) / (xx2 ** 2 + yy2 ** 2)

    # データセット1のメッシュデータと色の取得
    mesh_ds1, face_values1, vertices1, faces1 = create_mesh_data(
        xx1, yy1, zz1, data1)

    # データセット2のメッシュデータと色の取得
    mesh_ds2, face_values2, vertices2, faces2 = create_mesh_data(
        xx2, yy2, zz2, data2)

    # データセットを連結
    vertices = np.concatenate((vertices1, vertices2))
    faces = np.concatenate(
        (faces1, faces2 + len(vertices1)))  # faces2のインデックスを調整
    face_values = np.concatenate((face_values1, face_values2))

    # Noneを含むデータを処理するために、Noneを最小値に置き換える
    face_values[np.isnan(face_values)] = np.min(
        face_values[np.isfinite(face_values)])

    # メッシュデータソースの作成
    mesh_ds = MeshDS_DataSource(vertices, faces)

    # 面の色の計算
    cmap = colormaps["jet"]

    # カラーマップの範囲を制限するために、データ範囲を設定
    data_min = max(np.min(face_values), -10)  # 最小値を-10に制限（例）
    data_max = min(np.max(face_values), 10)   # 最大値を10に制限（例）
    data_ptp = data_max - data_min

    face_colors = [cmap((value - data_min) / data_ptp)[:3]
                   for value in face_values]

    # メッシュの視覚化の設定
    mesh_vs = MeshVS_Mesh()
    mesh_vs.SetDataSource(mesh_ds)

    # 面の色を設定するビルダーの作成
    element_builder = MeshVS_ElementalColorPrsBuilder(
        mesh_vs, MeshVS_DMF_ElementalColorDataPrs | MeshVS_DMF_OCCMask
    )

    # 面に色を設定
    for nFace in range(faces.shape[0]):
        if np.all(faces[nFace] < vertices.shape[0]):  # インデックスが範囲内か確認
            color = Quantity_Color(*face_colors[nFace], Quantity_TOC_RGB)
            element_builder.SetColor1(nFace + 1, color)  # 面のインデックスは1から始まる

    # ビルダーをメッシュに追加
    mesh_vs.AddBuilder(element_builder, True)

    # ノードの色を設定するビルダーの作成
    aColorMap = Aspect_SequenceOfColor()
    aColorMap.Append(Quantity_Color(Quantity_NOC_PURPLE))
    aColorMap.Append(Quantity_Color(Quantity_NOC_BLUE1))
    aColorMap.Append(Quantity_Color(Quantity_NOC_GREEN))
    aColorMap.Append(Quantity_Color(Quantity_NOC_ORANGE))

    aScaleMap = TColStd_DataMapOfIntegerReal()

    node_builder = MeshVS_NodalColorPrsBuilder(
        mesh_vs, MeshVS_DMF_NodalColorDataPrs | MeshVS_DMF_OCCMask
    )
    node_builder.UseTexture(True)
    node_builder.SetColorMap(aColorMap)
    node_builder.SetInvalidColor(Quantity_Color(Quantity_NOC_BLACK))
    node_builder.SetTextureCoords(aScaleMap)

    # ノード表示のON/OFFを切り替える関数
    display_nodes = [False]

    def toggle_nodes():
        display_nodes[0] = not display_nodes[0]
        if display_nodes[0]:
            mesh_vs.AddBuilder(node_builder, True)
        else:
            # node_builderのインデックスを計算して削除
            builder_index = mesh_vs.GetBuilder(1)
            mesh_vs.RemoveBuilder(1)
        display.Context.UpdateCurrentViewer()

    # add_menu("Node Display")
    # add_function_to_menu("Node Display", toggle_nodes)

    # エッジを非表示にする設定
    mesh_drawer = mesh_vs.GetDrawer()
    mesh_drawer.SetBoolean(MeshVS_DA_ShowEdges, False)
    mesh_vs.SetDrawer(mesh_drawer)

    # メッシュの表示
    display.Context.Display(mesh_vs, True)
    display.FitAll()
    start_display()
