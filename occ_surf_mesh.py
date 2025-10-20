from collections.abc import Sequence

import numpy as np
from OCC.Core.BRep import BRep_Builder, BRep_Tool
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Common
from OCC.Core.BRepMesh import BRepMesh_IncrementalMesh
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCC.Core.BRepTools import breptools
from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_REVERSED
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Core.TopoDS import TopoDS_Shape


def load_brep(file_path: str) -> TopoDS_Shape:
    shape = TopoDS_Shape()
    builder = BRep_Builder()
    breptools.Read(shape, file_path, builder)
    return shape


def show_shapes(shapes: Sequence[TopoDS_Shape]):
    from OCC.Display.SimpleGui import init_display

    display, start_display, _add_menu, _add_function_to_menu = init_display()

    for shape in shapes:
        display.DisplayShape(shape, update=True, transparency=0.5)
    start_display()


def intersect_with_box(shape: TopoDS_Shape):
    box = BRepPrimAPI_MakeBox(10, 10, 10).Shape()

    intersection = BRepAlgoAPI_Common(shape, box)
    intersection.Build()

    if intersection.IsDone():
        return intersection.Shape()
    else:
        raise RuntimeError("Intersection failed")


def write_stl(shape: TopoDS_Shape, filename: str) -> None:
    from OCC.Core.StlAPI import StlAPI_Writer

    writer = StlAPI_Writer()
    writer.SetASCIIMode(True)
    writer.Write(shape, filename)


def mesh_shape(shape: TopoDS_Shape):
    mesh = BRepMesh_IncrementalMesh(shape, 0.1)
    mesh.Perform()

    nodes: list[list[float]] = []
    vertices: list[Sequence[int]] = []

    # https://github.com/Open-Cascade-SAS/OCCT/blob/508700117cf4d41b99087deed2b05f93e751e5cf/src/DataExchange/TKDESTL/StlAPI/StlAPI_Writer.cxx#L69-L116
    exp = TopExp_Explorer(shape, TopAbs_FACE)
    while exp.More():
        face = exp.Current()
        if face.Orientation() == TopAbs_REVERSED:
            print("FIXME Reversed face")

        loc = TopLoc_Location()
        tri = BRep_Tool.Triangulation(face, loc)
        if tri:
            trsf = loc.Transformation()
            for i in range(tri.NbNodes()):
                p = tri.Node(1 + i)
                p.Transform(trsf)
                nodes.append([p.X(), p.Y(), p.Z()])

            vertices.extend(tri.Triangle(1 + j).Get() for j in range(tri.NbTriangles()))

        exp.Next()

    return np.array(nodes), np.array(vertices) - 1


def main():
    # shape = load_brep("sphere.brep")
    # intersection_shape = intersect_with_box(shape)
    # print(intersection_shape.DumpJson())

    box = BRepPrimAPI_MakeBox(10, 10, 10).Shape()

    nodes, connectivity = mesh_shape(box)

    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    ax = plt.figure().add_subplot(projection="3d")

    ax.add_collection3d(
        Poly3DCollection(
            [nodes[tri] for tri in connectivity],
            facecolors="skyblue",
            edgecolors="k",
            linewidths=0.5,
            alpha=0.8,
        )
    )

    plt.show()


if __name__ == "__main__":
    # https://github.com/tpaviot/pythonocc-core/issues/1444
    main()
