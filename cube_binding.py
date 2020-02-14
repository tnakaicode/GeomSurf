# https://gist.github.com/jay3sh/1351019
import sys

from OCC.Core.gp import gp_Pnt
from OCC.Core.GC import GC_MakeSegment
from OCC.Core.BRepBuilderAPI import \
    BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeWire,\
    BRepBuilderAPI_MakeShell, BRepBuilderAPI_MakeSolid
from OCC.Core.BRep import BRep_Builder
from OCC.Core.TopoDS import TopoDS_Shell, TopoDS_Solid
from OCC.Core.StlAPI import StlAPI_Writer
from OCC.Extend.DataExchange import write_stl_file

mesh = {
    "vertices": [[-0.2, -0.2, 0.2], [0.2, -0.2, 0.2], [0.2, 0.2, 0.2], [-0.2, 0.2, 0.2], [-0.2, -0.2, 0.6000000000000001], [0.2, -0.2, 0.6000000000000001], [0.2, 0.2, 0.6000000000000001], [-0.2, 0.2, 0.6000000000000001]],
    "faces": [[3, 2, 1, 0], [4, 5, 6, 7], [7, 6, 2, 3], [5, 4, 0, 1], [6, 5, 1, 2], [4, 7, 3, 0]]
}


def main():

    vertices = [gp_Pnt(p[0], p[1], p[2]) for p in mesh['vertices']]
    oFaces = []

    builder = BRep_Builder()
    shell = TopoDS_Shell()
    builder.MakeShell(shell)

    for face in mesh['faces']:
        edges = []
        face.reverse()
        for i in range(len(face)):
            cur = face[i]
            nxt = face[(i + 1) % len(face)]
            segment = GC_MakeSegment(vertices[cur], vertices[nxt])
            edges.append(BRepBuilderAPI_MakeEdge(segment.Value()))

        wire = BRepBuilderAPI_MakeWire()
        for edge in edges:
            wire.Add(edge.Edge())

        oFace = BRepBuilderAPI_MakeFace(wire.Wire())

        builder.Add(shell, oFace.Shape())

    write_stl_file(shell, "./cube_binding.stl")


if __name__ == '__main__':
    main()
