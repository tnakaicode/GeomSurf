import numpy as np

from src.base_occ import dispocc

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_XYZ
from OCC.Core.gp import gp_Lin, gp_Elips, gp_Pln
from OCC.Core.gp import gp_Mat, gp_GTrsf, gp_Trsf
from OCC.Core.gp import gp_Cylinder
from OCC.Core.TColgp import TColgp_Array1OfPnt, TColgp_Array2OfPnt
from OCC.Core.BRep import BRep_Tool
from OCC.Core.TopoDS import TopoDS_Wire, TopoDS_Vertex
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface, BRepAdaptor_Curve
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.Core.BRepProj import BRepProj_Projection
from OCC.Core.BRepAlgo import BRepAlgo_FaceRestrictor, BRepAlgo_NormalProjection
from OCC.Core.BOPAlgo import BOPAlgo_PaveFiller
from OCC.Core.ProjLib import ProjLib_ProjectOnSurface
from OCC.Core.GeomAPI import GeomAPI_PointsToBSplineSurface
from OCC.Core.GeomAbs import GeomAbs_G1, GeomAbs_G2
from OCC.Extend.DataExchange import write_step_file
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCCUtils.Common import vertex2pnt
from OCCUtils.Construct import make_face, make_polygon, make_plane, make_wire


def spl_face(px, py, pz):
    nx, ny = px.shape
    pnt_2d = TColgp_Array2OfPnt(1, nx, 1, ny)
    for row in range(pnt_2d.LowerRow(), pnt_2d.UpperRow() + 1):
        for col in range(pnt_2d.LowerCol(), pnt_2d.UpperCol() + 1):
            i, j = row - 1, col - 1
            pnt = gp_Pnt(px[i, j], py[i, j], pz[i, j])
            pnt_2d.SetValue(row, col, pnt)
            # print (i, j, px[i, j], py[i, j], pz[i, j])

    api = GeomAPI_PointsToBSplineSurface(pnt_2d, 3, 8, GeomAbs_G2, 0.001)
    api.Interpolate(pnt_2d)
    face = BRepBuilderAPI_MakeFace(api.Surface(), 1e-6).Face()
    return face


def make_trimmedface(poly, face, axs=gp_Ax3()):
    """make trimmed face by projected poly
    Args:
        poly (TopoDS_Wire):
        face (TopoDS_Face):
        axs (gp_Ax3, optional): Defaults to gp_Ax3().
    Returns:
        face (TopoDS_Face)
    """
    proj = BRepProj_Projection(poly, face, axs.Direction())
    # make hole
    # face = make_face(face, poly_proj)
    api_face = BRepAlgo_FaceRestrictor()
    api_face.Init(face, True, True)
    api_face.Add(poly)
    api_face.Perform()
    return api_face.Current()


def get_polygon_from_face(face):
    n = 10
    pts = []
    rim_tmp = gp_Pnt()
    for i, e in enumerate(TopologyExplorer(face).edges()):
        e_curve, u0, u1 = BRep_Tool.Curve(e)
        print(i, e, u0, u1, rim_tmp, e_curve.Value(u0), e_curve.Value(u1))
        if i != 0 and rim_tmp == e_curve.Value(u0):
            u_range = np.linspace(u0, u1, n)
            rim_tmp = e_curve.Value(u1)
            p = e_curve.Value(u0)
        elif i != 0 and rim_tmp == e_curve.Value(u1):
            u_range = np.linspace(u1, u0, n)
            rim_tmp = e_curve.Value(u0)
            p = e_curve.Value(u1)
        else:
            u_range = np.linspace(u1, u0, n)
            rim_tmp = e_curve.Value(u0)
            p = e_curve.Value(u1)
        pts.append(p)
        # for u in u_range[1:-1]:
        #    p = e_curve.Value(u)
        #    pts.append(p)
    # return make_polygon(pts, True)
    poly = make_wire([e for e in TopologyExplorer(face).edges()])
    return poly


def get_wire_ends(wire=TopoDS_Wire()):
    edges = [edge for edge in TopologyExplorer(wire).edges()]
    c0 = BRepAdaptor_Curve(edges[0])
    c1 = BRepAdaptor_Curve(edges[-1])
    return [v for v in TopologyExplorer(wire).vertices()]


if __name__ == "__main__":
    obj = dispocc()

    # make polygon on XY-plane
    pts = [
        gp_Pnt(0, -50, 0),
        gp_Pnt(50, -30, 0),
        gp_Pnt(70, -20, 0),
        gp_Pnt(25, 10, 0),
        gp_Pnt(70, 40, 0),
        gp_Pnt(0, 60, 0),
        gp_Pnt(-50, 40, 0),
        gp_Pnt(-45, -20, 0),
    ]
    pts.reverse()
    poly_2d = make_polygon(pts, True)
    obj.display.DisplayShape(poly_2d)

    # make plane on gp_Pnt(0,0,10) which size is over polygon on XY-plane
    plan = make_plane(
        gp_Pnt(0, 0, 0),
        extent_x_min=-60,
        extent_x_max=60,
        extent_y_min=-40,
        extent_y_max=55,
    )
    # obj.display.DisplayShape(plan, transparency=0.7)

    # make curved surface by spline on gp_Pnt(0,0,10) which size is same plane
    px = np.linspace(-60, 60, 100)
    py = np.linspace(-40, 50, 200)
    mesh = np.meshgrid(px, py)
    surf = mesh[0] ** 2 / 1000 + mesh[1] ** 2 / 2000 + 10.0
    face = spl_face(*mesh, surf)
    # face_pnts = get_polygon_from_face(face)
    # face_poly = get_polygon_from_face(face)
    obj.display.DisplayShape(face, transparency=0.9)
    # obj.display.DisplayShape(face_poly, color="GREEN")
    # for i, p in enumerate(get_polygon_from_face(plan)):
    #    display.DisplayShape(p, color="YELLOW")
    #    display.DisplayMessage(p, f"{i:2d}")

    # projection polygon to curved surface
    # the current TopoDS_Wire is not connectted each other beacause the size of curved surface is larger than polygon on XY-Plane.

    proj = BRepProj_Projection(poly_2d, face, gp_Dir(0, 0, 1))
    poly_face = []
    poly_edge = []
    while proj.More():
        poly = proj.Current()
        poly_face.append(poly)
        poly_edge += ([e for e in TopologyExplorer(poly).edges()])
        proj.Next()

    adap = BRepAdaptor_Surface(face)
    u0 = adap.FirstUParameter()
    v0 = adap.FirstVParameter()
    u1 = adap.LastVParameter()
    v1 = adap.LastVParameter()

    poly_plan = make_face(poly_2d)
    face_edge = make_polygon(
        [adap.Value(u, v)
         for u, v in [[u0, v0], [u1, v0], [u1, v1], [u0, v1]]],
        closed=True,
    )
    face_edge = make_wire([e for e in TopologyExplorer(face).edges()])
    proj = BRepProj_Projection(face_edge, make_plane(), gp_Dir(0, 0, 1))
    face_poly = []
    while proj.More():
        poly = proj.Current()
        # proj_face = BRepProj_Projection(poly, face, gp_Dir(0, 0, 1))
        # while proj_face.More():
        #    poly_face.append(proj_face.Current())
        #    proj_face.Next()
        face_poly.append(proj.Current())
        proj.Next()

#    wire = make_wire(poly_face)
    face_plan = make_face(face_poly[0])
    # obj.display.DisplayShape(face_plan)
    # obj.display.DisplayShape(poly_face)
    # obj.display.DisplayShape(poly_plan)

    from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Common
    poly_plan = BRepAlgoAPI_Common(face_plan, poly_plan).Shape()
    poly_comm = make_wire([e for e in TopologyExplorer(poly_plan).edges()])

    proj = BRepProj_Projection(poly_comm, face, gp_Dir(0, 0, 1))
    poly_face = []
    while proj.More():
        poly = proj.Current()
        poly_face.append(poly)
        proj.Next()


#    v0 = [v for v in TopologyExplorer(poly_face[0]).vertices()]
#    poly_face1 = [poly_face[0]]

#    for i, wire in enumerate(poly_face[1:]):
#        # c = BRepAdaptor_Curve(wire)
#        v1 = [v for v in TopologyExplorer(wire).vertices()]
#        print(i, wire, *get_wire_ends(wire))
#        for v in v1:
#            print(
#                vertex2pnt(v),
#                [v.IsEqual(p) for p in v0],
#                [vertex2pnt(v).IsEqual(vertex2pnt(p), 0.1e-3) for p in v0],
#            )
#            # gp_Pnt().IsEqual()
#            # TopoDS_Vertex().IsEqual()

#    idx = [0, 3, 1, 2]  # ok
#    idx = [0, 2, 1, 3]  # ok
#    idx = [0, 3, 2, 1]  # ok
#    idx = [0, 1, 3, 2]  # ng
#    idx = [0, 2, 3, 1]  # ng
#    idx = [0, 1, 2, 3]  # ng
#    poly_face1 = [poly_face[i] for i in idx]
#
    poly_face = make_wire(poly_face)
    trim_face = BRepAlgo_FaceRestrictor()
    trim_face.Init(face, True, True)
    trim_face.Add(poly_face)
    # for poly in poly_face:
    # trim_face.Add(poly)
    trim_face.Perform()
    print(trim_face.IsDone())
    while trim_face.More():
        print(1, trim_face.Current())
        obj.display.DisplayShape(trim_face.Current(), color="RED")
        obj.export_stp(trim_face.Current())
        trim_face.Next()
    obj.display.DisplayShape(poly_face, color="BLUE1")

    obj.ShowOCC()
