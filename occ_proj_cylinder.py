import numpy as np

from src.base_occ import dispocc, set_loc

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_XYZ
from OCC.Core.gp import gp_Lin, gp_Elips, gp_Pln
from OCC.Core.gp import gp_Mat, gp_GTrsf, gp_Trsf
from OCC.Core.gp import gp_Cylinder
from OCC.Core.TColgp import TColgp_Array1OfPnt, TColgp_Array2OfPnt
from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.Core.BRepProj import BRepProj_Projection
from OCC.Core.BRepAlgo import BRepAlgo_FaceRestrictor, BRepAlgo_NormalProjection
from OCC.Core.ProjLib import ProjLib_ProjectOnSurface
from OCC.Core.GeomAPI import GeomAPI_PointsToBSplineSurface
from OCC.Core.GeomAbs import GeomAbs_G1, GeomAbs_G2
from OCC.Extend.DataExchange import write_step_file
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCCUtils.Construct import make_face, make_polygon, make_plane, make_wire


def spl_face(px, py, pz):
    nx, ny = px.shape
    pnt_2d = TColgp_Array2OfPnt(1, nx, 1, ny)
    for row in range(pnt_2d.LowerRow(), pnt_2d.UpperRow() + 1):
        for col in range(pnt_2d.LowerCol(), pnt_2d.UpperCol() + 1):
            i, j = row - 1, col - 1
            pnt = gp_Pnt(px[i, j], py[i, j], pz[i, j])
            pnt_2d.SetValue(row, col, pnt)
            #print (i, j, px[i, j], py[i, j], pz[i, j])

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


if __name__ == '__main__':

    obj = dispocc()

    # make polygon on XY-plane
    pts = [
        gp_Pnt(-50, -50, 0),
        gp_Pnt(50, -50, 0),
        gp_Pnt(50, 50, 0),
        gp_Pnt(0, 60, 0),
        gp_Pnt(-50, 50, 0)
    ]
    axs = gp_Ax3(
        gp_Pnt(0, -200, 0),
        gp_Dir(0,0,1)
    )
    #axs.SetDirection(gp_Dir(0,1,0))
    poly_2d = make_polygon(pts, True)
    poly_2d.Move(set_loc(axs, gp_Ax3()))
    obj.display.DisplayShape(poly_2d)
    
    surf = make_face(gp_Cylinder(gp_Ax3(), 100),0,2*np.pi,0, 100)
    obj.display.DisplayShape(surf)
    
    obj.ShowOCC()
