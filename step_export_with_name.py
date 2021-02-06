import numpy as np
import matplotlib.pyplot as plt
import sys
import time
import os

from OCC.Core.gp import gp_Pnt
from OCC.Core.IFSelect import IFSelect_RetError
from OCC.Core.Interface import Interface_Static_SetCVal
from OCC.Core.STEPConstruct import stepconstruct_FindEntity
from OCC.Core.STEPControl import STEPControl_AsIs, STEPControl_Writer
from OCC.Core.TCollection import TCollection_HAsciiString
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCC.Core.TopoDS import TopoDS_Iterator, topods_Vertex
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Compound
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Core.TColgp import TColgp_Array1OfPnt, TColgp_Array2OfPnt
from OCC.Core.TColgp import TColgp_HArray1OfPnt, TColgp_HArray2OfPnt
from OCC.Core.BRep import BRep_Builder
from OCC.Extend.ShapeFactory import make_vertex
from OCC.Extend.DataExchange import read_step_file_with_names_colors

from src.base import plotocc, spl_face


class Model (plotocc):

    def __init__(self):
        super().__init__()

        self.builder = BRep_Builder()
        self.compound = TopoDS_Compound()
        self.builder.MakeCompound(self.compound)

        schema = 'AP203'
        assembly_mode = 1

        self.writer = STEPControl_Writer()
        self.fp = self.writer.WS().TransferWriter().FinderProcess()
        Interface_Static_SetCVal('write.step.schema', schema)
        Interface_Static_SetCVal('write.step.unit', 'M')
        Interface_Static_SetCVal('write.step.assembly', str(assembly_mode))

    def AddShape(self, shp, name="shape"):
        Interface_Static_SetCVal('write.step.product.name', name)
        status = self.writer.Transfer(shp, STEPControl_AsIs)
        if int(status) > int(IFSelect_RetError):
            raise Exception('Some Error occurred')

        # This portion is not working as I hoped
        item = stepconstruct_FindEntity(self.fp, shp)
        if not item:
            raise Exception('Item not found')

    def export_stp_with_name(self):
        status = self.writer.Write(self.tempname + ".stp")
        if int(status) > int(IFSelect_RetError):
            raise Exception('Something bad happened')


if __name__ == '__main__':
    obj = Model()

    my_box1 = BRepPrimAPI_MakeBox(10., 20., 30.).Shape()
    my_box2 = BRepPrimAPI_MakeBox(20., 1., 30.).Shape()

    px = np.linspace(-1, 1, 5) * 100 + 50
    py = np.linspace(-1, 1, 7) * 200 - 50
    mesh = np.meshgrid(px, py)
    surf = mesh[0]**2 / 1000
    face = spl_face(*mesh, surf)

    obj.AddShape(my_box1, "box1")
    obj.AddShape(my_box2, "box2")
    obj.AddShape(face, "face")
    for idx, val in np.ndenumerate(surf):
        name = "pnt_{:d}_{:d}".format(*idx)
        pnt = gp_Pnt(mesh[0][idx], mesh[1][idx], surf[idx])
        obj.AddShape(make_vertex(pnt), name)
    obj.export_stp_with_name()
