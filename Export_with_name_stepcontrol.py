import numpy as np
import sys
import time
import os

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pnt
from OCC.Core.STEPConstruct import stepconstruct_FindEntity
from OCC.Core.STEPControl import STEPControl_Writer
from OCC.Core.STEPControl import STEPControl_AsIs
from OCC.Core.Interface import Interface_Static_SetCVal
from OCC.Core.IFSelect import IFSelect_RetDone, IFSelect_RetError
from OCC.Core.TDataStd import TDataStd_Name, TDataStd_Name_GetID
from OCC.Core.TCollection import TCollection_HAsciiString
from OCC.Core.TCollection import TCollection_AsciiString
from OCC.Core.TCollection import TCollection_ExtendedString
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCC.Extend.DataExchange import write_step_file, read_step_file
from OCCUtils.Construct import make_plane, make_vertex, make_circle


# https://www.opencascade.com/doc/occt-7.4.0/refman/html/class_t_collection___extended_string.html
# https://www.opencascade.com/doc/occt-7.4.0/refman/html/class_x_c_a_f_app___application.html
# https://www.opencascade.com/doc/occt-7.4.0/refman/html/class_t_doc_std___document.html
class ExportCAFMethod (object):

    def __init__(self, name="name", tol=1.0E-10):
        self.name = name
        self.writer = STEPControl_Writer()
        self.fp = self.writer.WS().TransferWriter().FinderProcess()
        Interface_Static_SetCVal("write.step.schema", "AP214")
        Interface_Static_SetCVal('write.step.unit', 'mm')
        Interface_Static_SetCVal('write.step.assembly', str(1))

    def Add(self, shape, name="name"):
        """
        STEPControl_AsIs                   translates an Open CASCADE shape to its highest possible STEP representation.
        STEPControl_ManifoldSolidBrep      translates an Open CASCADE shape to a STEP manifold_solid_brep or brep_with_voids entity.
        STEPControl_FacetedBrep            translates an Open CASCADE shape into a STEP faceted_brep entity.
        STEPControl_ShellBasedSurfaceModel translates an Open CASCADE shape into a STEP shell_based_surface_model entity.
        STEPControl_GeometricCurveSet      translates an Open CASCADE shape into a STEP geometric_curve_set entity.
        """
        Interface_Static_SetCVal('write.step.product.name', name)
        status = self.writer.Transfer(shape, STEPControl_AsIs)
        if int(status) > int(IFSelect_RetError):
            raise Exception('Some Error occurred')

        # This portion is not working as I hoped
        item = stepconstruct_FindEntity(self.fp, shape)
        if not item:
            raise Exception('Item not found')

        item.SetName(TCollection_HAsciiString(name))

    def Write(self, filename=None):
        if not filename:
            filename = self.name
        path, ext = os.path.splitext(filename)
        if not ext:
            ext = ".stp"
        status = self.writer.Write(path + ext)
        assert(status == IFSelect_RetDone)


if __name__ == "__main__":
    display, start_display, add_menu, add_function_to_menu = init_display()
    display.DisplayShape(gp_Pnt())
    root = ExportCAFMethod(name="root")
    root.Add(make_vertex(gp_Pnt()), name="pnt")
    root.Add(make_plane(center=gp_Pnt(0, 0, 0)), name="pln0")
    root.Add(make_plane(center=gp_Pnt(0, 0, 100)), name="pln1")
    root.Add(make_circle(gp_Pnt(0, 0, 0), 100), name="circle")
    root.Write()

    display.FitAll()
    start_display()
