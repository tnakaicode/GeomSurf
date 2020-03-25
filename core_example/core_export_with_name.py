import numpy as np
import sys
import time
import os

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pnt
from OCC.Core.XSControl import XSControl_Writer, XSControl_WorkSession
from OCC.Core.XCAFApp import XCAFApp_Application
from OCC.Core.XCAFDoc import XCAFDoc_DocumentTool_ShapeTool
from OCC.Core.STEPCAFControl import STEPCAFControl_Writer
from OCC.Core.STEPControl import STEPControl_Writer, STEPControl_Reader
from OCC.Core.STEPControl import STEPControl_AsIs
from OCC.Core.Interface import Interface_Static_SetCVal
from OCC.Core.IFSelect import IFSelect_RetDone
from OCC.Core.TDF import TDF_LabelSequence, TDF_Label, TDF_Tool, TDF_Data
from OCC.Core.TDocStd import TDocStd_Document
from OCC.Core.TDataStd import TDataStd_Name, TDataStd_Name_GetID
from OCC.Core.TCollection import TCollection_AsciiString
from OCC.Core.TCollection import TCollection_ExtendedString
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCC.Extend.DataExchange import write_step_file, read_step_file
from OCC.Extend.DataExchange import read_step_file_with_names_colors
from OCCUtils.Construct import make_plane, make_vertex, make_circle


# https://www.opencascade.com/doc/occt-7.4.0/refman/html/class_t_collection___extended_string.html
# https://www.opencascade.com/doc/occt-7.4.0/refman/html/class_x_c_a_f_app___application.html
# https://www.opencascade.com/doc/occt-7.4.0/refman/html/class_t_doc_std___document.html
class ExportCAFMethod (object):

    def __init__(self, name="name", tol=1.0E-10):
        self.name = name
        self.writer = STEPCAFControl_Writer()
        self.writer.SetNameMode(True)
        self.doc = TDocStd_Document(
            TCollection_ExtendedString("pythonocc-doc"))
        #self.x_app = XCAFApp_Application.GetApplication()
        # self.x_app.NewDocument(
        #    TCollection_ExtendedString("MDTV-CAF"), self.doc)
        self.shape_tool = XCAFDoc_DocumentTool_ShapeTool(self.doc.Main())
        Interface_Static_SetCVal("write.step.schema", "AP214")

    def Add(self, shape, name="name"):
        """
        STEPControl_AsIs                   translates an Open CASCADE shape to its highest possible STEP representation.
        STEPControl_ManifoldSolidBrep      translates an Open CASCADE shape to a STEP manifold_solid_brep or brep_with_voids entity.
        STEPControl_FacetedBrep            translates an Open CASCADE shape into a STEP faceted_brep entity.
        STEPControl_ShellBasedSurfaceModel translates an Open CASCADE shape into a STEP shell_based_surface_model entity.
        STEPControl_GeometricCurveSet      translates an Open CASCADE shape into a STEP geometric_curve_set entity.
        """
        label = self.shape_tool.AddShape(shape)
        self.writer.Transfer(self.doc, STEPControl_AsIs)

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

    box = BRepPrimAPI_MakeBox(10, 10, 10).Solid()
    writer = STEPControl_Writer()
    fp = writer.WS().TransferWriter().FinderProcess()

    print(fp)

    # start_display()
