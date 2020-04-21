import numpy as np

from OCC.Core.gp import gp_Pnt
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCC.Extend.DataExchange import write_iges_file, write_step_file, write_stl_file
from OCC.Extend.ShapeFactory import make_box


shp = make_box(gp_Pnt(-50, -50, -50), 100, 100, 100)
write_iges_file(shp, "box.iges")
write_step_file(shp, "box.stp")
write_stl_file(shp, "box.stl")
