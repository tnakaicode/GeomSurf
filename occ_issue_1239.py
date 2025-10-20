import numpy as np

from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeCylinder, BRepPrimAPI_MakeBox
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
from OCC.Core.gp import gp_Trsf, gp_Vec, gp_Ax1, gp_Pnt, gp_Dir
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Extend.DataExchange import write_iges_file, write_step_file

from OCC.Display.SimpleGui import init_display

cil01 = BRepPrimAPI_MakeCylinder(10, 10).Shape()
cil02 = BRepPrimAPI_MakeCylinder(9, 16).Shape()
box01 = BRepPrimAPI_MakeBox(20, 20, 20).Shape()

transform = gp_Trsf()
transform.SetTranslationPart(gp_Vec(-10, 0, 0))
location = TopLoc_Location(transform)
box01.Location(location)

transform.SetTranslationPart(gp_Vec(0, 0, -3))
location = TopLoc_Location(transform)
cil02.Location(location)

solid01 = BRepAlgoAPI_Cut(cil01, cil02).Shape()
arc1 = BRepAlgoAPI_Cut(solid01, box01).Shape()
cil11 = BRepPrimAPI_MakeCylinder(10, 10).Shape()
cil12 = BRepPrimAPI_MakeCylinder(9, 16).Shape()
box11 = BRepPrimAPI_MakeBox(20, 20, 20).Shape()

transform = gp_Trsf()
transform.SetTranslationPart(gp_Vec(-10, 0, 0))
location = TopLoc_Location(transform)
box11.Location(location)

transform.SetTranslationPart(gp_Vec(0, 0, -3))
location = TopLoc_Location(transform)
cil12.Location(location)

solid11 = BRepAlgoAPI_Cut(cil11, cil12).Shape()
arc2 = BRepAlgoAPI_Cut(solid11, box11).Shape()
cil21 = BRepPrimAPI_MakeCylinder(10, 10).Shape()
cil22 = BRepPrimAPI_MakeCylinder(9, 16).Shape()
box21 = BRepPrimAPI_MakeBox(20, 20, 20).Shape()

transform = gp_Trsf()
transform.SetTranslationPart(gp_Vec(-10, 0, 0))
location = TopLoc_Location(transform)
box21.Location(location)

transform.SetTranslationPart(gp_Vec(0, 0, -3))
location = TopLoc_Location(transform)
cil22.Location(location)

solid21 = BRepAlgoAPI_Cut(cil21, cil22).Shape()
arc3 = BRepAlgoAPI_Cut(solid21, box21).Shape()
cil31 = BRepPrimAPI_MakeCylinder(10, 10).Shape()
cil32 = BRepPrimAPI_MakeCylinder(9, 16).Shape()
box31 = BRepPrimAPI_MakeBox(20, 20, 20).Shape()
box32 = BRepPrimAPI_MakeBox(10, 10, 10).Shape()

transform.SetTranslationPart(gp_Vec(-10, 0, 0))
location = TopLoc_Location(transform)
box31.Location(location)

transform.SetTranslationPart(gp_Vec(-10, -10, 0))
location = TopLoc_Location(transform)
box32.Location(location)

transform.SetTranslationPart(gp_Vec(0, 0, -3))
location = TopLoc_Location(transform)
cil32.Location(location)

solid31 = BRepAlgoAPI_Cut(cil31, cil32).Shape()
halfarc31 = BRepAlgoAPI_Cut(solid31, box31).Shape()
arc4 = BRepAlgoAPI_Cut(halfarc31, box32).Shape()
cil41 = BRepPrimAPI_MakeCylinder(10, 10).Shape()
cil42 = BRepPrimAPI_MakeCylinder(9, 16).Shape()
box41 = BRepPrimAPI_MakeBox(20, 20, 20).Shape()
box42 = BRepPrimAPI_MakeBox(10, 10, 10).Shape()

transform.SetTranslationPart(gp_Vec(-10, 0, 0))
location = TopLoc_Location(transform)
box41.Location(location)

transform.SetTranslationPart(gp_Vec(-10, -10, 0))
location = TopLoc_Location(transform)
box42.Location(location)

transform.SetTranslationPart(gp_Vec(0, 0, -3))
location = TopLoc_Location(transform)
cil42.Location(location)

solid41 = BRepAlgoAPI_Cut(cil41, cil42).Shape()
halfarc41 = BRepAlgoAPI_Cut(solid41, box41).Shape()
arc5 = BRepAlgoAPI_Cut(halfarc41, box42).Shape()
cil51 = BRepPrimAPI_MakeCylinder(10, 10).Shape()
cil52 = BRepPrimAPI_MakeCylinder(9, 16).Shape()
box51 = BRepPrimAPI_MakeBox(20, 20, 20).Shape()
box52 = BRepPrimAPI_MakeBox(10, 10, 10).Shape()

transform.SetTranslationPart(gp_Vec(-10, 0, 0))
location = TopLoc_Location(transform)
box51.Location(location)

transform.SetTranslationPart(gp_Vec(-10, -10, 0))
location = TopLoc_Location(transform)
box52.Location(location)

transform.SetTranslationPart(gp_Vec(0, 0, -3))
location = TopLoc_Location(transform)
cil52.Location(location)

solid51 = BRepAlgoAPI_Cut(cil51, cil52).Shape()
halfarc51 = BRepAlgoAPI_Cut(solid51, box51).Shape()
arc6 = BRepAlgoAPI_Cut(halfarc51, box52).Shape()
box101 = BRepPrimAPI_MakeBox(100, 10, 10).Shape()
box102 = BRepPrimAPI_MakeBox(10, 30, 10).Shape()

transform = gp_Trsf()
angle = np.deg2rad(120)
axis = gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1))

transform.SetRotation(axis, angle)
transform.SetTranslationPart(gp_Vec(52, -20, 0))
location = TopLoc_Location(transform)
box101.Location(location)

angle = np.deg2rad(0)
axis = gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1))
transform.SetRotation(axis, angle)
transform.SetTranslationPart(gp_Vec(0, 50, 0))
location = TopLoc_Location(transform)
box102.Location(location)

box1 = BRepAlgoAPI_Fuse(box101, box102).Shape()
box111 = BRepPrimAPI_MakeBox(100, 10, 10).Shape()
box112 = BRepPrimAPI_MakeBox(10, 30, 10).Shape()

transform = gp_Trsf()
angle = np.deg2rad(120)
axis = gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1))
transform.SetRotation(axis, angle)
transform.SetTranslationPart(gp_Vec(52, -20, 0))
location = TopLoc_Location(transform)
box111.Location(location)

angle = np.deg2rad(0)
axis = gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1))
transform.SetRotation(axis, angle)
transform.SetTranslationPart(gp_Vec(0, 50, 0))
location = TopLoc_Location(transform)
box112.Location(location)

box2 = BRepAlgoAPI_Fuse(box111, box112).Shape()
box = BRepPrimAPI_MakeBox(2, 50, 20).Shape()

transform = gp_Trsf()
angle = np.deg2rad(-180)
axis = gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0))
transform.SetRotation(axis, angle)
transform.SetTranslationPart(gp_Vec(19, 0, 10))
location = TopLoc_Location(transform)
arc2.Location(location)

angle = np.deg2rad(-180)
axis = gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0))
transform.SetRotation(axis, angle)
transform.SetTranslationPart(gp_Vec(-19, 0, 10))
location = TopLoc_Location(transform)
arc3.Location(location)

angle = np.deg2rad(0)
axis = gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1))
transform.SetRotation(axis, angle)
transform.SetTranslationPart(gp_Vec(-38, 0, 0))
location = TopLoc_Location(transform)
arc4.Location(location)

angle = np.deg2rad(-90)
axis = gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1))
transform.SetRotation(axis, angle)
transform.SetTranslationPart(gp_Vec(38, 0, 0))
location = TopLoc_Location(transform)
arc5.Location(location)

angle = np.deg2rad(-45)
axis = gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1))
transform.SetRotation(axis, angle)
transform.SetTranslationPart(gp_Vec(0, 50, 0))
location = TopLoc_Location(transform)
arc6.Location(location)

angle = np.deg2rad(180)
axis = gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 1, 0))
transform.SetRotation(axis, angle)
transform.SetTranslationPart(gp_Vec(0, 0, 10))
location = TopLoc_Location(transform)
box2.Location(location)

angle = np.deg2rad(0)
axis = gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0))
transform.SetRotation(axis, angle)
transform.SetTranslationPart(gp_Vec(-1, 45, -5))
location = TopLoc_Location(transform)
box.Location(location)

flame = BRepAlgoAPI_Fuse(box1, box2).Shape()
flame2 = BRepAlgoAPI_Cut(flame, box).Shape()
cp0 = BRepAlgoAPI_Fuse(flame2, arc1).Shape()
cp1 = BRepAlgoAPI_Fuse(cp0, arc2).Shape()
cp2 = BRepAlgoAPI_Fuse(cp1, arc3).Shape()
cp3 = BRepAlgoAPI_Fuse(cp2, arc4).Shape()
body = BRepAlgoAPI_Fuse(cp3, arc5).Shape()
picker = BRepAlgoAPI_Fuse(body, arc6).Shape()

write_step_file(picker, 'occ_issue_1239.stp')
write_iges_file(picker, 'occ_issue_1239.iges')

occ, start_display, add_menu, add_function_to_menu = init_display()
occ.DisplayShape(picker, update=True)
start_display()
# https://github.com/tpaviot/pythonocc-core/issues/1239
