from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Cut
from OCC.Core.gp import gp_Pnt, gp_Dir, gp_Ax2
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCC.Core.gp import gp_Trsf
from OCC.Core.STEPControl import STEPControl_Writer, STEPControl_AsIs
from OCC.Core.StlAPI import StlAPI_Writer
from OCC.Display.SimpleGui import init_display

# https://github.com/tpaviot/pythonocc-core/issues/1448

# Start viewer (optional)
display, start_display, add_menu, add_function_to_menu = init_display()

# Parameters
block_length = 60
block_width = 40
block_height = 20

pencil_radius = 4  # Ø8 mm potlood
pencil_depth = 20

offset_arm_length = 80
offset_arm_width = 10
offset_arm_thickness = 6
slot_length = 40
slot_width = 5

bolt_hole_radius = 2.5  # voor M5-schroef
bolt_hole_depth = 10

# 1. Basisblok
block = BRepPrimAPI_MakeBox(block_length, block_width, block_height).Shape()

# 2. Potloodgat (cilinder)
pencil_cyl = BRepPrimAPI_MakeCylinder(
    gp_Ax2(gp_Pnt(10, block_width / 2, block_height / 2), gp_Dir(1, 0, 0)),
    pencil_radius,
    pencil_depth,
).Shape()

block = BRepAlgoAPI_Cut(block, pencil_cyl).Shape()

# 3. Offset-arm
offset_arm = BRepPrimAPI_MakeBox(
    offset_arm_length, offset_arm_width, offset_arm_thickness
).Shape()

# 4. Sleuf in offset-arm (simulate as two small cylinders for now)
slot_cyl_1 = BRepPrimAPI_MakeCylinder(
    gp_Ax2(gp_Pnt(20, offset_arm_width / 2, offset_arm_thickness / 2), gp_Dir(0, 1, 0)),
    bolt_hole_radius,
    offset_arm_width,
).Shape()

slot_cyl_2 = BRepPrimAPI_MakeCylinder(
    gp_Ax2(gp_Pnt(60, offset_arm_width / 2, offset_arm_thickness / 2), gp_Dir(0, 1, 0)),
    bolt_hole_radius,
    offset_arm_width,
).Shape()

offset_arm = BRepAlgoAPI_Cut(offset_arm, slot_cyl_1).Shape()
offset_arm = BRepAlgoAPI_Cut(offset_arm, slot_cyl_2).Shape()

# 5. Transform offset-arm in juiste positie
trsf = gp_Trsf()
trsf.SetTranslation(gp_Pnt(0, 0, block_height), gp_Pnt(10, 10, block_height + 1))
offset_arm_moved = BRepBuilderAPI_Transform(offset_arm, trsf).Shape()

# Toon het model (optioneel)
display.DisplayShape(block, update=True, color="GRAY")
display.DisplayShape(offset_arm_moved, update=True, color="DARKGRAY")

# 6. Exporteren naar STL
writer = StlAPI_Writer()
writer.SetASCIIMode(False)
writer.Write(block, "scribe_jig_base.stl")
writer.Write(offset_arm_moved, "scribe_jig_offset_arm.stl")

print(
    "✅ STL-bestanden zijn geëxporteerd: 'scribe_jig_base.stl' & 'scribe_jig_offset_arm.stl'"
)
start_display()
