from OCC.Core.IFSelect import IFSelect_RetError
from OCC.Core.Interface import Interface_Static_SetCVal
from OCC.Core.STEPConstruct import stepconstruct_FindEntity
from OCC.Core.STEPControl import (STEPControl_AsIs, STEPControl_Writer)
from OCC.Core.TCollection import TCollection_HAsciiString
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCC.Extend.DataExchange import read_step_file_with_names_colors

schema = 'AP203'
assembly_mode = 1

writer = STEPControl_Writer()
fp = writer.WS().TransferWriter().FinderProcess()
Interface_Static_SetCVal('write.step.schema', schema)
Interface_Static_SetCVal('write.step.unit', 'M')
Interface_Static_SetCVal('write.step.assembly', str(assembly_mode))

my_box1 = BRepPrimAPI_MakeBox(10., 20., 30.).Shape()
my_box2 = BRepPrimAPI_MakeBox(20., 1., 30.).Shape()

components = [my_box1, my_box2]
comp_names = ['PartA', 'PartB']
for i, comp in enumerate(components):
    Interface_Static_SetCVal('write.step.product.name', comp_names[i])
    status = writer.Transfer(comp, STEPControl_AsIs)
    if int(status) > int(IFSelect_RetError):
        raise Exception('Some Error occurred')

    # This portion is not working as I hoped
    item = stepconstruct_FindEntity(fp, comp)
    if not item:
        raise Exception('Item not found')

    item.SetName(TCollection_HAsciiString(comp_names[i]))

status = writer.Write('step_export_with_name.stp')
if int(status) > int(IFSelect_RetError):
    raise Exception('Something bad happened')

read_step_file_with_names_colors('step_export_with_name.stp')
