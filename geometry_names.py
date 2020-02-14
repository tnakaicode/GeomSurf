import numpy as np
from OCC.Core.STEPControl import STEPControl_Reader
from OCC.Core.TopAbs import TopAbs_FACE
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.StepRepr import StepRepr_RepresentationItem
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
from OCC.Core.GeomAbs import GeomAbs_Cylinder
from OCC.Core.TopoDS import topods

# Read the file and get the shape
reader = STEPControl_Reader()
tr = reader.WS().TransferReader()
reader.ReadFile('geometry_names.stp')
reader.TransferRoots()
shape = reader.OneShape()

# Explore the faces of the shape (these are known to be named)
exp = TopExp_Explorer(shape, TopAbs_FACE)
while exp.More():
    s = exp.Current()
    exp.Next()

    # Converting TopoDS_Shape to TopoDS_Face
    face = topods.Face(s)

    # Working on the geometry
    face_adaptor = BRepAdaptor_Surface(face)
    face_type = face_adaptor.GetType()
    if face_type == GeomAbs_Cylinder:
        cylinder = face_adaptor.Cylinder()
        entity = {}
        entity['type'] = 'cylinder'
        entity['location'] = cylinder.Axis().Location().Coord()
        entity['direction'] = cylinder.Axis().Direction().Coord()
        entity['radius'] = cylinder.Radius()
        entity['coefficients'] = cylinder.Coefficients()
        print('cylinder:', entity)

    # Working on the name
    item = tr.EntityFromShapeResult(s, 1)
    # if item.IsNull():
    #    continue
    item = StepRepr_RepresentationItem.DownCast(item)
    name = item.Name().ToCString()
    if name:
        print('Found entity named: {}: {}.'.format(name, s))
