#!/usr/bin/env python
# coding: utf-8

from OCC.Display.SimpleGui import init_display

import numpy as np
from math import pi

from OCC.Core.gp import gp_Pnt2d, gp_Pnt, gp_Lin2d, gp_Ax3, gp_Dir2d, gp_Dir, gp_Ax2d
from OCC.Core.TColgp import TColgp_Array1OfPnt2d
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge
from OCC.Core.Geom import (
    Geom_CylindricalSurface,
    Geom_ToroidalSurface,
    Geom_Curve,
    Geom_Surface,
    Geom_ConicalSurface,
)
from OCC.Core.Geom2d import Geom2d_Circle, Geom2d_Curve, Geom2d_BSplineCurve
from OCC.Core.Geom2dAPI import Geom2dAPI_Interpolate, Geom2dAPI_PointsToBSpline
from OCC.Core.GCE2d import GCE2d_MakeSegment
from OCCUtils.Common import vertex2pnt
from OCCUtils.Construct import make_face, make_edge, make_wire

display, start_display, add_menu, add_function_to_menu = init_display()


def make_helix_wire(uvcurve=Geom_Curve, surface=Geom_Surface, t0=0, t1=1):
    print(uvcurve, surface)
    helix_api = BRepBuilderAPI_MakeEdge(uvcurve, surface, t0, t1)
    helix_edge = helix_api.Edge()
    helix_vert = [vertex2pnt(helix_api.Vertex1()), vertex2pnt(helix_api.Vertex2())]
    print(helix_vert)
    for t in [t0, (t0 + t1) / 2, t1]:
        u, v = uvcurve.Value(t).X(), uvcurve.Value(t).Y()
        print(t, u, v, surface.Value(u, v))
        # display.DisplayShape(surface.Value(u, v))
    if helix_vert[0].IsEqual(helix_vert[1], 0.1):
        helix_wire = make_wire(helix_edge)
    else:
        helix_wire = make_wire([helix_edge, make_edge(*helix_vert)])
    return helix_wire


aCylinder = Geom_CylindricalSurface(gp_Ax3(), 3.0)
# display.DisplayShape(make_face(aCylinder.Cylinder(), 0, 2 * pi, 0, 10))

# Build an helix
u0, v0 = 0.0, 0.1
ur, vr = 1.0, 1.1
ut, vt = 1 / ur, 1 / vr
uv = np.sqrt(ut**2 + vt**2)
aLine2d = gp_Lin2d(gp_Pnt2d(u0, v0), gp_Dir2d(ut, vt))
aSegment = GCE2d_MakeSegment(aLine2d, 0, 1).Value()

# Cylinder
t0, t1 = 0.0, uv * 2 * pi * ur
helix_wire = make_helix_wire(aSegment, aCylinder, t0, t1)
# face = make_face(aCylinder.Cylinder(), helix_wire)
display.DisplayShape(helix_wire)
# display.DisplayShape(face)

# Build an helix
u0, v0 = 0.0, 0.1
ur, vr = 1.0, -1.1
ut, vt = 1 / ur, 1 / vr
uv = np.sqrt(ut**2 + vt**2)
aLine2d = gp_Lin2d(gp_Pnt2d(u0, v0), gp_Dir2d(ut, vt))
aSegment = GCE2d_MakeSegment(aLine2d, 0, 1).Value()

# Cylinder
t0, t1 = 0.0, -uv * 1.9 * pi * ur
helix_wire = make_helix_wire(aSegment, aCylinder, t0, t1)
face = make_face(aCylinder.Cylinder(), helix_wire)
print(face)
display.DisplayShape(helix_wire)
# display.DisplayShape(face)

# Build an helix
pt = np.linspace(0, 1, 100) * 2 * np.pi
pts = TColgp_Array1OfPnt2d(1, 100)
for i, t in enumerate(pt):
    pts.SetValue(i + 1, gp_Pnt2d(t, np.sin(t) + t))
curve = Geom2dAPI_PointsToBSpline(pts).Curve()
helix_wire = make_helix_wire(curve, aCylinder, 0, 1)
face = make_face(make_face(aCylinder.Cylinder(), 0, 2 * pi, 0, 10), helix_wire)
print(face)
display.DisplayShape(helix_wire, color="BLUE1")
display.DisplayShape(face)

display.FitAll()
start_display()
