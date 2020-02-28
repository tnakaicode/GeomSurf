#!/usr/bin/env python
##
# This file is part of pythonOCC.
##
# pythonOCC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
##
# pythonOCC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
##
# You should have received a copy of the GNU Lesser General Public License
# along with pythonOCC.  If not, see <http://www.gnu.org/licenses/>.

import random

from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCC.Core.gp import gp_Pnt, gp_Ax2, gp_Dir, gp_XYZ
from OCC.Core.BRepBndLib import brepbndlib_AddOBB
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeVertex
from OCC.Core.Bnd import Bnd_OBB
from OCC.Extend.DataExchange import write_step_file

from OCC.Display.SimpleGui import init_display

display, start_display, add_menu, add_function_to_menu = init_display()


def ConvertBndToShape(theBox):
    aBaryCenter = theBox.Center()
    aXDir = theBox.XDirection()
    aYDir = theBox.YDirection()
    aZDir = theBox.ZDirection()
    aHalfX = theBox.XHSize()
    aHalfY = theBox.YHSize()
    aHalfZ = theBox.ZHSize()

    ax = gp_XYZ(aXDir.X(), aXDir.Y(), aXDir.Z())
    ay = gp_XYZ(aYDir.X(), aYDir.Y(), aYDir.Z())
    az = gp_XYZ(aZDir.X(), aZDir.Y(), aZDir.Z())
    p = gp_Pnt(aBaryCenter.X(), aBaryCenter.Y(), aBaryCenter.Z())
    anAxes = gp_Ax2(p, gp_Dir(aZDir), gp_Dir(aXDir))
    anAxes.SetLocation(
        gp_Pnt(p.XYZ() - ax * aHalfX - ay * aHalfY - az * aHalfZ))
    aBox = BRepPrimAPI_MakeBox(
        anAxes, 2.0 * aHalfX, 2.0 * aHalfY, 2.0 * aHalfZ).Shape()
    return aBox


obb = Bnd_OBB()

# choose n random vertices
n = 10
for _ in range(n):
    x = random.uniform(100, 1000)
    y = random.uniform(100, 1000)
    z = random.uniform(100, 1000)
    p = BRepBuilderAPI_MakeVertex(gp_Pnt(x, y, z)).Shape()
    display.DisplayShape(p)
    brepbndlib_AddOBB(p, obb)
obb_shape = ConvertBndToShape(obb)
display.DisplayShape(obb_shape)

# a ref box
b = BRepPrimAPI_MakeBox(10, 10, 10).Shape()
display.DisplayShape(b, update=True)
write_step_file(b, "./tmp/S_Bois.stp")

start_display()
