import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import json
import sys
import time
import os
import glob
import shutil
import datetime
from optparse import OptionParser

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Pnt2d
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCC.Core.BRepFilletAPI import BRepFilletAPI_MakeChamfer, BRepFilletAPI_MakeFillet
from OCC.Core.ChFi3d import ChFi3d_Rational
from OCC.Core.TColgp import TColgp_Array1OfPnt2d
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_EDGE
from OCC.Extend.DataExchange import read_step_file, write_step_file
from OCCUtils.Construct import make_box
from OCCUtils.Construct import make_line, make_wire, make_edge


from src.base import plotocc, create_tempdir


if __name__ == '__main__':
    argvs = sys.argv
    parser = OptionParser()
    opt, argc = parser.parse_args(argvs)
    print(opt, argc)

    #
    # https://www.opencascade.com/doc/occt-6.9.1/overview/html/occt_user_guides__modeling_algos.html#occt_modalg_6
    #

    obj = plotocc()

    axs = gp_Ax3()
    box = make_box(200, 200, 200)
    chf = BRepFilletAPI_MakeChamfer(box)
    # chf.Build()
    fil = BRepFilletAPI_MakeFillet(box)
    fil.SetFilletShape(ChFi3d_Rational)
    par = TColgp_Array1OfPnt2d(1, 2)
    par.SetValue(1, gp_Pnt2d(-1000, 10))
    par.SetValue(2, gp_Pnt2d(1000, 10))
    top = TopExp_Explorer(box, TopAbs_EDGE)

    fil.Add(par, top.Current())
    top.Next()
    fil.Add(par, top.Current())
    top.Next()
    fil.Add(par, top.Current())

    write_step_file(box, obj.tmpdir + "box.stp")
    write_step_file(fil.Shape(), obj.tmpdir + "box_fillet.stp", "AP214IS")

    obj.display.DisplayShape(fil.Shape())
    obj.display.DisplayShape(axs.Location())
    obj.show()
