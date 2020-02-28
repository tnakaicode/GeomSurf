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


def create_tempdir(flag=1):
    print(datetime.date.today())
    datenm = "{0:%Y%m%d}".format(datetime.date.today())
    dirnum = len(glob.glob("./temp_" + datenm + "*/"))
    if flag == -1 or dirnum == 0:
        tmpdir = "./temp_{}{:03}/".format(datenm, dirnum)
        os.makedirs(tmpdir)
        fp = open(tmpdir + "not_ignore.txt", "w")
        fp.close()
    else:
        tmpdir = "./temp_{}{:03}/".format(datenm, dirnum - 1)
    print(tmpdir)
    return tmpdir


if __name__ == '__main__':
    argvs = sys.argv
    parser = OptionParser()
    opt, argc = parser.parse_args(argvs)
    print(opt, argc)
    tmpdir = create_tempdir(1)

    #
    # https://www.opencascade.com/doc/occt-6.9.1/overview/html/occt_user_guides__modeling_algos.html#occt_modalg_6
    #

    axs = gp_Ax3()
    box = make_box(200, 200, 200)
    chf = BRepFilletAPI_MakeChamfer(box)
    # chf.Build()
    fil = BRepFilletAPI_MakeFillet(box)
    fil.SetFilletShape(ChFi3d_Rational)
    par = TColgp_Array1OfPnt2d(1, 6)
    par.SetValue(1, gp_Pnt2d(0, 10))
    par.SetValue(2, gp_Pnt2d(50, 20))
    par.SetValue(3, gp_Pnt2d(70, 20))
    par.SetValue(4, gp_Pnt2d(130, 60))
    par.SetValue(5, gp_Pnt2d(160, 30))
    par.SetValue(6, gp_Pnt2d(200, 20))
    top = TopExp_Explorer(box, TopAbs_EDGE)

    fil.Add(par, top.Current())
    top.Next()
    fil.Add(par, top.Current())

    write_step_file(box, tmpdir + "box.stp")
    write_step_file(fil.Shape(), tmpdir + "box_fillet.stp", "AP214IS")

    display, start_display, add_menu, add_functionto_menu = init_display()

    display.DisplayShape(fil.Shape())
    display.DisplayShape(axs.Location())

    display.FitAll()
    start_display()
