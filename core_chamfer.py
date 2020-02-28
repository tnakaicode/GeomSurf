"""
#include <TopoDS_Shape.hxx> 
#include <TopoDS.hxx> 
#include <BRepPrimAPI_MakeBox.hxx> 
#include <TopoDS_Solid.hxx> 
#include <BRepFilletAPI_MakeFillet.hxx> 
#include <TopExp_Explorer.hxx> 
TopoDS_Shape FilletedBox(const Standard_Real a, 
                        const Standard_Real  b, 
                        const Standard_Real  c, 
                        const Standard_Real  r) 
{ 
    TopoDS_Solid Box =  BRepPrimAPI_MakeBox(a,b,c); 
    BRepFilletAPI_MakeFillet  MF(Box); 
    
    // add all the edges  to fillet 
    TopExp_Explorer  ex(Box,TopAbs_EDGE); 
    while (ex.More()) 
    { 
    MF.Add(r,TopoDS::Edge(ex.Current())); 
    ex.Next(); 
    } 
    return MF.Shape(); 
} 
void CSampleTopologicalOperationsDoc::OnEvolvedblend1() 
{ 
    TopoDS_Shape theBox  = BRepPrimAPI_MakeBox(200,200,200); 
    BRepFilletAPI_MakeFillet  Rake(theBox); 
    ChFi3d_FilletShape  FSh = ChFi3d_Rational; 
    Rake.SetFilletShape(FSh); 
    TColgp_Array1OfPnt2d  ParAndRad(1, 6); 
    ParAndRad(1).SetCoord(0.,  10.); 
    ParAndRad(1).SetCoord(50.,  20.); 
    ParAndRad(1).SetCoord(70.,  20.); 
    ParAndRad(1).SetCoord(130.,  60.); 
    ParAndRad(1).SetCoord(160.,  30.); 
    ParAndRad(1).SetCoord(200.,  20.); 
    TopExp_Explorer  ex(theBox,TopAbs_EDGE); 
    Rake.Add(ParAndRad, TopoDS::Edge(ex.Current())); 
    TopoDS_Shape  evolvedBox = Rake.Shape(); 
} 
"""

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
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCC.Core.BRepFilletAPI import BRepFilletAPI_MakeChamfer, BRepFilletAPI_MakeFillet
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
    box = make_box(100, 100, 100)

    display, start_display, add_menu, add_functionto_menu = init_display()

    display.DisplayShape(box)
    display.DisplayShape(axs.Location())

    display.FitAll()
    start_display()
