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
import argparse

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Pnt2d
from OCC.Core.MeshVS import MeshVS_Mesh, MeshVS_PrsBuilder
from OCC.Core.IMeshTools import IMeshTools_MeshAlgoFactory
from OCCUtils.Construct import make_box
from OCCUtils.Construct import make_line, make_wire, make_edge

from src.base_occ import dispocc, write_stl_file_mesh1, write_stl_file_mesh2


if __name__ == '__main__':
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    opt = parser.parse_args()
    print(opt, argvs)

    obj = dispocc()
    axs = gp_Ax3()
    box = make_box(axs.Ax2(), 200, 200, 200)
    obj.display.DisplayShape(box)
    obj.show_axs_pln(axs, scale=250)
    obj.ShowOCC()
