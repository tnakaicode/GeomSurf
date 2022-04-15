import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import argparse

sys.path.append(os.path.join("./"))
from src.base_occ import dispocc

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

from OCC.Core.gp import gp_Pnt
from OCC.Core.BRepGProp import brepgprop_SurfaceProperties
from OCC.Core.BRepGProp import brepgprop_VolumeProperties
from OCC.Core.BRepGProp import brepgprop_LinearProperties
from OCC.Core.GProp import GProp_GProps

if __name__ == '__main__':
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument("--pxyz", dest="pxyz",
                      default=[0.0, 0.0, 0.0], type=float, nargs=3)
    opt = parser.parse_args()
    print(opt, argvs)

    prop = GProp_GProps()

    obj = dispocc(touch=True)
    shp = obj.show_ellipsoid(rxyz=[10.0, 15.0, 20.0])

    brepgprop_VolumeProperties(shp, prop)
    mas = prop.Mass()
    print(mas)
    obj.display.DisplayMessage(gp_Pnt(), "{:.3f}".format(mas), 30.0)
    obj.show_axs_pln(scale=10.0)
    obj.ShowOCC()
