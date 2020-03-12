import numpy as np
import matplotlib.pyplot as plt
import sys
import time
import os

from base import plotocc

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pnt
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCC.Core.AIS import AIS_Manipulator, AIS_Dimension

obj = plotocc()
obj.SaveMenu()
obj.AddManipulator()
obj.display.View.TriedronErase()

box = BRepPrimAPI_MakeBox(gp_Pnt(0, 0, 10), 10., 20., 30.).Shape()
obj.display.DisplayShape(box, update=True)
#manip = AIS_Manipulator()
# manip.Attach(ais)

obj.show_axs_pln(scale=10)
obj.show()
