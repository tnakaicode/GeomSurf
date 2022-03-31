import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import json
import shutil
import urllib.request as urllib2
import argparse

sys.path.append(os.path.join("../"))
from src.base import plot2d, plot3d
from AirFoil.selig import coord_database, uiuc_database

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

if __name__ == '__main__':
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument("--name", dest="name", default="dae51")
    parser.add_argument("--data", dest="data", default="coord")
    opt = parser.parse_args()
    print(opt, argvs)

    obj = plot2d(aspect="equal")
    cfg = json.load(open("./cfg.json", "r"))
    dat_name = opt.name
    if opt.data == "coord":
        dat1, dat2 = coord_database(dat_name)
    elif opt.data == "uiuc":
        dat1, dat2 = uiuc_database(dat_name)
    else:
        dat1, dat2 = coord_database(dat_name)
    print(dat_name, dat1.shape)
    obj.axs.plot(dat1[:, 0], dat1[:, 1])
    obj.axs.plot(dat2[:, 0], dat2[:, 1])
    obj.axs.set_title(dat_name)
    obj.SavePng("./coord_dat/{}.png".format(dat_name))
