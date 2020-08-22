import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import json
import shutil
import urllib.request as urllib2
from optparse import OptionParser

sys.path.append(os.path.join("../"))
from base import plot2d, plot3d

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)


def coord_database(name="dae51"):
    dat_file = "./coord_dat/{}.dat".format(name)
    if os.path.exists(dat_file):
        fp = open(dat_file, "r")
        fp_lines = fp.readlines()
    else:
        uiuc_url = 'http://m-selig.ae.illinois.edu/ads/coord/'
        foil_dat_url = uiuc_url + '{}.dat'.format(name)
        fp = urllib2.urlopen(foil_dat_url)
        fp_lines = fp.readlines()
        qp = open(dat_file, "w")
        for line in fp_lines:
            txt = ""
            dat = line.split()
            for t in dat:
                txt += str(t.decode("utf-8")) + "\t"
            qp.write(txt + "\n")

    upp_data, bot_data = [], []
    nx, ny = [int(float(v)) for v in fp_lines[1].split()]
    print(nx, ny)
    # print(fp_lines[3].split())
    for idx, line in enumerate(fp_lines[3:nx + 3]):
        upp_data.append([float(v) for v in line.split()])
    # print(fp_lines[nx+3+1].split())
    for idx, line in enumerate(fp_lines[nx + 3 + 1:nx + ny + 3 + 2]):
        bot_data.append([float(v) for v in line.split()])
    return np.array(upp_data), np.array(bot_data)


if __name__ == '__main__':
    argvs = sys.argv
    parser = OptionParser()
    parser.add_option("--dir", dest="dir", default="./")
    parser.add_option("--name", dest="name", default="dae51")
    opt, argc = parser.parse_args(argvs)
    print(opt, argc)

    obj = plot2d(aspect="equal")
    cfg = json.load(open("./cfg.json", "r"))
    dat_name = opt.name
    dat1, dat2 = coord_database(dat_name)
    print(dat_name)
    obj.axs.plot(dat1[:, 0], dat1[:, 1])
    obj.axs.plot(dat2[:, 0], dat2[:, 1])
    obj.axs.set_title(dat_name)
    obj.SavePng("./coord_dat/{}.png".format(dat_name))
