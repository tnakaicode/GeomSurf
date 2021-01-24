import numpy as np
import os
import urllib.request as urllib2


def uiuc_database(name="dae51"):
    dat_file = "./uiuc_dat/{}.dat".format(name)
    if os.path.exists(dat_file):
        fp = open(dat_file, "r")
        fp_lines = fp.readlines()
    else:
        uiuc_url = 'http://m-selig.ae.illinois.edu/ads/coord_seligFmt/'
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
    data = []
    for idx, line in enumerate(fp_lines[1:]):
        data.append([float(v) for v in line.split()])
    return np.array(data)


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
