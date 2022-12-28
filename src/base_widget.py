import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os
import time
import shutil
import subprocess
import csv
import argparse
from datetime import date, datetime
from linecache import getline, clearcache, updatecache

basename = os.path.dirname(__file__) + "/"

from PyQt5 import Qt, QtCore, QtGui, QtWidgets, sip
from PyQt5.QtWidgets import qApp, QWidget, QApplication
from PyQt5.QtWidgets import QMenu, QAction, QFileDialog
from PyQt5.QtWidgets import QFormLayout, QHBoxLayout, QVBoxLayout, QGridLayout, QBoxLayout
from PyQt5.QtWidgets import QTableWidget, QTableWidgetItem, QTableView
from PyQt5.QtWidgets import QListWidget, QListWidgetItem
from PyQt5.QtWidgets import QComboBox, QLabel, QLineEdit, QGroupBox, QPushButton
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

sys.path.append(os.path.join("./"))
from base_qtApp import MainWindow
from base import plot2d, create_tempnum

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)


class UVWidget(QtWidgets.QWidget):

    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent)
        self.left = 200
        self.top = 200
        self.width = 125
        self.height = 125
        self.setGeometry(self.left, self.top, self.width, self.height)

        self.setting = QtCore.QSettings(basename + "temp_setting/UVWidget.ini",
                                        QtCore.QSettings.IniFormat)
        self.setting.setFallbacksEnabled(False)
        self.move(self.setting.value("pos", self.pos()))
        self.resize(self.setting.value("size", self.size()))
        font = self.font()
        font.setPointSize(self.setting.value("font", 9, int))
        self.setFont(font)

        self.closeButton = QtWidgets.QPushButton("Set & Close", self)
        self.uval_text = QtWidgets.QLineEdit("", self)
        self.vval_text = QtWidgets.QLineEdit("", self)
        self.uval_labl = QtWidgets.QLabel("", self)
        self.vval_labl = QtWidgets.QLabel("", self)
        self.uval_text.setText(self.setting.value("u", "0.0", str))
        self.vval_text.setText(self.setting.value("v", "0.0", str))

        self.layout = QtWidgets.QGridLayout(self)
        self.layout.addWidget(self.uval_labl, 1, 0)
        self.layout.addWidget(self.vval_labl, 2, 0)
        self.layout.addWidget(self.uval_text, 1, 1)
        self.layout.addWidget(self.vval_text, 2, 1)
        self.layout.addWidget(self.closeButton, 3, 0)
        self.setLayout(self.layout)

    def close_self(self):
        if self.setting != None:
            # Write window size and position to config file
            self.setting.setValue("size", self.size())
            self.setting.setValue("pos", self.pos())
            self.setting.setValue("font", self.font().pointSize())
            self.setting.setValue("u", self.uval_text.text())
            self.setting.setValue("v", self.vval_text.text())
        self.close()


if __name__ == '__main__':
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument("--pxyz", dest="pxyz",
                        default=[0.0, 0.0, 0.0], type=float, nargs=3)
    opt = parser.parse_args()
    print(opt)

    print("Qt version:", QtCore.QT_VERSION_STR)
    print("PyQt version:", QtCore.PYQT_VERSION_STR)
    print("SIP version:", sip.SIP_VERSION_STR)
