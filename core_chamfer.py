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
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCC.Core.BRepFilletAPI import BRepFilletAPI_MakeChamfer, BRepFilletAPI_MakeFillet
from OCC.Core.BRepMesh import BRepMesh_IncrementalMesh
from OCC.Core.BRepMeshData import BRepMeshData_Face
from OCC.Core.ChFi3d import ChFi3d_Rational
from OCC.Core.TColgp import TColgp_Array1OfPnt2d
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_EDGE
from OCC.Core.StlAPI import StlAPI_Reader, StlAPI_Writer
from OCC.Core.IGESControl import IGESControl_Reader, IGESControl_Writer
from OCC.Extend.DataExchange import read_step_file, write_step_file, write_stl_file
from OCCUtils.Construct import make_box
from OCCUtils.Construct import make_line, make_wire, make_edge

from src.base_occ import dispocc


def write_stl_file_mesh1(a_shape, filename, mode="ascii", linear_deflection=0.9, angular_deflection=0.5):
    """ 
    export the shape to a STL file
    Be careful, the shape first need to be explicitely meshed using BRepMesh_IncrementalMesh

    a_shape: the topods_shape to export
    filename: the filename
    mode: optional, "ascii" by default. Can either be "binary"
    linear_deflection: optional, default to 0.001. Lower, more occurate mesh
    angular_deflection: optional, default to 0.5. Lower, more accurate_mesh
    """
    if a_shape.IsNull():
        raise AssertionError("Shape is null.")
    if mode not in ["ascii", "binary"]:
        raise AssertionError("mode should be either ascii or binary")
    if os.path.isfile(filename):
        print("Warning: %s file already exists and will be replaced" % filename)
    # first mesh the shape
    mesh = BRepMesh_IncrementalMesh(
        a_shape, linear_deflection, True, angular_deflection, True)
    mesh.SetParallelDefault(True)
    mesh.Perform()
    if not mesh.IsDone():
        raise AssertionError("Mesh is not done.")

    stl_exporter = StlAPI_Writer()
    if mode == "ascii":
        stl_exporter.SetASCIIMode(True)
    else:  # binary, just set the ASCII flag to False
        stl_exporter.SetASCIIMode(False)
    stl_exporter.Write(a_shape, filename)

    if not os.path.isfile(filename):
        raise IOError("File not written to disk.")


def write_stl_file_mesh2(a_shape, filename, mode="ascii", linear_deflection=0.9, angular_deflection=0.5):
    """ 
    export the shape to a STL file
    Be careful, the shape first need to be explicitely meshed using BRepMesh_IncrementalMesh

    a_shape: the topods_shape to export
    filename: the filename
    mode: optional, "ascii" by default. Can either be "binary"
    linear_deflection: optional, default to 0.001. Lower, more occurate mesh
    angular_deflection: optional, default to 0.5. Lower, more accurate_mesh
    """
    if a_shape.IsNull():
        raise AssertionError("Shape is null.")
    if mode not in ["ascii", "binary"]:
        raise AssertionError("mode should be either ascii or binary")
    if os.path.isfile(filename):
        print("Warning: %s file already exists and will be replaced" % filename)
    # first mesh the shape
    mesh = BRepMesh_IncrementalMesh(
        a_shape, linear_deflection, False, angular_deflection, True)
    mesh.SetParallelDefault(True)
    mesh.Perform()
    if not mesh.IsDone():
        raise AssertionError("Mesh is not done.")

    stl_exporter = StlAPI_Writer()
    if mode == "ascii":
        stl_exporter.SetASCIIMode(True)
    else:  # binary, just set the ASCII flag to False
        stl_exporter.SetASCIIMode(False)
    stl_exporter.Write(a_shape, filename)

    if not os.path.isfile(filename):
        raise IOError("File not written to disk.")


if __name__ == '__main__':
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    opt = parser.parse_args()
    print(opt, argvs)

    obj = dispocc()
    #
    # https://www.opencascade.com/doc/occt-7.4.0/overview/html/occt_user_guides__modeling_algos.html#occt_modalg_6
    # https://www.opencascade.com/doc/occt-7.5.0/overview/html/occt_user_guides__modeling_algos.html#occt_modalg_6
    #

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
    write_stl_file(fil.Shape(), obj.tmpdir + "box_fillet.stl")
    write_stl_file_mesh1(fil.Shape(), obj.tmpdir + "box_fillet_mesh1.stl",
                         linear_deflection=0.1E-1, angular_deflection=0.1E-1)
    write_stl_file_mesh2(fil.Shape(), obj.tmpdir + "box_fillet_mesh2.stl",
                         linear_deflection=0.1E-1, angular_deflection=0.1E-1)

    obj.display.DisplayShape(fil.Shape())
    obj.display.DisplayShape(axs.Location())
    obj.ShowOCC()
