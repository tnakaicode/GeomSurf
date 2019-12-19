import numpy as np
import matplotlib.pyplot as plt
import sys

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Lin, gp_Sphere
from OCC.Core.BRepMesh import BRepMesh_DelaunayBaseMeshAlgo
from OCC.Extend.DataExchange import read_step_file, write_step_file
from OCCUtils.Construct import make_n_sided, make_n_sections
from OCCUtils.Construct import make_edge, make_polygon
from OCCUtils.Construct import vec_to_dir, dir_to_vec
from OCCUtils.Topology import Topo

from base import plotocc, gen_ellipsoid

# BRepMesh_DelaunayDeflectionControlMeshAlgo
# https://dev.opencascade.org/doc/refman/html/class_b_rep_mesh___delaunay_base_mesh_algo.html

if __name__ == '__main__':
    print("ok")
