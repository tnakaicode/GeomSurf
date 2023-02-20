import numpy as np
import time
from scipy.linalg import lu, svd

from OCC.Core.math import math_Matrix, math_Vector
from OCC.Core.math import math_MultipleVarFunction
from OCC.Core.math import math_GaussLeastSquare, math_SVD
from OCC.Core.math import math_PSO, math_NewtonMinimum, math_Function

# import MultipleVarFunction: GCPnts_DistFunction
from OCC.Core.GCPnts import GCPnts_DistFunction, GCPnts_DistFunctionMV
# GCPnts_DistFunction not wrappeed

# import MultipleVarFunction: Extrema_GlobOptFuncCQuadric
from OCC.Core.Extrema import Extrema_GlobOptFuncCQuadric, Extrema_GlobOptFuncCS
from OCC.Core.Adaptor3d import Adaptor3d_Curve, Adaptor3d_Surface
from OCC.Core.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Surface
from OCC.Core.GeomAdaptor import GeomAdaptor_Curve, GeomAdaptor_Surface

from OCC.Core.gp import gp_Pnt, gp_Vec
from OCC.Extend.ShapeFactory import make_edge
from OCCUtils.Construct import make_plane

edge = make_edge(gp_Pnt(-100, -100, 10), gp_Pnt(100, 100, 10))
face = make_plane(gp_Pnt(0, 0, 0),
                  gp_Vec(0, 0, 1),
                  -100, 100, -100, 100)
edge_adap = BRepAdaptor_Curve(edge)
face_adap = BRepAdaptor_Surface(face)
func = Extrema_GlobOptFuncCQuadric(edge_adap, face_adap)
func = Extrema_GlobOptFuncCS(edge_adap, face_adap)
n = func.NbVariables()
n = 201
m = 204
math_vecl = math_Vector(1, n)
math_vecu = math_Vector(1, n)
math_vecs = math_Vector(1, m)
for i in range(n):
    math_vecl.SetValue(i + 1, 0)
    math_vecu.SetValue(i + 1, i)
for i in range(m):
    math_vecs.SetValue(i + 1, i/2)
#func.QuadricParameters(math_vecl, math_vecu)
# vec_txt = ""
# for i in range(math_vecl.Lower(), math_vecl.Upper() + 1):
#    vec_txt += f"{math_vecl.Value(i)}\t"
# print(vec_txt)
print(func.Value(math_vecs))
math_PSO(func, math_vecl, math_vecu, math_vecs, 1, 100)
# RuntimeError: Standard_DimensionErrormath_Vector::Initialized() - 
#   input vector has wrong dimensions raised from method math_PSO of class math_PSO
