import numpy as np
import time
from scipy.linalg import lu, svd

from OCC.Core.math import math_Matrix, math_Vector
from OCC.Core.math import math_MultipleVarFunction
from OCC.Core.math import math_GaussLeastSquare, math_SVD
from OCC.Core.math import math_PSO

# import MultipleVarFunction: GCPnts_DistFunction
from OCC.Core.GCPnts import GCPnts_DistFunction, GCPnts_DistFunctionMV
# GCPnts_DistFunction not wrappeed

# import MultipleVarFunction: Extrema_GlobOptFuncCQuadric
from OCC.Core.Extrema import Extrema_GlobOptFuncCQuadric
from OCC.Core.Adaptor3d import Adaptor3d_Curve, Adaptor3d_Surface
from OCC.Core.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Surface
from OCC.Core.GeomAdaptor import GeomAdaptor_Curve, GeomAdaptor_Surface
