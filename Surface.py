import numpy as np
import matplotlib.pyplot as plt
import sys

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Lin, gp_Sphere
from OCC.Core.Geom import Geom_Circle, Geom_Conic, Geom_ConicalSurface
from OCC.Core.IMeshData import IMeshData_Curve

"""
#include <IMeshData_Status.hxx>
#include <IMeshTools_Parameters.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
Standard_Boolean meshing_explicit_parameters()
{
  BRepMesh_IncrementalMesh aMesher (aShape, 0.1, Standard_False, 0.5, Standard_True);
  const Standard_Integer aStatus = aMesher.GetStatusFlags();
  return !aStatus;
}
Standard_Boolean meshing_new()
{
  IMeshTools_Parameters aMeshParams;
  aMeshParams.Deflection               = 0.1;
  aMeshParams.Angle                    = 0.5;
  aMeshParams.Relative                 = Standard_False;
  aMeshParams.InParallel               = Standard_True;
  aMeshParams.MinSize                  = Precision::Confusion();
  aMeshParams.InternalVerticesMode     = Standard_True;
  aMeshParams.ControlSurfaceDeflection = Standard_True;
  BRepMesh_IncrementalMesh aMesher (aShape, aMeshParams);
  const Standard_Integer aStatus = aMesher.GetStatusFlags();
  return !aStatus;
}
"""