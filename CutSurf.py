import numpy as np
import matplotlib.pyplot as plt
import sys

from OCC.Display.SimpleGui import init_display
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.gp import gp_Lin, gp_Sphere
from OCC.Core.gp import gp_Mat, gp_XYZ
from OCC.Core.gp import gp_Trsf, gp_GTrsf
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.Extend.DataExchange import write_step_file
from OCCUtils.Construct import make_edge

from base import plotocc, spl_face

# https://www.opencascade.com/doc/occt-7.4.0/refman/html/class_b_rep_builder_a_p_i___make_face.html
# Adds the wire <W> in the face <F> A general method to create a face is to give.
#
# a surface S as the support (the geometric domain) of the face, and a wire W to bound it.
# The bounds of the face can also be defined by four parameter values umin, umax, vmin, vmax
# which determine isoparametric limitations on the parametric space of the surface.
#
# In this way, a patch is defined.
# The parameter values are optional.
# If they are omitted, the natural bounds of the surface are used.
# A wire is automatically built using the defined bounds.
# Up to four edges and four vertices are created with this wire (no edge is created when the corresponding parameter value is infinite).
# Wires can then be added using the function Add to define other restrictions on the face.
# These restrictions represent holes.
# More than one wire may be added by this way,
# provided that the wires do not cross each other and that they define only one area on the surface.
# (Be careful, however, as this is not checked).
# Forbidden addition of wires Note that in this schema,
# the third case is valid if edges of the wire W are declared internal to the face.
# As a result, these edges are no longer bounds of the face.
# A default tolerance (Precision::Confusion()) is given to the face,
# this tolerance may be increased during construction of the face using various algorithms.
# Rules applied to the arguments For the surface:
#
# The surface must not be a 'null handle'.
# If the surface is a trimmed surface, the basis surface is used.
#
# For the wire: the wire is composed of connected edges, each edge having a parametric curve description in the parametric domain of the surface;
# in other words, as a pcurve. For the parameters:
# The parameter values must be in the parametric range of the surface (or the basis surface, if the surface is trimmed).
# If this condition is not satisfied, the face is not built,
# and the Error function will return BRepBuilderAPI_ParametersOutOfRange.
# The bounding parameters p1 and p2 are adjusted on a periodic surface in a given parametric direction by adding or subtracting the period to obtain p1 in the parametric range of the surface and such p2, that p2 - p1 <= Period, where Period is the period of the surface in this parametric direction.
# A parameter value may be infinite. There will be no edge and no vertex in the corresponding direction.

if __name__ == '__main__':
    obj = plotocc()
    obj.show_axs_pln(scale=1)
    obj.display.DisplayShape(obj.make_PolyWire())
    obj.display.DisplayShape(obj.make_PolyWire(num=10))

    px = np.linspace(-1, 1, 10) * 2
    py = np.linspace(-1, 1, 10) * 2
    mesh = np.meshgrid(px, py)
    surf = mesh[0]**2 / 10 + mesh[1]**2 / 20
    face = spl_face(*mesh, surf)
    poly = obj.make_EllipWire()
    api = BRepBuilderAPI_MakeFace(face, poly)
    bound_face = api.Face()
    obj.display.DisplayShape(bound_face)
    obj.display.DisplayShape(poly)

    obj.show()
