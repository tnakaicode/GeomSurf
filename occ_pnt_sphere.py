import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time
import argparse
from linecache import getline, clearcache, updatecache

basename = os.path.dirname(__file__)

sys.path.append(os.path.join("./"))
from src.base_occ import dispocc

import logging

logging.getLogger("matplotlib").setLevel(logging.ERROR)

from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.gp import gp_Ax1, gp_Ax2, gp_Ax3
from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Compound
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Extend.TopologyUtils import TopologyExplorer
from OCCUtils.Common import vertex2pnt
from OCCUtils.Construct import make_box, make_line, make_wire, make_edge
from OCCUtils.Construct import make_plane, make_polygon
from OCCUtils.Construct import point_to_vector, vector_to_point
from OCCUtils.Construct import dir_to_vec, vec_to_dir


def r8mat_uniform_ab(m, n, a, b, seed):
    # *****************************************************************************80
    #
    # R8MAT_UNIFORM_AB returns a scaled pseudorandom R8MAT.
    #
    #  Discussion:
    #
    #    An R8MAT is an array of R8's.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    08 April 2013
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Reference:
    #
    #    Paul Bratley, Bennett Fox, Linus Schrage,
    #    A Guide to Simulation,
    #    Second Edition,
    #    Springer, 1987,
    #    ISBN: 0387964673,
    #    LC: QA76.9.C65.B73.
    #
    #    Bennett Fox,
    #    Algorithm 647:
    #    Implementation and Relative Efficiency of Quasirandom
    #    Sequence Generators,
    #    ACM Transactions on Mathematical Software,
    #    Volume 12, Number 4, December 1986, pages 362-376.
    #
    #    Pierre L'Ecuyer,
    #    Random Number Generation,
    #    in Handbook of Simulation,
    #    edited by Jerry Banks,
    #    Wiley, 1998,
    #    ISBN: 0471134031,
    #    LC: T57.62.H37.
    #
    #    Peter Lewis, Allen Goodman, James Miller,
    #    A Pseudo-Random Number Generator for the System/360,
    #    IBM Systems Journal,
    #    Volume 8, Number 2, 1969, pages 136-143.
    #
    #  Parameters:
    #
    #    Input, integer M, N, the number of rows and columns in the array.
    #
    #    Input, real A, B, the range of the pseudorandom values.
    #
    #    Input, integer SEED, the integer "seed" used to generate
    #    the output random number.
    #
    #    Output, real R(M,N), an array of random values between 0 and 1.
    #
    #    Output, integer SEED, the updated seed.  This would
    #    normally be used as the input seed on the next call.
    #

    i4_huge = 2147483647

    seed = int(seed)

    if seed < 0:
        seed = seed + i4_huge

    if seed == 0:
        print("")
        print("R8MAT_UNIFORM_AB - Fatal error!")
        print("  Input SEED = 0!")
        exit("R8MAT_UNIFORM_AB - Fatal error!")

    r = np.zeros([m, n])

    for j in range(0, n):
        for i in range(0, m):
            k = seed // 127773

            seed = 16807 * (seed - k * 127773) - k * 2836
            seed = seed % i4_huge

            if seed < 0:
                seed = seed + i4_huge

            r[i, j] = a + (b - a) * seed * 4.656612875e-10

    if m == 1 and n == 1:
        r = r[0, 0]

    return r, seed


def polygon_sample(nv, v, n, seed):
    # *****************************************************************************80
    #
    # POLYGON_SAMPLE uniformly samples a polygon.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    18 October 2015
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Parameters:
    #
    #    Input, integer NV, the number of vertices.
    #
    #    Input, real V(NV,2), the vertices of the polygon, listed in
    #    counterclockwise order.
    #
    #    Input, integer N, the number of points to create.
    #
    #    Input/output, integer SEED, a seed for the random
    #    number generator.
    #
    #    Output, real S(2,N), the points.
    #

    #  Triangulate the polygon.
    x = np.zeros(nv)
    y = np.zeros(nv)
    for j in range(0, nv):
        x[j] = v[j, 0]
        y[j] = v[j, 1]

    triangles = polygon_triangulate(nv, x, y)

    #  Determine the areas of each triangle.
    area_triangle = np.zeros(nv - 2)

    area_polygon = 0.0
    for i in range(0, nv - 2):
        area_triangle[i] = triangle_area(
            v[triangles[i, 0], 0],
            v[triangles[i, 0], 1],
            v[triangles[i, 1], 0],
            v[triangles[i, 1], 1],
            v[triangles[i, 2], 0],
            v[triangles[i, 2], 1],
        )
        area_polygon = area_polygon + area_triangle[i]

    #  Normalize the areas.
    area_relative = np.zeros(nv - 1)

    for i in range(0, nv - 2):
        area_relative[i] = area_triangle[i] / area_polygon

    #  Replace each area by the sum of itself and all previous ones.
    area_cumulative = np.zeros(nv - 2)

    area_cumulative[0] = area_relative[0]
    for i in range(1, nv - 2):
        area_cumulative[i] = area_relative[i] + area_cumulative[i - 1]

    s = np.zeros([n, 2])
    for j in range(0, n):
        #  Choose triangle I at random, based on areas.
        area_percent, seed = r8mat_uniform_ab(1, 1, seed)

        for k in range(0, nv - 2):
            i = k
            if area_percent <= area_cumulative[k]:
                break

        #  Now choose a point at random in triangle I.
        r, seed = r8mat_uniform_ab(1, 2, 0, 1, seed)

        if 1.0 < r[0] + r[1]:
            r[0] = 1.0 - r[0]
            r[1] = 1.0 - r[1]

        s[j, 0] = (
            (1.0 - r[0] - r[1]) * v[triangles[i, 0], 0]
            + r[0] * v[triangles[i, 1], 0]
            + r[1] * v[triangles[i, 2], 0]
        )

        s[j, 1] = (
            (1.0 - r[0] - r[1]) * v[triangles[i, 0], 1]
            + r[0] * v[triangles[i, 1], 1]
            + r[1] * v[triangles[i, 2], 1]
        )

    return s, seed


if __name__ == "__main__":
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument(
        "--pxyz", dest="pxyz", default=[0.0, 0.0, 0.0], type=float, nargs=3
    )
    opt = parser.parse_args()
    print(opt)

    obj = dispocc(touch=True)
    axs = gp_Ax3()

    n = 100
    xyz = np.random.normal(0.0, 1.0, [n, 3])
    for p in xyz:
        obj.display.DisplayShape(gp_Pnt(*p))

    obj.show_axs_pln(scale=1)
    obj.ShowOCC()
