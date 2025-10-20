# ------------------------------------------------------------------------------
#
#  Gmsh Python tutorial 18
#
#  Periodic meshes
#
# ------------------------------------------------------------------------------

# Periodic meshing constraints can be imposed on surfaces and curves.

import gmsh
import math
import os
import sys
import ctypes

from OCCUtils.Construct import make_box

gmsh.initialize()

gmsh.model.add("occ_t18")

# Let's use the OpenCASCADE geometry kernel to build two geometries.

# The first geometry is very simple: a unit cube with a non-uniform mesh size
# constraint (set on purpose to be able to verify visually that the periodicity
# constraint works!):

box = make_box(1,1,1)

gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1, 1)
#gmsh.model.occ.importShapesNativePointer(ctypes.py_object(box))
#gmsh.model.occ.importShapesNativePointer(box.__class__)
gmsh.model.occ.synchronize()


gmsh.model.mesh.setSize(gmsh.model.getEntities(0), 0.1)
gmsh.model.mesh.setSize([(0, 1)], 0.02)

# To impose that the mesh on surface 2 (the right side of the cube) should
# match the mesh from surface 1 (the left side), the following periodicity
# constraint is set:
translation = [1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
gmsh.model.mesh.setPeriodic(2, [2], [1], translation)

# The periodicity transform is provided as a 4x4 affine transformation matrix,
# given by row.

# During mesh generation, the mesh on surface 2 will be created by copying
# the mesh from surface 1.

# Multiple periodicities can be imposed in the same way:
gmsh.model.mesh.setPeriodic(2, [6], [5],
                            [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1])
gmsh.model.mesh.setPeriodic(2, [4], [3],
                            [1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1])

# For more complicated cases, finding the corresponding surfaces by hand can
# be tedious, especially when geometries are created through solid
# modelling. Let's construct a slightly more complicated geometry.

# Ask OpenCASCADE to compute more accurate bounding boxes of entities using
# the STL mesh:
gmsh.option.setNumber("Geometry.OCCBoundsUseStl", 1)

# We then retrieve all the volumes in the bounding box of the original cube,
# and delete all the parts outside it:
eps = 1e-3
vin = gmsh.model.getEntitiesInBoundingBox(2 - eps, -eps, -eps, 2 + 1 + eps,
                                          1 + eps, 1 + eps, 3)

# We now set a non-uniform mesh size constraint (again to check results
# visually):
p = gmsh.model.getBoundary(vin, False, False, True)  # Get all points
gmsh.model.mesh.setSize(p, 0.1)
p = gmsh.model.getEntitiesInBoundingBox(2 - eps, -eps, -eps, 2 + eps, eps, eps,
                                        0)
gmsh.model.mesh.setSize(p, 0.001)

# We now identify corresponding surfaces on the left and right sides of the
# geometry automatically.

# First we get all surfaces on the left:
sxmin = gmsh.model.getEntitiesInBoundingBox(2 - eps, -eps, -eps, 2 + eps,
                                            1 + eps, 1 + eps, 2)

for i in sxmin:
    # Then we get the bounding box of each left surface
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(i[0], i[1])
    # We translate the bounding box to the right and look for surfaces inside
    # it:
    sxmax = gmsh.model.getEntitiesInBoundingBox(xmin - eps + 1, ymin - eps,
                                                zmin - eps, xmax + eps + 1,
                                                ymax + eps, zmax + eps, 2)
    # For all the matches, we compare the corresponding bounding boxes...
    for j in sxmax:
        xmin2, ymin2, zmin2, xmax2, ymax2, zmax2 = gmsh.model.getBoundingBox(
            j[0], j[1])
        xmin2 -= 1
        xmax2 -= 1
        # ...and if they match, we apply the periodicity constraint
        if (abs(xmin2 - xmin) < eps and abs(xmax2 - xmax) < eps and
            abs(ymin2 - ymin) < eps and abs(ymax2 - ymax) < eps and
            abs(zmin2 - zmin) < eps and abs(zmax2 - zmax) < eps):
            gmsh.model.mesh.setPeriodic(2, [j[1]], [i[1]], translation)

gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.model.mesh.generate(3)
gmsh.write("occ_t18.msh")

#gmsh.fltk.run()
#gmsh.finalize()
