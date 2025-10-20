import numpy as np
from OCC.Core.TColgp import TColgp_Array2OfPnt
from OCC.Core.gp import gp_Pnt
from OCC.Core.GeomAPI import GeomAPI_PointsToBSplineSurface
from OCC.Core.GeomAPI import GeomAPI_ProjectPointOnSurf

points = [
    [gp_Pnt(0, 0, 0), gp_Pnt(1, 0, 1), gp_Pnt(2, 0, 0)],
    [gp_Pnt(0, 1, 1), gp_Pnt(1, 1, 2), gp_Pnt(2, 1, 1)],
    [gp_Pnt(0, 2, 0), gp_Pnt(1, 2, 1), gp_Pnt(2, 2, 0)],
]

# Create TColgp_Array2OfPnt
num_u = len(points)
num_v = len(points[0])
array = TColgp_Array2OfPnt(1, num_u, 1, num_v)

for i in range(1, num_u + 1):
    for j in range(1, num_v + 1):
        array.SetValue(i, j, points[i - 1][j - 1])

# Create the B-spline surface
bspline_surface = GeomAPI_PointsToBSplineSurface(array).Surface()

# Function to project points onto the B-spline surface and get (u, v) parameters


def find_uv_on_surface(surface, point):
    projector = GeomAPI_ProjectPointOnSurf(point, surface)
    if projector.NbPoints() > 0:
        u, v = projector.LowerDistanceParameters()
        return u, v
    return None, None


# Evaluate the accuracy
differences = []
for i in range(num_u):
    for j in range(num_v):
        orig_point = points[i][j]
        u, v = find_uv_on_surface(bspline_surface, orig_point)
        if u is not None and v is not None:
            approx_point = bspline_surface.Value(u, v)
            diff = orig_point.Distance(approx_point)
            differences.append(diff)
        else:
            print(f"Failed to project point {orig_point}")

# Print or analyze the differences
mean_diff = np.mean(differences)
max_diff = np.max(differences)
print(f"Mean difference: {mean_diff}")
print(f"Max difference: {max_diff}")
