from geomdl import CPGen
from geomdl import BSpline
from geomdl import utilities
from geomdl.visualization import VisMPL
from matplotlib import cm
import matplotlib
print(matplotlib.get_backend())

# Generate a plane with the dimensions 50x100
surfgrid = CPGen.Grid(50, 100)

# Generate a grid of 25x30
surfgrid.generate(50, 60)

# Generate bumps on the grid
surfgrid.bumps(num_bumps=5, bump_height=20, base_extent=8)

# Create a BSpline surface instance
surf = BSpline.Surface()

# Set degrees
surf.degree_u = 3
surf.degree_v = 3

# Get the control points from the generated grid
surf.ctrlpts2d = surfgrid.grid

# Set knot vectors
surf.knotvector_u = utilities.generate_knot_vector(
    surf.degree_u, surf.ctrlpts_size_u)
surf.knotvector_v = utilities.generate_knot_vector(
    surf.degree_v, surf.ctrlpts_size_v)

# Set sample size
surf.sample_size = 100

# Set visualization component
surf.vis = VisMPL.VisSurface(ctrlpts=False, legend=False)

# Plot the surface
surf.render(colormap=cm.terrain)
