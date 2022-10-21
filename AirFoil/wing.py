import gmsh
import os
import sys

gmsh.initialize()
path = os.path.dirname(os.path.abspath(__file__))
gmsh.merge(os.path.join(path, 'wing001.step'))

gmsh.model.mesh.generate(3)
gmsh.model.setTag(2.0, 0, 8)
gmsh.write("wing001_2.msh")

# Launch the GUI to see the results:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()
gmsh.finalize()
