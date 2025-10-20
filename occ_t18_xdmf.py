import meshio

msh = meshio.read("occ_t18.msh")
for cell in msh.cells:
    if cell.type == "triangle":
        triangle_cells = cell.data
    elif cell.type == "tetra":
        tetra_cells = cell.data

for key in msh.cell_data_dict["gmsh:physical"].keys():
    if key == "triangle":
        triangle_data = msh.cell_data_dict["gmsh:physical"][key]
    elif key == "tetra":
        tetra_data = msh.cell_data_dict["gmsh:physical"][key]

tetra_mesh = meshio.Mesh(points=msh.points, cells={"tetra": tetra_cells})
meshio.write("occ_t18_mesh.xdmf", tetra_mesh)

triangle_mesh = meshio.Mesh(points=msh.points,
                            cells=[("triangle", triangle_cells)],
                            cell_data={"name_to_read": [triangle_data]})
meshio.write("occ_t18_mf.xdmf", triangle_mesh)
