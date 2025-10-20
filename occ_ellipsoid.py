import sys
from OCC.Core.gp import (
    gp_Pnt,
    gp_Ax2,
    gp_Dir,
    gp_Pln,
    gp_Ax3,
    gp_Mat,
    gp_XYZ,
    gp_GTrsf,
    gp_Trsf,
)
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeSphere
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_GTransform
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Cut, BRepAlgoAPI_Section
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.Display.SimpleGui import init_display
from OCC.Core.TopoDS import TopoDS_Compound
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_EDGE
from OCC.Core.BRep import BRep_Builder


def create_ellipsoid(
    radius=1, a1=1, a2=1, a3=1, center=gp_Pnt(0, 0, 0), direction=gp_Dir(0, 0, 1)
):
    """Create an ellipsoid by scaling a sphere."""
    axs = gp_Ax3(center, direction)
    sphere = BRepPrimAPI_MakeSphere(axs.Ax2(), radius).Shape()
    scale_transform = gp_GTrsf(
        gp_Mat(a1, 0, 0, 0, a2, 0, 0, 0, a3), gp_XYZ(0.0, 0.0, 0.0)
    )
    trf = gp_Trsf()
    trf.SetTransformation(axs)
    ellipsoid = BRepBuilderAPI_GTransform(sphere, scale_transform).Shape()
    return ellipsoid


if __name__ == "__main__":
    # Initialize the display
    display, start_display, add_menu, add_function_to_menu = init_display()

    # Create a large ellipsoid
    radius_large = 1
    a1_large, a2_large, a3_large = 50, 30, 20
    center_large = gp_Pnt(0, 0, 0)
    direction_large = gp_Dir(0, 0, 1)
    ellipsoid_large = create_ellipsoid(
        radius_large, a1_large, a2_large, a3_large, center_large, direction_large
    )

    # Create a small ellipsoid
    radius_small = 0.5
    a1_small, a2_small, a3_small = 25, 15, 10
    center_small = gp_Pnt(0, 0, 0)
    direction_small = gp_Dir(1, 1, 1)
    ellipsoid_small = create_ellipsoid(
        radius_small, a1_small, a2_small, a3_small, center_small, direction_small
    )

    # Subtract the small ellipsoid from the large ellipsoid to create a hollow shape
    hollow_shape = BRepAlgoAPI_Cut(ellipsoid_large, ellipsoid_small).Shape()

    # Create a plane for the cross-section
    plane = gp_Pln(gp_Ax3(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0)))

    # Compute the cross-section
    section = BRepAlgoAPI_Section(hollow_shape, plane).Shape()

    # Create a compound to hold the edges
    compound = TopoDS_Compound()
    builder = BRep_Builder()
    builder.MakeCompound(compound)

    # Extract edges from the section and add them to the compound
    explorer = TopExp_Explorer(section, TopAbs_EDGE)
    while explorer.More():
        edge = explorer.Current()
        builder.Add(compound, edge)
        explorer.Next()
    display.DisplayShape(compound)

    # Create a face from the compound
    # face = BRepBuilderAPI_MakeFace(compound).Face()

    # Display the hollow shape with transparency
    display.DisplayShape(hollow_shape, transparency=0.9)

    # Display the cross-section face
    # display.DisplayShape(face, update=True)

    # Start the display loop
    display.FitAll()
    start_display()
