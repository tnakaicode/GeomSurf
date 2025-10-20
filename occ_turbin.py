import numpy as np
from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepBuilderAPI import (
    BRepBuilderAPI_MakeEdge,
    BRepBuilderAPI_MakeWire,
    BRepBuilderAPI_MakeFace,
)
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_ThruSections
from OCC.Core.gp import gp_Trsf, gp_Pnt, gp_Vec, gp_Ax1, gp_Dir
from OCC.Display.SimpleGui import init_display
from OCC.Core.Quantity import Quantity_Color, Quantity_TOC_RGB


def naca4_profile(m, p, t, c):
    x = np.linspace(0, c, 100)
    yt = (
        5
        * t
        * (
            0.2969 * np.sqrt(x / c)
            - 0.1260 * (x / c)
            - 0.3516 * (x / c) ** 2
            + 0.2843 * (x / c) ** 3
            - 0.1036 * (x / c) ** 4
        )
    )
    yc = np.where(
        x <= p * c,
        m * (x / p**2) * (2 * p - x / c),
        m * ((c - x) / (1 - p) ** 2) * (1 + x / c - 2 * p),
    )
    dyc_dx = np.where(
        x <= p * c, 2 * m / p**2 * (p - x / c), 2 * m / (1 - p) ** 2 * (p - x / c)
    )
    theta = np.arctan(dyc_dx)
    xu = x - yt * np.sin(theta)
    yu = yc + yt * np.cos(theta)
    xl = x + yt * np.sin(theta)
    yl = yc - yt * np.cos(theta)
    return xu, yu, xl, yl


def create_twisted_sections():
    m = 0.02
    p = 0.4
    t = 0.12
    c = 1.0
    height = 5.0
    twist_angle = 30.0

    xu, yu, xl, yl = naca4_profile(m, p, t, c)

    sections = []
    num_sections = 10  # 断面の数

    for i in range(num_sections):
        z = (height / num_sections) * i
        angle = twist_angle * (i / num_sections)

        edges = []
        trsf = gp_Trsf()
        trsf.SetRotation(gp_Ax1(gp_Pnt(0, 0, z), gp_Dir(0, 0, 1)), np.radians(angle))

        for j in range(len(xu) - 1):
            p1 = gp_Pnt(xu[j], yu[j], z)
            p2 = gp_Pnt(xu[j + 1], yu[j + 1], z)
            p1.Transform(trsf)
            p2.Transform(trsf)
            edge = BRepBuilderAPI_MakeEdge(p1, p2).Edge()
            edges.append(edge)

        for j in range(len(xl) - 1):
            p1 = gp_Pnt(xl[j], yl[j], z)
            p2 = gp_Pnt(xl[j + 1], yl[j + 1], z)
            p1.Transform(trsf)
            p2.Transform(trsf)
            edge = BRepBuilderAPI_MakeEdge(p1, p2).Edge()
            edges.append(edge)

        wire = BRepBuilderAPI_MakeWire()
        for edge in edges:
            wire.Add(edge)

        sections.append(wire.Wire())

    return sections


def create_twisted_blade():
    sections = create_twisted_sections()

    brep_offset = BRepOffsetAPI_ThruSections()
    for wire in sections:
        brep_offset.AddWire(wire)

    brep_offset.CheckCompatibility(True)  # 必要に応じて互換性チェックを行います
    blade = brep_offset.Shape()

    return blade


def display_blade():
    display, start_display, add_menu, add_function_to_menu = init_display()

    blade = create_twisted_blade()
    display.DisplayShape(blade, update=True)
    start_display()


if __name__ == "__main__":
    display_blade()
    # https://github.com/Mopolino8/cuIBM-FSI/blob/c84fde3cb218ccb962c2717fc58f78c5550ce94d/scripts/python/naca00xx.py
