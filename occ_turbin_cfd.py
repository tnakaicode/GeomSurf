import numpy as np
import matplotlib.pyplot as plt
from math import atan2, sin, radians
from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire
from OCC.Core.gp import gp_Pnt, gp_Trsf, gp_Ax1, gp_Dir
from OCC.Display.SimpleGui import init_display
from OCCUtils.Topology import Topo

def naca4_profile(m, p, t, c):
    x = np.linspace(0, c, 100)
    yt = 5 * t * (0.2969 * np.sqrt(x / c) - 0.1260 * (x / c) - 0.3516 * (x / c)**2 + 0.2843 * (x / c)**3 - 0.1036 * (x / c)**4)
    yc = np.where(x <= p * c, m * (x / p ** 2) * (2 * p - x / c), m * ((c - x) / (1 - p) ** 2) * (1 + x / c - 2 * p))
    dyc_dx = np.where(x <= p * c, 2 * m / p ** 2 * (p - x / c), 2 * m / (1 - p) ** 2 * (p - x / c))
    theta = np.arctan(dyc_dx)
    xu = x - yt * np.sin(theta)
    yu = yc + yt * np.cos(theta)
    xl = x + yt * np.sin(theta)
    yl = yc - yt * np.cos(theta)
    return xu, yu, xl, yl

def create_twisted_sections(num_sections=30, height=5.0, twist_angle=30.0):
    m = 0.02
    p = 0.4
    t = 0.12
    c = 1.0

    xu, yu, xl, yl = naca4_profile(m, p, t, c)
    
    sections = []
    for i in range(num_sections):
        z = (height / num_sections) * i
        angle = twist_angle * (i / num_sections)
        
        edges = []
        trsf = gp_Trsf()
        trsf.SetRotation(gp_Ax1(gp_Pnt(0, 0, z), gp_Dir(0, 0, 1)), radians(angle))
        
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

def calc_pressure_distribution(xu, yu, xl, yl, V_inf):
    num_points = len(xu)
    Cp_upper = np.zeros(num_points)
    Cp_lower = np.zeros(num_points)

    for i in range(1, num_points - 1):
        dy = yu[i + 1] - yu[i]
        dx = xu[i + 1] - xu[i]
        theta = atan2(dy, dx)
        V = V_inf * sin(theta)
        Cp_upper[i] = 1 - (V / V_inf)**2

        dy = yl[i + 1] - yl[i]
        dx = xl[i + 1] - xl[i]
        theta = atan2(dy, dx)
        V = V_inf * sin(theta)
        Cp_lower[i] = 1 - (V / V_inf)**2
    
    # 端部のデータ点を除外
    Cp_upper = Cp_upper[1:-1]
    Cp_lower = Cp_lower[1:-1]
    xu = xu[1:-1]
    xl = xl[1:-1]

    return Cp_upper, Cp_lower, xu, xl

def plot_pressure_distribution(x, Cp_upper, Cp_lower):
    plt.figure(figsize=(10, 5))
    plt.plot(x, Cp_upper, label='Upper Surface Cp', marker='o')
    plt.plot(x, Cp_lower, label='Lower Surface Cp', marker='x')
    plt.gca().invert_yaxis()
    plt.title('Pressure Distribution on Airfoil Surface')
    plt.xlabel('Chord Position')
    plt.ylabel('Pressure Coefficient (Cp)')
    plt.legend()
    plt.show()

def display_blade_and_pressure():
    sections = create_twisted_sections()

    V_inf = 1  # 流速 m/s

    # 複数断面での圧力分布を計算
    for section in sections:
        vertices = [BRep_Tool.Pnt(v) for v in Topo(section).vertices()]
        mid_idx = len(vertices) // 2
        xu_vertices = vertices[:mid_idx]
        xl_vertices = vertices[mid_idx:]

        xu = np.array([p.X() for p in xu_vertices])
        yu = np.array([p.Y() for p in xu_vertices])
        xl = np.array([p.X() for p in xl_vertices])
        yl = np.array([p.Y() for p in xl_vertices])
        
        Cp_upper, Cp_lower, xu, xl = calc_pressure_distribution(xu, yu, xl, yl, V_inf)
        x = np.linspace(0, 1.0, len(xu))
        plot_pressure_distribution(x, Cp_upper, Cp_lower)
    
    display, start_display, add_menu, add_function_to_menu = init_display()

    for section in sections:
        display.DisplayShape(section, update=False)
    
    mid_section = sections[len(sections) // 2]
    display.DisplayShape(mid_section, update=True)
    start_display()

if __name__ == "__main__":
    display_blade_and_pressure()
