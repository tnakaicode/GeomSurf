import numpy as np
import matplotlib.pyplot as plt
from linecache import getline, clearcache


def curvature(px, r, s):
    """( x + sx )**2 / 2*rx + ( y + sy )**2 / 2*ry"""
    if (r == 0):
        py = np.zeros_like(px + s)
    else:
        py = (px + s)**2 / (2 * r)
    return py


def fit_surf_2d_xyz(mesh, xyz, indx=5, idx_file="surf_idx.txt"):
    X = xyz[:, 0]
    Y = xyz[:, 1]
    A = [np.ones_like(X)]
    for i in range(1, indx):
        for j in range(i + 1):
            ix, iy = j, i - j
            A.append(X**(ix) * Y**(iy))

    A = np.array(A).T
    B = xyz[:, 2]
    coeff, r, rank, s = np.linalg.lstsq(A, B)
    print(coeff)
    print(r)
    print(s)

    X = mesh[0].flatten()
    Y = mesh[1].flatten()
    A = [np.ones_like(X)]
    for i in range(1, indx):
        for j in range(i + 1):
            ix, iy = j, i - j
            A.append(X**(ix) * Y**(iy))
    A = np.array(A).T
    surf = np.zeros_like(X)
    for i, c in enumerate(coeff):
        surf += A[:, i] * c

    xs, xe = mesh[0][0, 0], mesh[0][0, -1]
    ys, ye = mesh[1][0, 0], mesh[1][-1, 0]

    fp = open(idx_file, "w")
    fp.write('{:2d}\n'.format(indx - 1))
    fp.write('x\t{:.2f}\t{:.2f}\n'.format(xs, xe))
    fp.write('y\t{:.2f}\t{:.2f}\n'.format(ys, ye))
    fp.write('\n')
    fp.write('{}\t{}\t{}\n'.format("ix", "iy", "coeff"))
    for i in range(indx):
        for j in range(i + 1):
            ix, iy = j, i - j
            fp.write(
                '{:2d}\t{:2d}\t{:0.10E}\n'.format(
                    ix, iy, coeff[sum(range(i + 1)) + j]
                )
            )
    return surf.reshape(mesh[0].shape)


def fit_surf(mesh, func, indx=5, idx_file="surf_idx.txt"):
    X = mesh[0].flatten()
    Y = mesh[1].flatten()
    A = [np.ones_like(X)]
    for i in range(1, indx):
        for j in range(i + 1):
            ix, iy = j, i - j
            A.append(X**(ix) * Y**(iy))

    A = np.array(A).T
    B = func.flatten()
    coeff, r, rank, s = np.linalg.lstsq(A, B)
    print(coeff)
    print(r)
    print(s)

    xs, xe = mesh[0][0, 0], mesh[0][0, -1]
    ys, ye = mesh[1][0, 0], mesh[1][-1, 0]

    fp = open(idx_file, "w")
    fp.write('{:2d}\n'.format(indx - 1))
    fp.write('x\t{:.2f}\t{:.2f}\n'.format(xs, xe))
    fp.write('y\t{:.2f}\t{:.2f}\n'.format(ys, ye))
    fp.write('\n')
    fp.write('{}\t{}\t{}\n'.format("ix", "iy", "coeff"))
    for i in range(indx):
        for j in range(i + 1):
            ix, iy = j, i - j
            fp.write(
                '{:2d}\t{:2d}\t{:0.10E}\n'.format(
                    ix, iy, coeff[sum(range(i + 1)) + j]
                )
            )
    return coeff


def fit_surf_2d(mesh, func, div=1 / 4, sxy=[0, 0], indx=5, idx_file="surf_idx.txt"):
    indx = indx + 1
    X, Y, Z = mesh[0], mesh[1], func
    nx, ny = X.shape
    sx, sy = sxy
    xs, ys = mesh[0][0, 0], mesh[1][0, 0]
    xe, ye = mesh[0][0, -1], mesh[1][-1, 0]
    dx, dy = mesh[0][0, 1] - mesh[0][0, 0], mesh[1][1, 0] - mesh[1][0, 0]
    px, py = int((sy - ys) / dy), int((sx - xs) / dx)
    mx, my = int(nx * div), int(ny * div)
    X = X[mx:nx - mx, my:ny - my]
    Y = Y[mx:nx - mx, my:ny - my]
    Z = Z[mx:nx - mx, my:ny - my]
    print(sx, sy)
    print(X[0, 0], X[0, -1], X[-1, 0])
    print(Y[0, 0], Y[0, -1], Y[-1, 0])
    coeff = fit_surf([X, Y], Z, indx, idx_file)

    X = mesh[0].flatten()
    Y = mesh[1].flatten()
    A = [np.ones_like(X)]
    for i in range(1, indx):
        for j in range(i + 1):
            ix, iy = j, i - j
            A.append(X**(ix) * Y**(iy))
    A = np.array(A).T
    surf = np.zeros_like(X)
    for i, c in enumerate(coeff):
        surf += A[:, i] * c
    return surf.reshape(func.shape)


def fit_surf_2d_mat(mesh, func, div=1 / 4, indx=5):
    indx = indx + 1
    X, Y, Z = mesh[0], mesh[1], func
    nx, ny = X.shape
    mx, my = int(nx * div), int(ny * div)
    X = X[mx:nx - mx, my:ny - my]
    Y = Y[mx:nx - mx, my:ny - my]
    Z = Z[mx:nx - mx, my:ny - my]
    print(X[0, 0], X[0, -1], X[-1, 0])
    print(Y[0, 0], Y[0, -1], Y[-1, 0])
    coeff = fit_surf([X, Y], Z, indx, "surf_idx.txt")

    X = mesh[0].flatten()
    Y = mesh[1].flatten()
    A = [np.ones_like(X)]
    for i in range(1, indx):
        for j in range(i + 1):
            ix, iy = j, i - j
            A.append(X**(ix) * Y**(iy))
    A = np.array(A).T
    surf = np.zeros_like(X)
    for i, c in enumerate(coeff):
        surf += A[:, i] * c
    return surf.reshape(func.shape)


def fit_surf_1d(mesh, func, indx=5, filename="surf"):
    nx, ny = func.shape
    fx = np.poly1d(np.polyfit(
        mesh[0][int(nx / 2), :], func[int(nx / 2), :], indx))
    fy = np.poly1d(np.polyfit(
        mesh[1][:, int(ny / 2)], func[:, int(ny / 2)], indx))
    return fx(mesh[0]) + fy(mesh[1])


def gen_fit_surf_2d(mesh, func, g_mesh, div=1 / 4, indx=5, idx_file="surf_idx.txt"):
    indx = indx + 1
    X, Y, Z = mesh[0], mesh[1], func
    nx, ny = X.shape
    mx, my = int(nx * div), int(ny * div)
    X = X[mx:nx - mx, my:ny - my]
    Y = Y[mx:nx - mx, my:ny - my]
    Z = Z[mx:nx - mx, my:ny - my]
    coeff = fit_surf([X, Y], Z, indx, idx_file)

    X = g_mesh[0].flatten()
    Y = g_mesh[1].flatten()
    A = [np.ones_like(X)]
    for i in range(1, indx):
        for j in range(i + 1):
            ix, iy = j, i - j
            A.append(X**(ix) * Y**(iy))
    A = np.array(A).T
    surf = np.zeros_like(X)
    for i, c in enumerate(coeff):
        surf += A[:, i] * c
    g_surf = surf.reshape(g_mesh[0].shape)
    return g_surf


def grasp_sfc(mesh, surf, sfc_file="surf.sfc"):
    fp = open(sfc_file, "w")
    ny, nx = surf.shape
    xs, xe = mesh[0][0, 0], mesh[0][0, -1]
    ys, ye = mesh[1][0, 0], mesh[1][-1, 0]
    fp.write(" {} data \n".format(sfc_file))
    fp.write(" {:.2e} {:.2e} {:.2e} {:.2e}\n".format(xs, ys, xe, ye))
    fp.write(" {:d} {:d}\n".format(nx, ny))
    for ix in range(nx):
        for iy in range(ny):
            fp.write(" {:.5e} ".format(surf[iy, ix]))
        fp.write("\n")


def grad_pln(mesh, surf, sxy=[0, 0]):
    sx, sy = sxy
    xs, ys = mesh[0][0, 0], mesh[1][0, 0]
    xe, ye = mesh[0][0, -1], mesh[1][-1, 0]
    dx, dy = mesh[0][0, 1] - mesh[0][0, 0], mesh[1][1, 0] - mesh[1][0, 0]
    mx, my = int((sy - ys) / dy), int((sx - xs) / dx)
    grdx = np.gradient(mesh[0], axis=1)
    grdy = np.gradient(mesh[1], axis=0)
    surf_gx = np.gradient(surf, axis=1) / grdx
    surf_gy = np.gradient(surf, axis=0) / grdy
    ax, ay = surf_gx[mx, my], surf_gy[mx, my]
    dx, dy = np.arctan(ax), np.arctan(ay)
    return dx, dy


def rot_matrix(axs="x", deg=0.0):
    mat = np.zeros([3, 3])
    if axs == "x":
        ix, iy, iz = 1, 2, 0
    elif axs == "y":
        ix, iy, iz = 2, 0, 1
    elif axs == "z":
        ix, iy, iz = 0, 1, 2
    else:
        ix, iy, iz = 1, 2, 0
    mat[iz, iz] = 1
    mat[ix, ix] = np.cos(np.deg2rad(deg))
    mat[ix, iy] = -np.sin(np.deg2rad(deg))
    mat[iy, ix] = np.sin(np.deg2rad(deg))
    mat[iy, iy] = np.cos(np.deg2rad(deg))
    return mat


def idx_surf(mesh, idx_file="surf_idx.txt", indx=3):
    idx = int(getline(idx_file, 1).split()[0])
    mat = np.zeros((idx + 1, idx + 1))
    lines = open(idx_file, "r").readlines()
    for line in lines[5:]:
        dat = line.split()
        ix, iy = int(dat[0]), int(dat[1])
        mat[ix, iy] = float(dat[2])

    surf = np.zeros_like(mesh[0])
    for ix in range(indx):
        for iy in range(indx):
            surf += mat[ix, iy] * mesh[0]**ix * mesh[1]**iy
    return surf
