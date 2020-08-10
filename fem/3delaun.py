import numpy as np
import matplotlib.pyplot as plt
import sys
from linecache import getline, clearcache


def tetgen(ktj, node, x, y, z, kte, nelm, mtj, jac, vx, vy, vz, rv, kv, istack, map, err):
    #     ------   subroutine tetgen   -----------------------------------
    #     purpose : the delaunay triangulation in 3-dimensional space
    #     last modified : 21 Dec 2005
    xmin, ymin, zmin = np.min(x), np.min(y), np.min(z)
    xmax, ymax, zmax = np.max(x), np.max(y), np.max(z)
    rax = np.abs(xmax - xmin)
    ray = np.abs(ymax - ymin)
    raz = np.abs(zmax - zmin)
    rmax = np.max([rax, ray, raz])
    rrm = 1.0 / rmax

    #  -- normalization and shift to the positive region for x,y,z-coords
    x = rrm * (x - xmin)
    y = rrm * (y - ymin)
    z = rrm * (z - zmin)
    rax = rrm * rax
    ray = rrm * ray
    raz = rrm * raz

    #  -- compute delaunay triangulation
    delaun(ktj, kte, node, nelm, x, y, z, mtj, jac, vx,
           vy, vz, rv, kv, istack, map, err, rax, ray, raz)

    #  -- reset x,y,z-coords to original values
    x = rmax * x + xmin
    y = rmax * y + ymin
    z = rmax * z + zmin


def delaun(ktj, kte, node, nelm, x, y, z, mtj, jac, vx, vy, vz, rv, kv, istack, map, err, rax, ray, raz):
    #     ------   subroutine delaun   -----------------------------------
    #     purpose : compute 3d delaunay triangulation
    #     last modified : 21 Dec 2005
    #
    alpha = 2.0
    mmap = np.zeros(kte)

    #  -- translate to the original of this model
    xbar = 0.5 * (alpha - 1.0) * rax
    ybar = 0.5 * (alpha - 1.0) * ray
    zbar = 0.5 * (alpha - 1.0) * raz
    x = x + xbar
    y = y + ybar
    z = z + zbar

    #  -- prepare for six tetrahedra dividing a super-cubic
    box(ktj, kte, x, y, z, nelm, mtj, jac, vx, vy, vz, rv, alpha, rax, ray, raz)

    #  -- loop over each point

    for i in range(node):
        ip = i
        xp = x[ip]
        yp = y[ip]
        zp = z[ip]

        #  -- search for tetrahedron which includes new data point
        loc = locate(kte, xp, yp, zp, x, y, z, nelm, mtj, jac, err)

        #  -- search for all tetrahedra whose circumspheres
        #       enclose new data point
        iv = 0
        msk = 0

        iv = iv + 1
        kv[iv] = loc
        mmap[loc] = 1
        msk = msk + 1
        istack[msk] = loc

        if (msk != 0):
            isk = istack[msk]
            msk = msk - 1
            for j in range(4):
                jelm = jac[isk, j]
                if (jelm >= kte):
                    break
                if (jelm != 0):
                    if (mmap[jelm] != 1):
                        rad = rv[jelm] * (1.0 + err)
                        dst = vx[jelm] * vx[jelm] - \
                            2.0 * vx[jelm] * xp + xp * xp
                        if (dst >= rad):
                            break
                        dst = dst + vy[jelm] * vy[jelm] - \
                            2.0 * vy[jelm] * yp + yp * yp
                        if (dst >= rad):
                            break
                        dst = dst + vz[jelm] * vz[jelm] - \
                            2.0 * vz[jelm] * zp + zp * zp
                        if (dst >= rad):
                            break
                        iv = iv + 1
                        kv[iv] = jelm
                        mmap[jelm] = 1
                        msk = msk + 1
                        istack[msk] = jelm

        #  -- triangulation of the polyhedron formed by tetrahedra
        #                       whose circumspheres enclose new data point

        poly(ktj, kte, ip, iv, kv, nelm, mtj, jac,
             vx, vy, vz, rv, x, y, z, mmap, err)

    #  -- remove all tetrahedra which include the forming points
    #                                        of the six tetrahedra

    iv = 0
    kv = np.zeros(nelm)
    for i in range(nelm):
        if (mtj[i, 0] >= ktj or mtj[i, 1] >= ktj or mtj[i, 2] >= ktj or mtj[i, 3] >= ktj):
            iv += 1
            kv[iv] = i
    remove(kte, iv, kv, nelm, mtj, jac, vx, vy, vz, rv, mmap)

    #  -- examine results of mesh generation
    fill(kte, nelm, mtj, jac)

    #  -- return to original region
    x = x - xbar
    y = y - ybar
    z = z - zbar


if __name__ == '__main__':
    print("3delaun")
    fname = "3delaun.input.txt"
    node = int(getline(fname, 1).split()[0])
    data = np.loadtxt(fname, skiprows=1, delimiter=",")
    x, y, z = data[0, :], data[1, :], data[2, :]
