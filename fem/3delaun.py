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


def box(ktj, kte, x, y, z, nelm, mtj, jac, vx, vy, vz, rv, alpha, rax, ray, raz):
    # ==
    #     ------   subroutine box      -----------------------------------
    #   purpose : construct six tetrahedra dividing a variable-super-box
    #                                          instead of supertetrahedron
    #     last modified : 21 Dec 2005
    zero = 0.0
    nelm = 6
    xone = alpha * rax
    yone = alpha * ray
    zone = alpha * raz
    x[ktj + 1] = zero
    y[ktj + 1] = zero
    z[ktj + 1] = zero
    x[ktj + 2] = xone
    y[ktj + 2] = zero
    z[ktj + 2] = zero
    x[ktj + 3] = xone
    y[ktj + 3] = yone
    z[ktj + 3] = zero
    x[ktj + 4] = zero
    y[ktj + 4] = yone
    z[ktj + 4] = zero
    x[ktj + 5] = zero
    y[ktj + 5] = zero
    z[ktj + 5] = zone
    x[ktj + 6] = xone
    y[ktj + 6] = zero
    z[ktj + 6] = zone
    x[ktj + 7] = xone
    y[ktj + 7] = yone
    z[ktj + 7] = zone
    x[ktj + 8] = zero
    y[ktj + 8] = yone
    z[ktj + 8] = zone
    mtj[1, 1] = ktj + 2
    mtj[1, 2] = ktj + 7
    mtj[1, 3] = ktj + 5
    mtj[1, 4] = ktj + 6
    mtj[2, 1] = ktj + 1
    mtj[2, 2] = ktj + 2
    mtj[2, 3] = ktj + 3
    mtj[2, 4] = ktj + 5
    mtj[3, 1] = ktj + 2
    mtj[3, 2] = ktj + 3
    mtj[3, 3] = ktj + 5
    mtj[3, 4] = ktj + 7
    mtj[4, 1] = ktj + 5
    mtj[4, 2] = ktj + 4
    mtj[4, 3] = ktj + 8
    mtj[4, 4] = ktj + 7
    mtj[5, 1] = ktj + 1
    mtj[5, 2] = ktj + 3
    mtj[5, 3] = ktj + 4
    mtj[5, 4] = ktj + 5
    mtj[6, 1] = ktj + 4
    mtj[6, 2] = ktj + 3
    mtj[6, 3] = ktj + 7
    mtj[6, 4] = ktj + 5
    jac[1, 1] = 0
    jac[1, 2] = 0
    jac[1, 3] = 0
    jac[1, 4] = 3
    jac[2, 1] = 3
    jac[2, 2] = 5
    jac[2, 3] = 0
    jac[2, 4] = 0
    jac[3, 1] = 6
    jac[3, 2] = 1
    jac[3, 3] = 0
    jac[3, 4] = 2
    jac[4, 1] = 0
    jac[4, 2] = 0
    jac[4, 3] = 6
    jac[4, 4] = 0
    jac[5, 1] = 6
    jac[5, 2] = 0
    jac[5, 3] = 2
    jac[5, 4] = 0
    jac[6, 1] = 3
    jac[6, 2] = 4
    jac[6, 3] = 5
    jac[6, 4] = 0

    for i in range(6):
        ia = mtj[i, 1]
        ib = mtj[i, 2]
        ic = mtj[i, 3]
        ip = mtj[i, 4]
        v0 = volume[ia, ib, ic, ip, x, y, z]
        sphere(ia, ib, ic, ip, v0, x, y, z, vx[i], vy[i], vz[i], rv[i])


def locate(kte, xp, yp, zp, x, y, z, nelm, mtj, jac, err):
    # ==
    #     ------   function locate   -------------------------------------
    #     purpose : locate tetrahedron which encloses new data point
    #     last modified :  05 Sep 2005
    itet = nelm

    for n in range(4):
        i = mtj(itet, mod[n, 4] + 1)
        j = mtj(itet, 4 - (n - 1) / 2 * 2)
        k = mtj(itet, 3 - mod[int(n / 2), 2] * 2)

        a = y(i) * z(j) + y(j) * z(k) + y(k) * z(i) - \
            (y(i) * z(k) + y(j) * z(i) + y(k) * z(j))
        b = z(i) * x(j) + z(j) * x(k) + z(k) * x(i) - \
            (z(i) * x(k) + z(j) * x(i) + z(k) * x(j))
        c = x(i) * y(j) + x(j) * y(k) + x(k) * y(i) - \
            (x(i) * y(k) + x(j) * y(i) + x(k) * y(j))
        d = -(a * x(i) + b * y(i) + c * z(i))

        if (a * xp + b * yp + c * zp + d < -err):
            itet = jac[itet, n]
            break

    #  -- tetrahedron has been found
    locate = itet


def poly(ktj, kte, ip, iv, kv, nelm, mtj, jac, vx, vy, vz, rv, x, y, z, mmap, err):
    # ==
    #     ------   subroutine poly   -------------------------------------
    #     purpose : remove edges interior to the polyhedron
    #                  and connect its vertices to new data point
    #     last modified : 21 Dec 2005
    ix = 0
    for i in range(iv):
        ielm = kv[i]
        for j in range(4):
            jelm = jac[ielm, j]
            ia = mtj[ielm, mod[j, 4] + 1]
            ib = mtj[ielm, 4 - (j - 1) / 2 * 2]
            ic = mtj[ielm, 3 - mod[j / 2, 2] * 2]

            if (jelm == 0):
                ix = ix + 1
                imen[ix, 1] = ia
                imen[ix, 2] = ib
                imen[ix, 3] = ic
                jmen[ix] = 0
                kmen[ix] = 0
                vol[ix] = volume[ia, ib, ic, ip, x, y, z]
            elif (mmap[jelm] == 0):
                ix = ix + 1
                imen[ix, 1] = ia
                imen[ix, 2] = ib
                imen[ix, 3] = ic
                jmen[ix] = jelm
                kmen[ix] = iface[kte, jelm, ielm, jac]
                vol[ix] = volume[ia, ib, ic, ip, x, y, z]
                if (vol[ix] <= err):
                    iv = iv + 1
                    kv[iv] = jelm
                    mmap[jelm] = 1
                    break

    #  -- connect vertices of the surfaces to new data point
    ibound = ix
    for i in range(iv + 1, ibound):
        nelm = nelm + 1
        kv[i] = nelm
        mmap[nelm] = 1
    for i in range(1, ibound):
        mmap[kv[i]] = 0

    for i in range(1, ibound):
        ielm = kv[i]
        determ = vol[i]
        ia = imen[i, 1]
        ib = imen[i, 2]
        ic = imen[i, 3]
        mtj[ielm, 1] = ia
        mtj[ielm, 2] = ib
        mtj[ielm, 3] = ic
        mtj[ielm, 4] = ip
        jac[ielm, 4] = jmen[i]
        if (jmen(i) != 0):
            jac[jmen[i], kmen[i]] = ielm
        sphere(ia, ib, ic, ip, determ, x, y, z, xv, yv, zv, rr)
        vx[ielm] = xv
        vy[ielm] = yv
        vz[ielm] = zv
        rv[ielm] = rr

    #  -- connect tetrahedra in polyhedron each other
    ix = 0
    for i in range(1, ibound):
        ielm = kv(i)
        for j in range(3):
            ia = mtj(ielm, mod(j, 3) + 1)
            ib = mtj(ielm, mod(mod(j, 3) + 1, 3) + 1)
            for k in range(1, ix):
                ja = imen(k, 1)
                jb = imen(k, 2)
                if (ia == ja and ib == jb):
                    jac(ielm, j) = jmen(k)
                    jac(jmen(k), kmen(k)) = ielm
                    imen(k, 1) = imen(ix, 1)
                    imen(k, 2) = imen(ix, 2)
                    jmen(k) = jmen(ix)
                    kmen(k) = kmen(ix)
                    ix = ix - 1
                    break
            ix = ix + 1
            imen(ix, 1) = ib
            imen(ix, 2) = ia
            jmen(ix) = ielm
            kmen(ix) = j

    #  -- in case that new tetrahedra is less than old tetrahedra
    if (iv > ibound):
        ir = iv - ibound
        for i in range(1, ir):
            kv(i) = kv(ibound + i)
            mmap(kv(i)) = kv(i)

        qsorti(kte, ir, kv, map)

        for i in range(1, ir):
            ielm = kv(ir - i + 1)
            mmap(ielm) = 0

            if (ielm != nelm):
                vx(ielm) = vx(nelm)
                vy(ielm) = vy(nelm)
                vz(ielm) = vz(nelm)
                rv(ielm) = rv(nelm)
                for j in range(1, 4):
                    mtj(ielm, j) = mtj(nelm, j)
                    jelm = jac(nelm, j)
                    jac(ielm, j) = jelm
                    if (jelm != 0):
                        jac(jelm, iface(kte, jelm, nelm, jac)) = ielm

            nelm = nelm - 1


if __name__ == '__main__':
    print("3delaun")
    fname = "3delaun.input.txt"
    node = int(getline(fname, 1).split()[0])
    data = np.loadtxt(fname, skiprows=1, delimiter=",")
    x, y, z = data[0, :], data[1, :], data[2, :]
