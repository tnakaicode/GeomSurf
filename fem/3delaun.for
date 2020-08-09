
c     ------   program 3delaun  --------------------------------------

c     purpose : the delaunay triangulation in 3d space
c                                      using variable super box
c     last modified : 21 Dec 2005

      implicit real*8(a-h,o-z)

      parameter(ktj=50,kte=200,err=1.d-14)

      dimension x(ktj+8),y(ktj+8),z(ktj+8),
     &          mtj(kte,4),jac(kte,4),vx(kte),vy(kte),vz(kte),
     &          rv(kte),kv(kte),istack(kte),map(kte)


c  -- input of data

      call input(ktj,node,x,y,z)

c  -- delaunay triangulation

      call tetgen(ktj,node,x,y,z,
     &            kte,nelm,mtj,jac,vx,vy,vz,rv,kv,istack,map,err)

c  -- output of results

      call data(kte,nelm,node,mtj,jac,x,y,z)


      stop
      end
c==
c     ------   subroutine input   -----------------------------------

c     purpose : input of data
c     last modified :  09 Sep 2005

      subroutine input(ktj,node,x,y,z)
      implicit real*8(a-h,o-z)

      dimension x(*),y(*),z(*)
      character fname*30


      write(6,600) ' Input File Name = ? '
      read(5,600) fname

      open(8,file=fname)

      read(8,*) node
      if ( node.gt.ktj ) then
        write(6,600) ' *** ERROR(sub.input) node > ktj'
        stop 900
      elseif ( node.lt.1 ) then
        write(6,600) ' *** ERROR(sub.input) node < 1'
        stop 900
      endif
      read(8,*) (x(i),y(i),z(i),i=1,node)

      close(8)


      return
  600 format(a)
      end
c==
c     ------   subroutine tetgen   -----------------------------------

c     purpose : the delaunay triangulation in 3-dimensional space
c     last modified : 21 Dec 2005

      subroutine tetgen(ktj,node,x,y,z,
     &                  kte,nelm,mtj,jac,vx,vy,vz,rv,kv,istack,map,err)
      implicit real*8(a-h,o-z)

      dimension x(*),y(*),z(*),
     &          mtj(kte,4),jac(kte,4),vx(*),vy(*),vz(*),rv(*),
     &          kv(*),istack(*),map(*)


c  -- computation of max & min coords for x,y,z

      xmin = x(1)
      xmax = xmin
      ymin = y(1)
      ymax = ymin
      zmin = z(1)
      zmax = zmin
      do 10 i=2,node
        xmin = dmin1(xmin,x(i))
        xmax = dmax1(xmax,x(i))
        ymin = dmin1(ymin,y(i))
        ymax = dmax1(ymax,y(i))
        zmin = dmin1(zmin,z(i))
        zmax = dmax1(zmax,z(i))
   10 continue
      rax  = dabs(xmax-xmin)
      ray  = dabs(ymax-ymin)
      raz  = dabs(zmax-zmin)
      rmax = dmax1(rax,ray,raz)
      if ( rmax.le.err ) then
        write(6,600) ' *** ERROR(sub.tetgen) rmax is zero'
        stop 900
      endif

c  -- normalization and shift to the positive region for x,y,z-coords

      rrm = 1.d0/rmax
      do 20 i=1,node
        x(i) = rrm*(x(i)-xmin)
        y(i) = rrm*(y(i)-ymin)
        z(i) = rrm*(z(i)-zmin)
   20 continue

      rax = rrm*rax
      ray = rrm*ray
      raz = rrm*raz

c  -- compute delaunay triangulation

      call delaun(ktj,kte,node,nelm,x,y,z,mtj,jac,vx,vy,vz,
     &            rv,kv,istack,map,err,rax,ray,raz)

c  -- reset x,y,z-coords to original values

      do 40 i=1,node
        x(i) = rmax*x(i)+xmin
        y(i) = rmax*y(i)+ymin
        z(i) = rmax*z(i)+zmin
   40 continue


      return
  600 format(a)
      end
c==
c     ------   subroutine delaun   -----------------------------------

c     purpose : compute 3d delaunay triangulation
c     last modified : 21 Dec 2005

      subroutine delaun(ktj,kte,node,nelm,x,y,z,mtj,jac,vx,vy,vz,
     &                  rv,kv,istack,map,err,rax,ray,raz)
      implicit real*8(a-h,o-z)

      dimension x(*),y(*),z(*),mtj(kte,4),jac(kte,4),
     &          vx(*),vy(*),vz(*),rv(*),kv(*),istack(*),map(*)

      data alpha/2.d0/


c  -- initialization

      do 10 i=1,kte
        map(i) = 0
   10 continue

c  -- translate to the original of this model

      xbar = 0.5d0*(alpha-1.d0)*rax
      ybar = 0.5d0*(alpha-1.d0)*ray
      zbar = 0.5d0*(alpha-1.d0)*raz

      do 100 i=1,node
        x(i) = x(i)+xbar
        y(i) = y(i)+ybar
        z(i) = z(i)+zbar
  100 continue

c  -- prepare for six tetrahedra dividing a super-cubic

      call box(ktj,kte,x,y,z,
     &         nelm,mtj,jac,vx,vy,vz,rv,alpha,
     &         rax,ray,raz)

c  -- loop over each point

      do 20 i=1,node

        ip = i
        xp = x(ip)
        yp = y(ip)
        zp = z(ip)

c  -- search for tetrahedron which includes new data point

        loc = locate(kte,xp,yp,zp,x,y,z,nelm,mtj,jac,err)

c  -- search for all tetrahedra whose circumspheres
c       enclose new data point

        iv  = 0
        msk = 0

        iv = iv+1
        kv(iv) = loc
        map(loc) = 1
        msk = msk+1
        istack(msk) = loc

   30   continue
        if ( msk.ne.0 ) then

          isk = istack(msk)
          msk = msk-1
          do 40 j=1,4
            jelm = jac(isk,j)
            if ( jelm.gt.kte ) then
              write(6,600) ' *** ERROR(sub.delaun) jelm > kte'
              stop 900
            endif
            if ( jelm.ne.0 ) then
              if ( map(jelm).ne.1 ) then
                rad = rv(jelm)*(1.d0+err)
                dst = vx(jelm)*vx(jelm)-2.d0*vx(jelm)*xp+xp*xp
                if ( dst.ge.rad ) goto 40
                dst = dst
     &               +vy(jelm)*vy(jelm)-2.d0*vy(jelm)*yp+yp*yp
                if ( dst.ge.rad ) goto 40
                dst = dst
     &               +vz(jelm)*vz(jelm)-2.d0*vz(jelm)*zp+zp*zp
                if ( dst.ge.rad ) goto 40
                iv = iv+1
                kv(iv) = jelm
                map(jelm) = 1
                msk = msk+1
                istack(msk) = jelm
              endif
            endif
   40     continue

          goto 30

        endif

c  -- triangulation of the polyhedron formed by tetrahedra
c                       whose circumspheres enclose new data point

        call poly(ktj,kte,ip,iv,kv,nelm,mtj,jac,vx,vy,vz,rv,x,y,z,
     &            map,err)

   20 continue

      if ( nelm.gt.kte ) then
        write(6,600) ' *** ERROR(sub.delaun) nelm > kte'
        stop 900
      endif

c  -- remove all tetrahedra which include the forming points
c                                        of the six tetrahedra

      iv = 0
      do 70 i=1,nelm
        kv(i) = 0
   70 continue

      do 80 i=1,nelm
        if ( mtj(i,1).gt.ktj .or.
     &       mtj(i,2).gt.ktj .or.
     &       mtj(i,3).gt.ktj .or.
     &       mtj(i,4).gt.ktj ) then
          iv = iv+1
          kv(iv) = i
        endif
   80 continue

      call remove(kte,iv,kv,nelm,mtj,jac,vx,vy,vz,rv,map)

c  -- examine results of mesh generation

      call fill(kte,nelm,mtj,jac)

c  -- return to original region

      do 110 i=1,node
        x(i) = x(i)-xbar
        y(i) = y(i)-ybar
        z(i) = z(i)-zbar
  110 continue


      return
  600 format(a)
      end
c==
c     ------   subroutine box      -----------------------------------

c     purpose : construct six tetrahedra dividing a variable-super-box
c                                          instead of supertetrahedron
c     last modified : 21 Dec 2005

      subroutine box(ktj,kte,x,y,z,
     &               nelm,mtj,jac,vx,vy,vz,rv,alpha,
     &               rax,ray,raz)
      implicit real*8(a-h,o-z)

      dimension x(*),y(*),z(*),mtj(kte,4),jac(kte,4),
     &          vx(*),vy(*),vz(*),rv(*)

      data zero/0.d0/


      nelm = 6

      xone = alpha*rax
      yone = alpha*ray
      zone = alpha*raz

      x(ktj+1) = zero
      y(ktj+1) = zero
      z(ktj+1) = zero
      x(ktj+2) = xone
      y(ktj+2) = zero
      z(ktj+2) = zero
      x(ktj+3) = xone
      y(ktj+3) = yone
      z(ktj+3) = zero
      x(ktj+4) = zero
      y(ktj+4) = yone
      z(ktj+4) = zero
      x(ktj+5) = zero
      y(ktj+5) = zero
      z(ktj+5) = zone
      x(ktj+6) = xone
      y(ktj+6) = zero
      z(ktj+6) = zone
      x(ktj+7) = xone
      y(ktj+7) = yone
      z(ktj+7) = zone
      x(ktj+8) = zero
      y(ktj+8) = yone
      z(ktj+8) = zone

      mtj(1,1) = ktj+2
      mtj(1,2) = ktj+7
      mtj(1,3) = ktj+5
      mtj(1,4) = ktj+6
      mtj(2,1) = ktj+1
      mtj(2,2) = ktj+2
      mtj(2,3) = ktj+3
      mtj(2,4) = ktj+5
      mtj(3,1) = ktj+2
      mtj(3,2) = ktj+3
      mtj(3,3) = ktj+5
      mtj(3,4) = ktj+7
      mtj(4,1) = ktj+5
      mtj(4,2) = ktj+4
      mtj(4,3) = ktj+8
      mtj(4,4) = ktj+7
      mtj(5,1) = ktj+1
      mtj(5,2) = ktj+3
      mtj(5,3) = ktj+4
      mtj(5,4) = ktj+5
      mtj(6,1) = ktj+4
      mtj(6,2) = ktj+3
      mtj(6,3) = ktj+7
      mtj(6,4) = ktj+5

      jac(1,1) = 0
      jac(1,2) = 0
      jac(1,3) = 0
      jac(1,4) = 3
      jac(2,1) = 3
      jac(2,2) = 5
      jac(2,3) = 0
      jac(2,4) = 0
      jac(3,1) = 6
      jac(3,2) = 1
      jac(3,3) = 0
      jac(3,4) = 2
      jac(4,1) = 0
      jac(4,2) = 0
      jac(4,3) = 6
      jac(4,4) = 0
      jac(5,1) = 6
      jac(5,2) = 0
      jac(5,3) = 2
      jac(5,4) = 0
      jac(6,1) = 3
      jac(6,2) = 4
      jac(6,3) = 5
      jac(6,4) = 0

      do 10 i=1,6
        ia = mtj(i,1)
        ib = mtj(i,2)
        ic = mtj(i,3)
        ip = mtj(i,4)
        v0 = volume(ia,ib,ic,ip,x,y,z)
        call sphere(ia,ib,ic,ip,v0,x,y,z,vx(i),vy(i),vz(i),rv(i))
   10 continue


      return
      end
c==
c     ------   function locate   -------------------------------------

c     purpose : locate tetrahedron which encloses new data point
c     last modified :  05 Sep 2005

      function locate(kte,xp,yp,zp,x,y,z,nelm,mtj,jac,err)
      implicit real*8(a-h,o-z)

      dimension x(*),y(*),z(*),mtj(kte,4),jac(kte,4)


      itet = nelm

   10 continue

      do 20 n=1,4

        i = mtj(itet,mod(n,4)+1)
        j = mtj(itet,4-(n-1)/2*2)
        k = mtj(itet,3-mod(n/2,2)*2)

        a = y(i)*z(j)+y(j)*z(k)+y(k)*z(i)
     &    -(y(i)*z(k)+y(j)*z(i)+y(k)*z(j))
        b = z(i)*x(j)+z(j)*x(k)+z(k)*x(i)
     &    -(z(i)*x(k)+z(j)*x(i)+z(k)*x(j))
        c = x(i)*y(j)+x(j)*y(k)+x(k)*y(i)
     &    -(x(i)*y(k)+x(j)*y(i)+x(k)*y(j))
        d = -(a*x(i)+b*y(i)+c*z(i))

        if ( a*xp+b*yp+c*zp+d.lt.-err ) then
          itet = jac(itet,n)
          goto 10
        endif

   20 continue

c  -- tetrahedron has been found

      locate = itet


      return
      end
c==
c     ------   subroutine poly   -------------------------------------

c     purpose : remove edges interior to the polyhedron
c                  and connect its vertices to new data point
c     last modified : 21 Dec 2005

      subroutine poly(ktj,kte,ip,iv,kv,nelm,mtj,jac,vx,vy,vz,rv,x,y,z,
     &                map,err)
      implicit real*8(a-h,o-z)
      parameter(ksf=50)

      dimension mtj(kte,4),jac(kte,4),vx(*),vy(*),vz(*),rv(*),
     &          kv(*),x(*),y(*),z(*),map(*),
     &          imen(ksf,3),jmen(ksf),kmen(ksf),vol(ksf)


c  -- search for boundary surfaces of the polyhedron

   10 continue
      ix = 0

      do 20 i=1,iv
        ielm = kv(i)

        do 30 j=1,4

          jelm = jac(ielm,j)
          ia = mtj(ielm,mod(j,4)+1)
          ib = mtj(ielm,4-(j-1)/2*2)
          ic = mtj(ielm,3-mod(j/2,2)*2)

          if ( jelm.eq.0 ) then
            ix = ix+1
            if ( ix.gt.ksf ) then
              write(6,600) ' *** ERROR(sub.poly) ix > ksf'
              stop 900
            endif
            imen(ix,1) = ia
            imen(ix,2) = ib
            imen(ix,3) = ic
            jmen(ix) = 0
            kmen(ix) = 0
            vol(ix)  = volume(ia,ib,ic,ip,x,y,z)
          elseif ( map(jelm).eq.0 ) then
            ix = ix+1
            if ( ix.gt.ksf ) then
              write(6,600) ' *** ERROR(sub.poly) ix > ksf'
              stop 900
            endif
            imen(ix,1) = ia
            imen(ix,2) = ib
            imen(ix,3) = ic
            jmen(ix) = jelm
            kmen(ix) = iface(kte,jelm,ielm,jac)
            vol(ix)  = volume(ia,ib,ic,ip,x,y,z)
            if ( vol(ix).le.err ) then
              iv = iv+1
              kv(iv) = jelm
              map(jelm) = 1
              goto 10
            endif
          endif

   30   continue

   20 continue

c  -- connect vertices of the surfaces to new data point

      ibound = ix
      do 40 i=iv+1,ibound
        nelm = nelm+1
        kv(i) = nelm
        map(nelm) = 1
   40 continue

      do 50 i=1,ibound
        map(kv(i)) = 0
   50 continue

      do 60 i=1,ibound

        ielm = kv(i)
        determ = vol(i)
        if ( dabs(determ).eq.0.d0 ) then
          write(6,600) ' *** ERROR(sub.poly) determ is zero'
          stop 900
        endif
        ia = imen(i,1)
        ib = imen(i,2)
        ic = imen(i,3)
        mtj(ielm,1) = ia
        mtj(ielm,2) = ib
        mtj(ielm,3) = ic
        mtj(ielm,4) = ip
        jac(ielm,4) = jmen(i)
        if ( jmen(i).ne.0 ) then
          jac(jmen(i),kmen(i)) = ielm
        endif
        call sphere(ia,ib,ic,ip,determ,x,y,z,xv,yv,zv,rr)
        vx(ielm) = xv
        vy(ielm) = yv
        vz(ielm) = zv
        rv(ielm) = rr

   60 continue

c  -- connect tetrahedra in polyhedron each other

      ix = 0
      do 70 i=1,ibound
        ielm = kv(i)

        do 80 j=1,3
          ia = mtj(ielm,mod(j,3)+1)
          ib = mtj(ielm,mod(mod(j,3)+1,3)+1)

          do 90 k=1,ix
            ja = imen(k,1)
            jb = imen(k,2)
            if ( ia.eq.ja .and. ib.eq.jb ) then
              jac(ielm,j) = jmen(k)
              jac(jmen(k),kmen(k)) = ielm
              imen(k,1) = imen(ix,1)
              imen(k,2) = imen(ix,2)
              jmen(k) = jmen(ix)
              kmen(k) = kmen(ix)
              ix = ix-1
              goto 80
            endif
   90     continue

          ix = ix+1
          imen(ix,1) = ib
          imen(ix,2) = ia
          jmen(ix) = ielm
          kmen(ix) = j
   80   continue

   70 continue

c  -- in case that new tetrahedra is less than old tetrahedra

      if ( iv.gt.ibound ) then
        ir = iv-ibound
        do 100 i=1,ir
          kv(i) = kv(ibound+i)
          map(kv(i)) = kv(i)
  100   continue

        call qsorti(kte,ir,kv,map)

        do 110 i=1,ir
          ielm = kv(ir-i+1)
          map(ielm) = 0

          if ( ielm.ne.nelm ) then
            vx(ielm) = vx(nelm)
            vy(ielm) = vy(nelm)
            vz(ielm) = vz(nelm)
            rv(ielm) = rv(nelm)
            do 120 j=1,4
              mtj(ielm,j) = mtj(nelm,j)
              jelm = jac(nelm,j)
              jac(ielm,j) = jelm
              if ( jelm.ne.0 ) then
                jac(jelm,iface(kte,jelm,nelm,jac)) = ielm
              endif
  120       continue
          endif

          nelm = nelm-1
  110   continue

      endif


      return
  600 format(a)
      end
c==
c     ------   subroutine qsorti   -------------------------------------

c     purpose : quick sorting
c     last modified :  21 Nov 2005

c     maxstk: maximum stack size

      subroutine qsorti(k,n,list,key)
      implicit real*8(a-h,o-z)

      parameter(maxstk=32)

      dimension list(k),key(k),ilst(maxstk),irst(maxstk)

      ll = 1
      lr = n
      istk = 0
   10 continue
      if ( ll.lt.lr ) then
        nl = ll
        nr = lr
        lm = (ll+lr)/2
        iguess = key(list(lm))

c --  find keys for exchange

   20   continue
        if ( key(list(nl)).lt.iguess ) then
          nl = nl+1
          goto 20
        endif

   30   continue
        if ( iguess.lt.key(list(nr)) ) then
          nr = nr-1
          goto 30
        endif

        if ( nl.lt.(nr-1) ) then
          ltemp = list(nl)
          list(nl) = list(nr)
          list(nr) = ltemp
          nl = nl+1
          nr = nr-1
          goto 20
        endif

c --  deal with crossing of pointers

        if ( nl.le.nr ) then
          if ( nl.lt.nr ) then
            ltemp = list(nl)
            list(nl) = list(nr)
            list(nr) = ltemp
          endif
          nl = nl+1
          nr = nr-1
        endif

c --  select sub-list to be processed next

        istk = istk+1
        if ( istk.gt.maxstk ) then
          write(6,600) ' *** ERROR(sub.qsorti) istk > maxstk'
          stop 900
        endif
        if ( nr.lt.lm ) then
          ilst(istk) = nl
          irst(istk) = lr
          lr = nr
        else
          ilst(istk) = ll
          irst(istk) = nr
          ll = nl
        endif
        goto 10

      endif

c --  process any stacked sub-lists

      if ( istk.ne.0 ) then
        ll = ilst(istk)
        lr = irst(istk)
        istk = istk-1
        goto 10
      endif


      return
  600 format(a)
      end
c==
c     ------   function volume   -------------------------------------

c     purpose : computation of volume of tetrahedron
c     last modified : 05 Sep 2005

      function volume(ia,ib,ic,ip,x,y,z)
      implicit real*8(a-h,o-z)

      dimension x(*),y(*),z(*)


      xa = x(ia)
      ya = y(ia)
      za = z(ia)
      xb = x(ib)
      yb = y(ib)
      zb = z(ib)
      xc = x(ic)
      yc = y(ic)
      zc = z(ic)
      xp = x(ip)
      yp = y(ip)
      zp = z(ip)

      va = xb*yc*zp+xa*ya*zp+xb*ya*za+xa*yc*za
     &   -(xb*yc*za+xa*ya*za+xb*ya*zp+xa*yc*zp)
      vb = yb*zc*xp+ya*za*xp+yb*za*xa+ya*zc*xa
     &   -(yb*zc*xa+ya*za*xa+yb*za*xp+ya*zc*xp)
      vc = zb*xc*yp+za*xa*yp+zb*xa*ya+za*xc*ya
     &   -(zb*xc*ya+za*xa*ya+zb*xa*yp+za*xc*yp)

      wa = xb*zc*ya+xa*za*ya+xb*za*yp+xa*zc*yp
     &   -(xb*zc*yp+xa*za*yp+xb*za*ya+xa*zc*ya)
      wb = yb*xc*za+ya*xa*za+yb*xa*zp+ya*xc*zp
     &   -(yb*xc*zp+ya*xa*zp+yb*xa*za+ya*xc*za)
      wc = zb*yc*xa+za*ya*xa+zb*ya*xp+za*yc*xp
     &   -(zb*yc*xp+za*ya*xp+zb*ya*xa+za*yc*xa)

      volume = va+vb+vc+wa+wb+wc


      return
      end
c==
c     ------   subroutine sphere   -----------------------------------

c     purpose : computation of circumsphere of tetrahedron
c     last modified : 05 Sep 2005

      subroutine sphere(ia,ib,ic,ip,determ,x,y,z,xv,yv,zv,rr)
      implicit real*8(a-h,o-z)

      dimension x(*),y(*),z(*)


      xa = x(ia)
      ya = y(ia)
      za = z(ia)
      xb = x(ib)
      yb = y(ib)
      zb = z(ib)
      xc = x(ic)
      yc = y(ic)
      zc = z(ic)
      xp = x(ip)
      yp = y(ip)
      zp = z(ip)

c  -- cofactor

      p11 = yc*zp+ya*za+yp*za+ya*zc-(yc*za+ya*zp+yp*zc+ya*za)
      p12 = xp*zc+xa*za+xc*za+xa*zp-(xp*za+xa*zc+xc*zp+xa*za)
      p13 = xc*yp+xa*ya+xp*ya+xa*yc-(xc*ya+xa*yp+xp*yc+xa*ya)
      p21 = yp*zb+ya*za+yb*za+ya*zp-(yp*za+ya*zb+yb*zp+ya*za)
      p22 = xb*zp+xa*za+xp*za+xa*zb-(xb*za+xa*zp+xp*zb+xa*za)
      p23 = xp*yb+xa*ya+xb*ya+xa*yp-(xp*ya+xa*yb+xb*yp+xa*ya)
      p31 = yb*zc+ya*za+yc*za+ya*zb-(yb*za+ya*zc+yc*zb+ya*za)
      p32 = xc*zb+xa*za+xb*za+xa*zc-(xc*za+xa*zb+xb*zc+xa*za)
      p33 = xb*yc+xa*ya+xc*ya+xa*yb-(xb*ya+xa*yc+xc*yb+xa*ya)

      xyza = xa*xa+ya*ya+za*za
      aa = 0.5d0*(xb*xb+yb*yb+zb*zb-xyza)
      bb = 0.5d0*(xc*xc+yc*yc+zc*zc-xyza)
      cc = 0.5d0*(xp*xp+yp*yp+zp*zp-xyza)
      xx = p11*aa+p21*bb+p31*cc
      yy = p12*aa+p22*bb+p32*cc
      zz = p13*aa+p23*bb+p33*cc
      xv = xx/determ
      yv = yy/determ
      zv = zz/determ

      rr = xa*xa+xv*xv+ya*ya+yv*yv+za*za+zv*zv
     &    -2.d0*(xa*xx+ya*yy+za*zz)/determ


      return
      end
c==
c     ------   subroutine remove   -----------------------------------

c     purpose : remove unnecessary tetrahedra
c     last modified :  31 Aug 2005

      subroutine remove(kte,iv,kv,nelm,mtj,jac,vx,vy,vz,rv,map)
      implicit real*8(a-h,o-z)

      dimension kv(*),mtj(kte,4),jac(kte,4),vx(*),vy(*),vz(*),
     &          rv(*),map(*)


c  -- initialization

      m = 0
      n = 0
      do 10 i=1,nelm
        map(i) = 1
   10 continue
      do 20 i=1,iv
        map(kv(i)) = 0
   20 continue

      do 30 i=1,nelm
        if ( map(i).ne.0 ) then
          m = m+1
          map(i) = m
        endif
   30 continue

      do 40 i=1,nelm
        if ( map(i).ne.0 ) then
          n = n+1
          vx(n) = vx(i)
          vy(n) = vy(i)
          vz(n) = vz(i)
          rv(n) = rv(i)
          do 50 ia=1,4
            mtj(n,ia) = mtj(i,ia)
            if ( jac(i,ia).eq.0 ) then
              jac(n,ia) = 0
            else
              jac(n,ia) = map(jac(i,ia))
            endif
   50     continue
        endif
   40 continue

      do 60 i=n+1,nelm
        vx(i) = 0.d0
        vy(i) = 0.d0
        vz(i) = 0.d0
        rv(i) = 0.d0
        do 70 ia=1,4
          mtj(i,ia) = 0
          jac(i,ia) = 0
   70   continue
   60 continue

      nelm = nelm-iv


      return
      end
c==
c     ------   subroutine fill   -------------------------------------

c     purpose : check if the domain is filled by tetrahedra
c     last modified : 29 Aug 2005

      subroutine fill(kte,nelm,mtj,jac)
      implicit real*8(a-h,o-z)

      dimension mtj(kte,4),jac(kte,4)


      do 10 i=1,nelm
        ielm = i

        do 20 j=1,4
          ia = mtj(ielm,mod(j,4)+1)
          ib = mtj(ielm,4-(j-1)/2*2)
          ic = mtj(ielm,3-mod(j/2,2)*2)
          jelm = jac(ielm,j)
          if ( jelm.ne.0 .and. ielm.lt.jelm ) then
            k = iface(kte,jelm,ielm,jac)
            ja = mtj(jelm,mod(k,4)+1)
            jb = mtj(jelm,4-(k-1)/2*2)
            jc = mtj(jelm,3-mod(k/2,2)*2)
            if ( ia.eq.jc .and. ib.eq.jb .and. ic.eq.ja ) goto 20
            if ( ib.eq.jc .and. ic.eq.jb .and. ia.eq.ja ) goto 20
            if ( ic.eq.jc .and. ia.eq.jb .and. ib.eq.ja ) goto 20
            write(6,600)
     &        ' *** ERROR(sub.fill) no adjacent tetrahedron'
            stop 900
          endif
   20   continue

   10 continue


      return
  600 format(a)
      end
c==
c     ------   function iface   --------------------------------------

c     purpose : find surface in tet(l) which is adjacent to tet(k)
c     last modified : 29 Aug 2005

      function iface(kte,l,k,jac)
      implicit real*8(a-h,o-z)

      dimension jac(kte,4)


      do 10 i=1,4
        if ( jac(l,i).eq.k ) then
          iface = i
          goto 20
        endif
   10 continue

      write(6,600)
     &     ' *** ERROR(func.iface) no adjacent tetrahedron'
      stop 900

   20 continue


      return
  600 format(a)
      end
c==
c     ------   subroutine data   ---------------------------------------

c     purpose : print results on data file
c     last modified : 21 Dec 2005

      subroutine data(kte,nelm,node,mtj,jac,x,y,z)
      implicit real*8(a-h,o-z)

      dimension mtj(kte,4),jac(kte,4),x(*),y(*),z(*)
      character fname*30


      write(6,*)
      write(6,600) ' Output File Name = ?'
      read(5,600) fname

      open(9,file=fname)

      write(9,600) 'nelm,node'
      write(9,610) nelm,node

      write(9,600) 'e-id,mtj,jac'
      do 10 i=1,nelm
        write(9,620) i,(mtj(i,j),j=1,4),(jac(i,j),j=1,4)
   10 continue

      write(9,600) 'n-id,x,y,z'
      do 20 i=1,node
        write(9,630) i,x(i),y(i),z(i)
   20 continue

      close(9)


      return
  600 format(a)
  610 format(i7,',',i7)
  620 format(i7,8(',',i7))
  630 format(i7,',',1pe15.7,',',e15.7,',',e15.7)
      end
