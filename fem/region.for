c     ------   program region   --------------------------------------
c     purpose : remove tetrahedra out of the region
c       assumption : The results of delaunay does not break 
c                                           the boundary surface.
c     last modified : 2006.01.01

      implicit real*8(a-h,o-z)
      parameter(ktj=100,kte=500)


      dimension mtj(kte,4),jac(kte,4),x(ktj),y(ktj),z(ktj),
     &          mtri(kte,3),xtri(ktj),ytri(ktj),ztri(ktj),
     &          itet(kte),itri(kte),neid(kte)


c  -- data input

      call input(kte,ktj,nelm,node,mtj,jac,x,y,z,
     &           netri,ndtri,mtri,xtri,ytri,ztri)

c  -- decide the tetrahedron in/out

      call decide(kte,ktj,nelm,node,mtj,jac,x,y,z,
     &            netri,ndtri,mtri,xtri,ytri,ztri,
     &            itet,itri,new,neid)

c  -- data output

      call data(kte,nelm,node,mtj,jac,x,y,z,itet,new,neid)


      stop
      end
c==
c**********   subroutine input   ************************************

c     purpose : data input from 3delaun and the surface triangles
c     last modified : 2006.01.01

      subroutine input(kte,ktj,nelm,node,mtj,jac,x,y,z,
     &                 netri,ndtri,mtri,xtri,ytri,ztri)
      implicit real*8(a-h,o-z)

      dimension mtj(kte,4),jac(kte,4),x(*),y(*),z(*),
     &          mtri(kte,3),xtri(*),ytri(*),ztri(*)
      character fname*30


c  -- read tetrahedra data ( from 3delaun result )

      write(*,*)
      write(*,500)
      read(*,*) fname
      open(8,file=fname)

      read(8,*)
      read(8,*) nelm,node
      if ( nelm.gt.kte ) then
        write(6,600) ' *** ERROR(sub.input) nelm > kte'
        stop 900
      elseif ( node.gt.ktj ) then
        write(6,600) ' *** ERROR(sub.input) node > ktj'
        stop 900
      endif
      read(8,*)
      do 10 i=1,nelm
        read(8,*) id,(mtj(i,j),j=1,4),(jac(i,j),j=1,4)
   10 continue
      read(8,*)
      do 20 i=1,node
        read(8,*) id,x(i),y(i),z(i)
   20 continue

      close(8)

c  -- read triangles data

      write(*,510)
      read(*,*) fname
      open(8,file=fname)

      read(8,*) netri,ndtri
      if ( netri.gt.kte ) then
        write(6,600) ' *** ERROR(sub.input) netri > kte'
        stop 900
      elseif ( ndtri.gt.ktj ) then
        write(6,600) ' *** ERROR(sub.input) ndtri > ktj'
        stop 900
      endif
      do 50 i=1,netri
        read(8,*) id,(mtri(i,j),j=1,3)
   50 continue
      do 60 i=1,ndtri
        read(8,*) id,xtri(i),ytri(i),ztri(i)
   60 continue

      close(8)


      return
  500 format(' Input tetrahedra file name = ? ')
  510 format(/' Input triangles file name = ? ')
  600 format(a)
      end
c==
c**********   subroutine decide   ************************************

c     purpose : decide whether a tetrahedron is in the region or not
c     last modified : 2006.01.01

      subroutine decide(kte,ktj,nelm,node,mtj,jac,x,y,z,
     &                  netri,ndtri,mtri,xtri,ytri,ztri,
     &                  itet,itri,new,neid)
      implicit real*8(a-h,o-z)

      dimension mtj(kte,4),jac(kte,4),x(*),y(*),z(*),
     &          mtri(kte,3),xtri(*),ytri(*),ztri(*),
     &          itet(*),itri(*),neid(*)

c -- itet(i) =  1; outside
c -- itet(i) = -1; inside


      do 10 i=1,nelm
        itet(i) = 0
   10 continue

      do 20 i=1,netri
        itri(i) = 0
   20 continue

      call coinface(kte,nelm,mtj,netri,mtri,itet,itri)

      npick = 0
      do 300 i=1,netri
        if ( itri(i).ne.0 ) npick = npick+1
  300 continue

      if ( npick.lt.netri ) then
        call coinedge(kte,nelm,mtj,x,y,z,
     &                netri,mtri,xtri,ytri,ztri,itet,itri)
      endif

      npick = 0
      do 310 i=1,netri
        if ( itri(i).ne.0 ) npick = npick+1
  310 continue

      if ( npick.lt.netri ) then
        call coinnode(kte,nelm,mtj,x,y,z,
     &                netri,mtri,xtri,ytri,ztri,itet,itri)
      endif

  120 continue

      jcnt = 0
      do 200 i=1,nelm
        if ( itet(i).eq.0 ) then
          do 210 j=1,4
            if ( jac(i,j).ne.0 ) then
              if ( itet(jac(i,j)).ne.0 ) then
                itet(i) = itet(jac(i,j))
                jcnt = jcnt+1
                goto 200
              endif
            endif
  210     continue
        else
          jcnt = jcnt+1
        endif
  200 continue

      if ( jcnt.lt.nelm ) goto 120

      jchk = 0
      do 30 i=1,nelm
        if ( itet(i).eq.0 ) then
          write(6,600) i
          jchk = 1
        endif
   30 continue

      if ( jchk.eq.1 ) stop 900

      new = 0
      do 40 i=1,nelm
        if ( itet(i).eq.-1 ) then
          new = new+1
          neid(i) = new
        else
          neid(i) = 0
        endif
   40 continue


      return
  600 format(' *** ERROR in sub.decide.',
     &       ' The tetrahedron cannot find adjacent surfaces.'/
     &       '  Tet-id = ',i6)
      end
c==
c**********   subroutine coinface ************************************

c     purpose : pick up the tetrahedron which coincides 
c                                           with the surface triangle
c     last modified : 2006.01.01

      subroutine coinface(kte,nelm,mtj,netri,mtri,itet,itri)
      implicit real*8(a-h,o-z)

      dimension mtj(kte,4),mtri(kte,3),itet(*),itri(*)


      do 100 i=1,netri

        i1 = mtri(i,1)
        i2 = mtri(i,2)
        i3 = mtri(i,3)

        do 110 j=1,nelm

          if ( i1.eq.mtj(j,1) ) then
            if ( i2.eq.mtj(j,2) ) then
              if ( i3.eq.mtj(j,3) ) then
                itet(j) = 1
                itri(i) = itri(i)+1
              elseif ( i3.eq.mtj(j,4) ) then
                itet(j) = -1
                itri(i) = itri(i)+1
              endif
            elseif ( i2.eq.mtj(j,3) ) then
              if ( i3.eq.mtj(j,2) ) then
                itet(j) = -1
                itri(i) = itri(i)+1
              elseif ( i3.eq.mtj(j,4) ) then
                itet(j) = 1
                itri(i) = itri(i)+1
              endif
            elseif ( i2.eq.mtj(j,4) ) then
              if ( i3.eq.mtj(j,2) ) then
                itet(j) = 1
                itri(i) = itri(i)+1
              elseif ( i3.eq.mtj(j,3) ) then
                itet(j) = -1
                itri(i) = itri(i)+1
              endif
            endif
          endif
          if ( i1.eq.mtj(j,2) ) then
            if ( i2.eq.mtj(j,1) ) then
              if ( i3.eq.mtj(j,4) ) then
                itet(j) = 1
                itri(i) = itri(i)+1
              elseif ( i3.eq.mtj(j,3) ) then
                itet(j) = -1
                itri(i) = itri(i)+1
              endif
            elseif ( i2.eq.mtj(j,3) ) then
              if ( i3.eq.mtj(j,1) ) then
                itet(j) = 1
                itri(i) = itri(i)+1
              elseif ( i3.eq.mtj(j,4) ) then
                itet(j) = -1
                itri(i) = itri(i)+1
              endif
            elseif ( i2.eq.mtj(j,4) ) then
              if ( i3.eq.mtj(j,1) ) then
                itet(j) = -1
                itri(i) = itri(i)+1
              elseif ( i3.eq.mtj(j,3) ) then
                itet(j) = 1
                itri(i) = itri(i)+1
              endif
            endif
          endif
          if ( i1.eq.mtj(j,3) ) then
            if ( i2.eq.mtj(j,1) ) then
              if ( i3.eq.mtj(j,2) ) then
                itet(j) = 1
                itri(i) = itri(i)+1
              elseif ( i3.eq.mtj(j,4) ) then
                itet(j) = -1
                itri(i) = itri(i)+1
              endif
            elseif ( i2.eq.mtj(j,2) ) then
              if ( i3.eq.mtj(j,1) ) then
                itet(j) = -1
                itri(i) = itri(i)+1
              elseif ( i3.eq.mtj(j,4) ) then
                itet(j) = 1
                itri(i) = itri(i)+1
              endif
            elseif ( i2.eq.mtj(j,4) ) then
              if ( i3.eq.mtj(j,1) ) then
                itet(j) = 1
                itri(i) = itri(i)+1
              elseif ( i3.eq.mtj(j,2) ) then
                itet(j) = -1
                itri(i) = itri(i)+1
              endif
            endif
          endif
          if ( i1.eq.mtj(j,4) ) then
            if ( i2.eq.mtj(j,1) ) then
              if ( i3.eq.mtj(j,2) ) then
                itet(j) = -1
                itri(i) = itri(i)+1
              elseif ( i3.eq.mtj(j,3) ) then
                itet(j) = 1
                itri(i) = itri(i)+1
              endif
            elseif ( i2.eq.mtj(j,2) ) then
              if ( i3.eq.mtj(j,1) ) then
                itet(j) = 1
                itri(i) = itri(i)+1
              elseif ( i3.eq.mtj(j,3) ) then
                itet(j) = -1
                itri(i) = itri(i)+1
              endif
            elseif ( i2.eq.mtj(j,3) ) then
              if ( i3.eq.mtj(j,1) ) then
                itet(j) = -1
                itri(i) = itri(i)+1
              elseif ( i3.eq.mtj(j,2) ) then
                itet(j) = 1
                itri(i) = itri(i)+1
              endif
            endif
          endif

  110   continue

  100 continue


      return
      end
c==
c**********   subroutine coinedge   **********************************

c     purpose : pick up the tetrahedron which coincides 
c                                       with the surface triangle edge
c     last modified : 2006.01.01

      subroutine coinedge(kte,nelm,mtj,x,y,z,
     &                    netri,mtri,xtri,ytri,ztri,itet,itri)
      implicit real*8(a-h,o-z)

      dimension mtj(kte,4),x(*),y(*),z(*),
     &          mtri(kte,3),xtri(*),ytri(*),ztri(*),
     &          itet(*),itri(*)

      dimension jp(4)

      data tole/1.d-3/,actol/1.d-10/,dth/1.745329252d-3/
c                                           --> 0.1 deg

      do 100 i=1,netri

        if ( itri(i).eq.0 ) then

          call trieq(xtri(mtri(i,1)),ytri(mtri(i,1)),ztri(mtri(i,1)),
     &               xtri(mtri(i,2)),ytri(mtri(i,2)),ztri(mtri(i,2)),
     &               xtri(mtri(i,3)),ytri(mtri(i,3)),ztri(mtri(i,3)),
     &               a0,b0,c0,d0)

          do 110 ie=1,3

            i1 = mtri(i,ie)
            i2 = mtri(i,mod(ie,3)+1)

            do 200 j=1,nelm

              jp(1) = 0
              jp(2) = 0
              jp(3) = 0
              jp(4) = 0

              do 210 jn=1,4
                if ( i1.eq.mtj(j,jn) ) jp(1) = jn
                if ( i2.eq.mtj(j,jn) ) jp(2) = jn
  210         continue

              if ( jp(1).eq.0 .or. jp(2).eq.0 ) goto 200

              if ( jp(1)+jp(2).eq.3 ) then
                jp(3) = 3
                jp(4) = 4
              elseif ( jp(1)+jp(2).eq.4 ) then
                jp(3) = 2
                jp(4) = 4
              elseif ( jp(1)+jp(2).eq.6 ) then
                jp(3) = 1
                jp(4) = 3
              elseif ( jp(1)+jp(2).eq.7 ) then
                jp(3) = 1
                jp(4) = 2
              else
                if ( iabs(jp(1)-jp(2)).eq.3 ) then
                  jp(3) = 2
                  jp(4) = 3
                else
                  jp(3) = 1
                  jp(4) = 4
                endif
              endif

              delh = hmin(x(mtj(j,1)),y(mtj(j,1)),z(mtj(j,1)),
     &                    x(mtj(j,2)),y(mtj(j,2)),z(mtj(j,2)),
     &                    x(mtj(j,3)),y(mtj(j,3)),z(mtj(j,3)),
     &                    x(mtj(j,4)),y(mtj(j,4)),z(mtj(j,4)))

              if ( dabs(a0*x(mtj(j,jp(3)))+b0*y(mtj(j,jp(3)))
     &                 +c0*z(mtj(j,jp(3)))+d0).le.tole*delh ) then
                j3 = 3
                j4 = 4
              elseif ( dabs(a0*x(mtj(j,jp(4)))+b0*y(mtj(j,jp(4)))
     &                +c0*z(mtj(j,jp(4)))+d0).le.tole*delh ) then
                j3 = 4
                j4 = 3
              else
                goto 200
              endif

              call trieq(x(mtj(j,jp(1))),y(mtj(j,jp(1))),
     &                   z(mtj(j,jp(1))),
     &                   x(mtj(j,jp(2))),y(mtj(j,jp(2))),
     &                   z(mtj(j,jp(2))),
     &                   x(mtj(j,jp(j3))),y(mtj(j,jp(j3))),
     &                   z(mtj(j,jp(j3))),
     &                   a,b,c,d)

              pr = a0*a+b0*b+c0*c
              u0 = dsqrt(a0*a0+b0*b0+c0*c0)
              v0 = dsqrt(a*a+b*b+c*c)
              cs = pr/(u0*v0)
              if ( dabs(cs).le.1.d0 ) then
                th = dacos(cs)
              elseif ( dabs(cs).le.1.d0+actol ) then
                th = dacos(1.d0)
              else
                write(6,600) i,j,cs
                stop 900
              endif

              if ( th.lt.dth ) then

                if ( itet(j).eq.0 ) then

                  if ( a0*x(mtj(j,jp(j4)))+b0*y(mtj(j,jp(j4)))
     &                +c0*z(mtj(j,jp(j4)))+d0.gt.0.d0 ) then
                    itet(j) =  1
                  else
                    itet(j) = -1
                  endif

                endif

                itri(i) = itri(i)+1

              endif

  200       continue
  110     continue
        endif
  100 continue


      return
  600 format(' *** ERROR in sub.coinedge. The cosine is more than 1.0.'/
     &       '  Tri-ID = ',i6,', Tet-ID = ',i6,', Cos = ',1pe20.13)
      end
c==
c**********   subroutine coinnode   **********************************

c     purpose : pick up the tetrahedron which coincides 
c                                       with the surface triangle node
c     last modified : 2006.01.01

      subroutine coinnode(kte,nelm,mtj,x,y,z,
     &                    netri,mtri,xtri,ytri,ztri,itet,itri)
      implicit real*8(a-h,o-z)

      dimension mtj(kte,4),x(*),y(*),z(*),
     &          mtri(kte,3),xtri(*),ytri(*),ztri(*),
     &          itet(*),itri(*)

      dimension jp(9,4),vec(3,3),jpl(2)

      data tole/1.d-3/,actol/1.d-10/,
     &     jp/2,3,4,2,4,3,3,4,2, 1,3,4,1,4,3,3,4,1,
     &        1,2,4,1,4,2,2,4,1, 1,2,3,1,3,2,2,3,1/


      do 100 i=1,netri

        if ( itri(i).eq.0 ) then

          call trieq(xtri(mtri(i,1)),ytri(mtri(i,1)),ztri(mtri(i,1)),
     &               xtri(mtri(i,2)),ytri(mtri(i,2)),ztri(mtri(i,2)),
     &               xtri(mtri(i,3)),ytri(mtri(i,3)),ztri(mtri(i,3)),
     &               a0,b0,c0,d0)

          do 110 in=1,3

            ip0 = mtri(i,in)

            do 200 j=1,nelm

              if ( itet(j).eq.0 ) then

                do 210 jn=1,4
                  if ( ip0.eq.mtj(j,jn) ) then
                    jn0 = jn
                    goto 220
                  endif
  210           continue
                goto 200

  220           continue
                delh = hmin(x(mtj(j,1)),y(mtj(j,1)),z(mtj(j,1)),
     &                      x(mtj(j,2)),y(mtj(j,2)),z(mtj(j,2)),
     &                      x(mtj(j,3)),y(mtj(j,3)),z(mtj(j,3)),
     &                      x(mtj(j,4)),y(mtj(j,4)),z(mtj(j,4)))

                do 300 k=1,7,3

                  jpl(1) = jp(k  ,jn0)
                  jpl(2) = jp(k+1,jn0)

                  if ( (dabs(a0*x(mtj(j,jpl(1)))+b0*y(mtj(j,jpl(1)))
     &                   +c0*z(mtj(j,jpl(1)))+d0).le.tole*delh)
     &                 .and.
     &                 (dabs(a0*x(mtj(j,jpl(2)))+b0*y(mtj(j,jpl(2)))
     &                   +c0*z(mtj(j,jpl(2)))+d0).le.tole*delh) ) then
                    k0 = k
                    goto 310
                  endif

  300           continue
                goto 200

  310           continue
                ip1 = mtri(i,mod(in,3)+1)
                ip2 = mtri(i,mod(in+1,3)+1)

                vec(1,1) = xtri(ip1)-xtri(ip0)
                vec(2,1) = ytri(ip1)-ytri(ip0)
                vec(3,1) = ztri(ip1)-ztri(ip0)
                vec(1,2) = xtri(ip2)-xtri(ip0)
                vec(2,2) = ytri(ip2)-ytri(ip0)
                vec(3,2) = ztri(ip2)-ztri(ip0)
                pr = vec(1,1)*vec(1,2)+vec(2,1)*vec(2,2)
     &               +vec(3,1)*vec(3,2)
                vlen1 = dsqrt(vec(1,1)*vec(1,1)+vec(2,1)*vec(2,1)+
     &                        vec(3,1)*vec(3,1))
                vlen2 = dsqrt(vec(1,2)*vec(1,2)+vec(2,2)*vec(2,2)+
     &                        vec(3,2)*vec(3,2))
                cs = pr/(vlen1*vlen2)
                if ( dabs(cs).le.1.d0 ) then
                  th12 = dacos(cs)
                elseif ( dabs(cs).le.1.d0+actol ) then
                  th12 = dacos(1.d0)
                else
                  write(6,600) i,j,cs
                  stop 900
                endif

                do 400 kp=1,2
                  vec(1,3) = x(mtj(j,jpl(kp)))-xtri(ip0)
                  vec(2,3) = y(mtj(j,jpl(kp)))-ytri(ip0)
                  vec(3,3) = z(mtj(j,jpl(kp)))-ztri(ip0)
                  pr = vec(1,1)*vec(1,3)+vec(2,1)*vec(2,3)
     &                 +vec(3,1)*vec(3,3)
                  vlen3 = dsqrt(vec(1,3)*vec(1,3)+vec(2,3)*vec(2,3)+
     &                          vec(3,3)*vec(3,3))
                  cs = pr/(vlen1*vlen3)
                  if ( dabs(cs).le.1.d0 ) then
                    th13 = dacos(cs)
                  elseif ( dabs(cs).le.1.d0+actol ) then
                    th13 = dacos(1.d0)
                  else
                    write(6,600) i,j,cs
                    stop 900
                  endif
                  pr = vec(1,2)*vec(1,3)+vec(2,2)*vec(2,3)
     &                 +vec(3,2)*vec(3,3)
                  cs = pr/(vlen2*vlen3)
                  if ( dabs(cs).le.1.d0 ) then
                    th23 = dacos(cs)
                  elseif ( dabs(cs).le.1.d0+actol ) then
                    th23 = dacos(1.d0)
                  else
                    write(6,600) i,j,cs
                    stop 900
                  endif

                  if ( dabs(th12-th13-th23).le.actol ) goto 410

  400           continue
                goto 200

  410           continue

                if ( a0*x(mtj(j,jp(k0+2,jn0)))
     &              +b0*y(mtj(j,jp(k0+2,jn0)))
     &              +c0*z(mtj(j,jp(k0+2,jn0)))+d0.gt.0.d0 ) then
                  itet(j) =  1
                else
                  itet(j) = -1
                endif

                itri(i) = itri(i)+1

              endif
  200       continue
  110     continue
        endif
  100 continue


      return
  600 format(' *** ERROR in sub.coinnode. The cosine is more than 1.0.'/
     &       '  Tri-ID = ',i6,', Tet-ID = ',i6,', Cos = ',1pe20.13)
      end
c==
c**********   subroutine trieq      *********************************

c     purpose : calculate a surface equation
c     last modified : 2006.01.01

      subroutine trieq(xi,yi,zi,xj,yj,zj,xk,yk,zk,a,b,c,d)
      implicit real*8(a-h,o-z)


      a = yi*zj+yj*zk+yk*zi-yi*zk-yj*zi-yk*zj
      b = zi*xj+zj*xk+zk*xi-zi*xk-zj*xi-zk*xj
      c = xi*yj+xj*yk+xk*yi-xi*yk-xj*yi-xk*yj
      d = -a*xi-b*yi-c*zi


      return
      end
c==
c**********   function hmin        *********************************

c     purpose : calculate the minimum hight with the tetrahedron
c     last modified : 2006.01.01


      function hmin(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4)
      implicit real*8(a-h,o-z)


      call trieq(x1,y1,z1,x2,y2,z2,x3,y3,z3,a,b,c,d)
      hmin = dabs(a*x4+b*y4+c*z4+d)/dsqrt(a*a+b*b+c*c)

      call trieq(x1,y1,z1,x2,y2,z2,x4,y4,z4,a,b,c,d)
      hmin = dmin1(hmin,dabs(a*x3+b*y3+c*z3+d)/dsqrt(a*a+b*b+c*c))

      call trieq(x1,y1,z1,x3,y3,z3,x4,y4,z4,a,b,c,d)
      hmin = dmin1(hmin,dabs(a*x2+b*y2+c*z2+d)/dsqrt(a*a+b*b+c*c))

      call trieq(x2,y2,z2,x3,y3,z3,x4,y4,z4,a,b,c,d)
      hmin = dmin1(hmin,dabs(a*x1+b*y1+c*z1+d)/dsqrt(a*a+b*b+c*c))


      return
      end
c==
c**********   subroutine data   **************************************

c     purpose : print results on data file
c     last modified : 2006.01.01

      subroutine data(kte,nelm,node,mtj,jac,x,y,z,itet,new,neid)
      implicit real*8(a-h,o-z)

      dimension mtj(kte,4),jac(kte,4),x(*),y(*),z(*),itet(*),neid(*)

      dimension jp(4)

      character fname*30


      write(6,*)
      write(6,600) ' Output File Name = ?'
      read(5,600) fname

      open(9,file=fname)

      write(9,600) 'nelm,node'
      write(9,610) new,node

      write(9,600) 'e-id,mtj,jac'
      do 10 i=1,nelm
        if ( itet(i).eq.-1 ) then
          do 15 j=1,4
            if ( jac(i,j).eq.0 ) then
              jp(j) = 0
            else
              jp(j) = neid(jac(i,j))
            endif
   15     continue
          write(9,620) neid(i),(mtj(i,j),j=1,4),(jp(j),j=1,4)
        endif
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
