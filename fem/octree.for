c     ------   program octree  --------------------------------------
c     purpose : generate new nodes using octree method
c     last modified : 2006.01.01

      implicit real*8(a-h,o-z)
      parameter(ktj=100,kte=500,ktb=100,
     &          ktc=50,ktp=100,ktl=6)


      dimension mtj(kte,4),jac(kte,4),x(ktj),y(ktj),z(ktj),
     &          naiten(ktc),xx(ktp),yy(ktp),zz(ktp),
     &          isankak(ktb,3),hentyou(ktb),icube(ktc,8),ibunk(ktc),
     &          level(ktl),kaisou(ktl,ktc),heimen(ktb,4),istack(ktp)

      common /unit/ku90,ku91


      ku90 = 90
      ku91 = 91

c  -- data input

      call input(kte,ktj,nelm,node,mtj,jac,x,y,z)

c  -- apply octree method

      call oct(kte,ktj,ktb,ktc,ktp,ktl,nelm,node,mtj,jac,
     &         x,y,z,ncube,naiten,num,xx,yy,zz,inode,
     &         isankak,hentyou,icube,ibunk,level,kaisou,
     &         heimen,istack)

c  -- hexhedra for octree data output

      call data1(ktp,ktc,ncube,num,xx,yy,zz)

c  -- original and new nodes output

      call data2(ktj,ktp,node,x,y,z,inode)


      stop
      end
c==
c**********   subroutine input   ************************************

c     purpose : data input from 3delaun
c     last modified : 2006.01.01

      subroutine input(kte,ktj,nelm,node,mtj,jac,x,y,z)
      implicit real*8(a-h,o-z)

      dimension mtj(kte,4),jac(kte,4),x(*),y(*),z(*)
      character fname*30


c  -- input file name which was constructed 
c                               by the 3delaun program

      write(*,500)
  500 format(' Input file name = ? ')
      read(*,*) fname
      open(8,file=fname)

c  -- read input data

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


      return
  600 format(a)
      end
c==
c**********   subtoutine oct   **************************************

c     purpose : generate new nodes using octree method
c     last modified : 2006.01.01

      subroutine oct(kte,ktj,ktb,ktc,ktp,ktl,nelm,node,mtj,jac,
     &               x,y,z,ncube,naiten,num,xx,yy,zz,inode,
     &               isankak,hentyou,icube,ibunk,level,kaisou,
     &               heimen,istack)
      implicit real*8(a-h,o-z)

      dimension mtj(kte,4),jac(kte,4),x(*),y(*),z(*),
     &          xx(*),yy(*),zz(*),isankak(ktb,3),hentyou(*),
     &          icube(ktc,8),ibunk(*),level(*),kaisou(ktl,*),
     &          heimen(ktb,4),naiten(*),istack(*)


      call fpoint(kte,ktb,nelm,node,mtj,jac,x,y,z,
     &            a,b,c,itri,heimen,twig,heikin,
     &            isankak,hentyou)

      call initcube(ktc,a,b,c,twig,num,xx,yy,zz,inel,icube)

      ib = 0
      ibunk(1) = 0
      level(1) = 1
      kaisou(1,1) = 1

      waru = twig

  100 continue

      iwar = 1
      do 30 i=1,ib
        iwar = iwar*2
   30 continue
      ib = ib+1
      twig = waru/iwar

      level(ib+1) = 0
      do 10 i=1,level(ib)
        ic = kaisou(ib,i)
        call bunkatu(ktp,ktc,ktl,ib,ic,inel,icube,num,
     &               xx,yy,zz,ibunk,level,kaisou,twig,node,
     &               x,y,z,naiten)
   10 continue

c  -- check the number of times for branch

      if ( ib.eq.ktl-1 ) then
        do 20 i=1,level(ib+1)
          ibunk(kaisou(ib+1,i))=1
   20   continue
        goto 200
      endif

c  -- check the length of branch

      if ( twig.lt.heikin/6.d0 ) then
        do 40 i=1,level(ib+1)
          ibunk(kaisou(ib+1,i)) = 1
   40   continue
        goto 200
      endif

c  -- check the number of cubes include plural nodes

      if ( level(ib+1).eq.0 ) goto 200

      goto 100

  200 continue

      call cube(ktj,ktp,ktc,inel,icube,num,xx,yy,zz,ibunk,
     &          ncube,naiten)

      call remove(kte,ktb,nelm,node,mtj,jac,x,y,z,
     &            num,xx,yy,zz,inode,itri,heimen,
     &            heikin,istack)


      return
      end
c==
c**********   subroutuine fpoint   **********************************

c     purpose : calculate the length of branch
c                                         and a plane equation
c     last modified : 2006.01.01

      subroutine fpoint(kte,ktb,nelm,node,mtj,jac,x,y,z,
     &                  a,b,c,itri,heimen,twig,heikin,
     &                  isankak,hentyou)
      implicit real*8(a-h,o-z)

      dimension mtj(kte,4),jac(kte,4),x(*),y(*),z(*),
     &          isankak(ktb,3),hentyou(*),heimen(ktb,4)


      a = 0.d0
      b = 0.d0
      c = 0.d0
      twig = 0.d0

      sa = 0.d0
      sb = 0.d0
      sc = 0.d0

      do 10 i=1,node
        sa = sa+x(i)
        sb = sb+y(i)
        sc = sc+z(i)
   10 continue

      a = sa/node
      b = sb/node
      c = sc/node

      do 20 i=1,node

        hx = x(i)-a
        hy = y(i)-b
        hz = z(i)-c

        hlong = hx*hx+hy*hy+hz*hz

        if ( sqrt(hlong).gt.twig ) then
          twig = sqrt(hlong)
        endif

   20 continue

      itri = 0

      do 30 i=1,nelm
        do 40 j=1,4
          if ( jac(i,j).eq.0 ) then
            itri = itri+1
            if ( itri.gt.ktb ) then
              write(6,600) ' *** ERROR(sub.fpoint) itri > ktb'
              stop 900
            endif
            isankak(itri,1) = mtj(i,mod(j,4)+1)
            isankak(itri,2) = mtj(i,4-(j-1)/2*2)
            isankak(itri,3) = mtj(i,3-mod(j/2,2)*2)
          endif
   40   continue
   30 continue

      sgoal = 0.d0

      do 50 i=1,itri

        xk1  = x(isankak(i,1))-x(isankak(i,2))
        yk1  = y(isankak(i,1))-y(isankak(i,2))
        zk1  = z(isankak(i,1))-z(isankak(i,2))
        hen1 = sqrt(xk1*xk1+yk1*yk1+zk1*zk1)

        xk2  = x(isankak(i,2))-x(isankak(i,3))
        yk2  = y(isankak(i,2))-y(isankak(i,3))
        zk2  = z(isankak(i,2))-z(isankak(i,3))
        hen2 = sqrt(xk2*xk2+yk2*yk2+zk2*zk2)

        xk3  = x(isankak(i,3))-x(isankak(i,1))
        yk3  = y(isankak(i,3))-y(isankak(i,1))
        zk3  = z(isankak(i,3))-z(isankak(i,1))
        hen3 = sqrt(xk3*xk3+yk3*yk3+zk3*zk3)

        hentyou(i) = hen1+hen2+hen3

   50 continue

      heikin = 0.d0

      do 60 i=1,itri

        x1 = x(isankak(i,1))
        y1 = y(isankak(i,1))
        z1 = z(isankak(i,1))
        x2 = x(isankak(i,2))
        y2 = y(isankak(i,2))
        z2 = z(isankak(i,2))
        x3 = x(isankak(i,3))
        y3 = y(isankak(i,3))
        z3 = z(isankak(i,3))

        aa = y1*z2+y2*z3+y3*z1-y1*z3-y2*z1-y3*z2
        bb = z1*x2+z2*x3+z3*x1-z1*x3-z2*x1-z3*x2
        cc = x1*y2+x2*y3+x3*y1-x1*y3-x2*y1-x3*y2
        dd = -(aa*x1+bb*y1+cc*z1)

        tyou = hentyou(i)
        heikin = heikin+tyou
        heimen(i,1) = aa
        heimen(i,2) = bb
        heimen(i,3) = cc
        heimen(i,4) = dd

   60 continue

      heikin = heikin/itri


      return
  600 format(a)
      end
c==
c**********   subroutuine initcube  **********************************

c     purpose : set an initial cube include the region
c     last modified : 2006.01.01

      subroutine initcube(ktc,a,b,c,twig,num,xx,yy,zz,
     &                    inel,icube)
      implicit real*8(a-h,o-z)

      dimension xx(*),yy(*),zz(*),icube(ktc,8)


      xx(1) = a+twig
      yy(1) = b+twig
      zz(1) = c+twig

      xx(2) = a+twig
      yy(2) = b-twig
      zz(2) = c+twig

      xx(3) = a-twig
      yy(3) = b-twig
      zz(3) = c+twig

      xx(4) = a-twig
      yy(4) = b+twig
      zz(4) = c+twig

      xx(5) = a+twig
      yy(5) = b+twig
      zz(5) = c-twig

      xx(6) = a+twig
      yy(6) = b-twig
      zz(6) = c-twig

      xx(7) = a-twig
      yy(7) = b-twig
      zz(7) = c-twig

      xx(8) = a-twig
      yy(8) = b+twig
      zz(8) = c-twig

      num  = 8
      inel = 1

      do 10 i=1,8
        icube(1,i) = i
   10 continue


      return
      end
c==
c**********   subroutine bunkatu   **********************************

c     purpose : branch out
c     last modified : 2006.01.01

      subroutine bunkatu(ktp,ktc,ktl,ib,ic,inel,icube,num,
     &                   xx,yy,zz,ibunk,level,kaisou,twig,node,
     &                   x,y,z,naiten)
      implicit real*8(a-h,o-z)
      parameter(err=1.0d-4)

      dimension x(*),y(*),z(*),xx(*),yy(*),zz(*),
     &          icube(ktc,8),ibunk(*),
     &          level(*),kaisou(ktl,*),naiten(*)

      dimension xoct(19),yoct(19),zoct(19),noct(19)


      sa = 0.d0
      sb = 0.d0
      sc = 0.d0

      do 10 i=1,8
        sa = sa+xx(icube(ic,i))
        sb = sb+yy(icube(ic,i))
        sc = sc+zz(icube(ic,i))
   10 continue

      a = sa/8
      b = sb/8
      c = sc/8

      xoct(1) = a
      yoct(1) = b
      zoct(1) = c

      xoct(2) = a
      yoct(2) = b
      zoct(2) = c-twig

      xoct(3) = a
      yoct(3) = b+twig
      zoct(3) = c

      xoct(4) = a+twig
      yoct(4) = b
      zoct(4) = c

      xoct(5) = a
      yoct(5) = b-twig
      zoct(5) = c

      xoct(6) = a-twig
      yoct(6) = b
      zoct(6) = c

      xoct(7) = a
      yoct(7) = b
      zoct(7) = c+twig

      xoct(8) = a
      yoct(8) = b+twig
      zoct(8) = c-twig

      xoct(9) = a+twig
      yoct(9) = b
      zoct(9) = c-twig

      xoct(10) = a
      yoct(10) = b-twig
      zoct(10) = c-twig

      xoct(11) = a-twig
      yoct(11) = b
      zoct(11) = c-twig

      xoct(12) = a-twig
      yoct(12) = b+twig
      zoct(12) = c

      xoct(13) = a+twig
      yoct(13) = b+twig
      zoct(13) = c

      xoct(14) = a+twig
      yoct(14) = b-twig
      zoct(14) = c

      xoct(15) = a-twig
      yoct(15) = b-twig
      zoct(15) = c

      xoct(16) = a
      yoct(16) = b+twig
      zoct(16) = c+twig

      xoct(17) = a+twig
      yoct(17) = b
      zoct(17) = c+twig

      xoct(18) = a
      yoct(18) = b-twig
      zoct(18) = c+twig

      xoct(19) = a-twig
      yoct(19) = b
      zoct(19) = c+twig

      do 20 i=1,19
        do 30 j=1,num
          if ( dabs(xx(j)-xoct(i)).le.err ) then
            if ( dabs(yy(j)-yoct(i)).le.err ) then
              if ( dabs(zz(j)-zoct(i)).le.err ) then
                noct(i) = j
                goto 20
              endif
            endif
          endif
   30   continue
        num = num+1
        if ( num.gt.ktp ) then
          write(6,600) ' *** ERROR(sub.bunkatu) num > ktp'
          stop 900
        endif
        noct(i) = num
        xx(num) = xoct(i)
        yy(num) = yoct(i)
        zz(num) = zoct(i)
   20 continue

c  -- cube 1

      inel = inel+1
      if ( inel.gt.ktc ) then
        write(6,600) ' *** ERROR(sub.bunkatu) inel > ktc'
        stop 900
      endif
      icube(inel,1) = icube(ic,1)
      icube(inel,2) = noct(17)
      icube(inel,3) = noct(7)
      icube(inel,4) = noct(16)
      icube(inel,5) = noct(13)
      icube(inel,6) = noct(4)
      icube(inel,7) = noct(1)
      icube(inel,8) = noct(3)
      ibunk(inel) = 1

      call setten(ktc,node,x,y,z,inel,icube,xx,yy,zz,iq,naiten)

      if ( iq.gt.1 ) then
        level(ib+1) = level(ib+1)+1
        kaisou(ib+1,level(ib+1)) = inel
        ibunk(inel)  = 0
        naiten(inel) = -iq
      endif

c  -- cube 2

      inel = inel+1
      if ( inel.gt.ktc ) then
        write(6,600) ' *** ERROR(sub.bunkatu) inel > ktc'
        stop 900
      endif
      icube(inel,1) = noct(17)
      icube(inel,2) = icube(ic,2)
      icube(inel,3) = noct(18)
      icube(inel,4) = noct(7)
      icube(inel,5) = noct(4)
      icube(inel,6) = noct(14)
      icube(inel,7) = noct(5)
      icube(inel,8) = noct(1)
      ibunk(inel) = 1

      call setten(ktc,node,x,y,z,inel,icube,xx,yy,zz,iq,naiten)

      if ( iq.gt.1 ) then
        level(ib+1) = level(ib+1)+1
        kaisou(ib+1,level(ib+1)) = inel
        ibunk(inel)  = 0
        naiten(inel) = -iq
      endif

c  -- cube 3

      inel = inel+1
      if ( inel.gt.ktc ) then
        write(6,600) ' *** ERROR(sub.bunkatu) inel > ktc'
        stop 900
      endif
      icube(inel,1) = noct(7)
      icube(inel,2) = noct(18)
      icube(inel,3) = icube(ic,3)
      icube(inel,4) = noct(19)
      icube(inel,5) = noct(1)
      icube(inel,6) = noct(5)
      icube(inel,7) = noct(15)
      icube(inel,8) = noct(6)
      ibunk(inel) = 1

      call setten(ktc,node,x,y,z,inel,icube,xx,yy,zz,iq,naiten)

      if ( iq.gt.1 ) then
        level(ib+1) = level(ib+1)+1
        kaisou(ib+1,level(ib+1)) = inel
        ibunk(inel)  = 0
        naiten(inel) = -iq
      endif

c  -- cube 4

      inel = inel+1
      if ( inel.gt.ktc ) then
        write(6,600) ' *** ERROR(sub.bunkatu) inel > ktc'
        stop 900
      endif
      icube(inel,1) = noct(16)
      icube(inel,2) = noct(7)
      icube(inel,3) = noct(19)
      icube(inel,4) = icube(ic,4)
      icube(inel,5) = noct(3)
      icube(inel,6) = noct(1)
      icube(inel,7) = noct(6)
      icube(inel,8) = noct(12)
      ibunk(inel) = 1

      call setten(ktc,node,x,y,z,inel,icube,xx,yy,zz,iq,naiten)

      if ( iq.gt.1 ) then
        level(ib+1) = level(ib+1)+1
        kaisou(ib+1,level(ib+1)) = inel
        ibunk(inel)  = 0
        naiten(inel) = -iq
      endif

c  -- cube 5

      inel = inel+1
      if ( inel.gt.ktc ) then
        write(6,600) ' *** ERROR(sub.bunkatu) inel > ktc'
        stop 900
      endif
      icube(inel,1) = noct(13)
      icube(inel,2) = noct(4)
      icube(inel,3) = noct(1)
      icube(inel,4) = noct(3)
      icube(inel,5) = icube(ic,5)
      icube(inel,6) = noct(9)
      icube(inel,7) = noct(2)
      icube(inel,8) = noct(8)
      ibunk(inel) = 1

      call setten(ktc,node,x,y,z,inel,icube,xx,yy,zz,iq,naiten)

      if ( iq.gt.1 ) then
        level(ib+1) = level(ib+1)+1
        kaisou(ib+1,level(ib+1)) = inel
        ibunk(inel)  = 0
        naiten(inel) = -iq
      endif

c  -- cube 6

      inel = inel+1
      if ( inel.gt.ktc ) then
        write(6,600) ' *** ERROR(sub.bunkatu) inel > ktc'
        stop 900
      endif
      icube(inel,1) = noct(4)
      icube(inel,2) = noct(14)
      icube(inel,3) = noct(5)
      icube(inel,4) = noct(1)
      icube(inel,5) = noct(9)
      icube(inel,6) = icube(ic,6)
      icube(inel,7) = noct(10)
      icube(inel,8) = noct(2)
      ibunk(inel) = 1

      call setten(ktc,node,x,y,z,inel,icube,xx,yy,zz,iq,naiten)

      if ( iq.gt.1 ) then
        level(ib+1) = level(ib+1)+1
        kaisou(ib+1,level(ib+1)) = inel
        ibunk(inel)  = 0
        naiten(inel) = -iq
      endif

c  -- cube 7

      inel = inel+1
      if ( inel.gt.ktc ) then
        write(6,600) ' *** ERROR(sub.bunkatu) inel > ktc'
        stop 900
      endif
      icube(inel,1) = noct(1)
      icube(inel,2) = noct(5)
      icube(inel,3) = noct(15)
      icube(inel,4) = noct(6)
      icube(inel,5) = noct(2)
      icube(inel,6) = noct(10)
      icube(inel,7) = icube(ic,7)
      icube(inel,8) = noct(11)
      ibunk(inel) = 1

      call setten(ktc,node,x,y,z,inel,icube,xx,yy,zz,iq,naiten)

      if ( iq.gt.1 ) then
        level(ib+1) = level(ib+1)+1
        kaisou(ib+1,level(ib+1)) = inel
        ibunk(inel)  = 0
        naiten(inel) = -iq
      endif

c  -- cube 8

      inel = inel+1
      if ( inel.gt.ktc ) then
        write(6,600) ' *** ERROR(sub.bunkatu) inel > ktc'
        stop 900
      endif
      icube(inel,1) = noct(3)
      icube(inel,2) = noct(1)
      icube(inel,3) = noct(6)
      icube(inel,4) = noct(12)
      icube(inel,5) = noct(8)
      icube(inel,6) = noct(2)
      icube(inel,7) = noct(11)
      icube(inel,8) = icube(ic,8)
      ibunk(inel) = 1

      call setten(ktc,node,x,y,z,inel,icube,xx,yy,zz,iq,naiten)

      if ( iq.gt.1 ) then
        level(ib+1) = level(ib+1)+1
        kaisou(ib+1,level(ib+1)) = inel
        ibunk(inel)  = 0
        naiten(inel) = -iq
      endif


      return
  600 format(a)
      end
c==
c**********   subroutine setten   ***********************************

c     purpose : search for a node in the cube
c     last modified : 2006.01.01

      subroutine setten(ktc,node,x,y,z,inel,icube,xx,yy,zz,iq,naiten)
      implicit real*8(a-h,o-z)

      dimension x(*),y(*),z(*),xx(*),yy(*),zz(*),
     &          icube(ktc,8),naiten(*)


      iq = 0
      naiten(inel) = 0

      do 10 i=1,node
        if ( x(i).le.xx(icube(inel,1)) ) then
          if ( x(i).ge.xx(icube(inel,7)) ) then
            if ( y(i).le.yy(icube(inel,1)) ) then
              if ( y(i).ge.yy(icube(inel,7)) ) then
                if ( z(i).le.zz(icube(inel,1)) ) then
                  if ( z(i).ge.zz(icube(inel,7)) ) then
                    iq = iq+1
                    naiten(inel) = i
                  endif
                endif
              endif
            endif
          endif
        endif
   10 continue


      return
      end
c==
c**********   subrouotine cube   ***********************************

c     purpose : make all cube
c     last modified : 2006.01.01

      subroutine cube(ktj,ktp,ktc,inel,icube,num,xx,yy,zz,ibunk,
     &                ncube,naiten)
      implicit real*8(a-h,o-z)

      dimension icube(ktc,8),ibunk(*),naiten(*)

      common /unit/ku90,ku91


      open(ku90,form='unformatted',status='scratch')

      ncube = 0
      do 10 i=1,inel
        if ( ibunk(i).ne.0 ) then
          ncube = ncube+1
          write(ku90) ncube,(icube(i,j),j=1,8),naiten(i)
        endif
   10 continue


      return
      end
c==
c**********   subroutine remove   ***********************************

c     purpose : remove points of outside region or unnecessary
c     last modified : 2006.01.01

      subroutine remove(kte,ktb,nelm,node,mtj,jac,x,y,z,
     &                  num,xx,yy,zz,inode,itri,heimen,
     &                  heikin,istack)
      implicit real*8(a-h,o-z)

      dimension mtj(kte,4),jac(kte,4),x(*),y(*),z(*),
     &          xx(*),yy(*),zz(*),
     &          istack(*),heimen(ktb,4)

      common /unit/ku90,ku91

      data coef1/0.25d0/,coef2/0.25d0/

c  -- dist1; tolerance of the distance
c                         from a cube node to an original node
c  -- dist2; tolerance of the distance
c                         from a cube node to a boundary face
c  -- the user must select the suitable value for dist1 and dist2


      dist1 = coef1*heikin
      dist2 = coef2*heikin

      do 10 i=1,num

        amin = 1.d12
        istack(i) = 1
        xi = xx(i)
        yi = yy(i)
        zi = zz(i)

        do 20 j=1,node
          xo  = x(j)-xi
          yo  = y(j)-yi
          zo  = z(j)-zi
          ans = sqrt(xo*xo+yo*yo+zo*zo)
          if ( ans.lt.amin ) then
            kinten = j
            amin = ans
          endif
   20   continue

c  -- check the distance from a cube node to an original node

        if ( amin.lt.dist1 ) then
          istack(i) = 0
          goto 10
        endif

        do 30 j=1,nelm
          do 40 k=1,4
            if ( mtj(j,k).eq.kinten ) then
              jtet = j
              goto 100
            endif
   40     continue
   30   continue

  100   continue

        call location(kte,xi,yi,zi,x,y,z,nelm,mtj,jac,jtet)

c  -- check a cube node which is outside the boundary surface

        if ( jtet.eq.0 ) then
          istack(i) = 0
          goto 10
        endif

        do 60 j=1,itri

          a1 = heimen(j,1)*xi
          b1 = heimen(j,2)*yi
          c1 = heimen(j,3)*zi
          d1 = heimen(j,4)
          away = dabs(a1+b1+c1+d1)
          away = away/sqrt(heimen(j,1)*heimen(j,1)
     &                    +heimen(j,2)*heimen(j,2)
     &                    +heimen(j,3)*heimen(j,3))

c  -- check the distance from a cube node to a boundary face

          if ( away.lt.dist2 ) then
            istack(i) = 0
            goto 10
          endif

   60   continue

   10 continue

      open(ku91,form='unformatted',status='scratch')

      inode = 0
      do 50 i=1,num
        if ( istack(i).ne.0 ) then
          inode = inode+1
          write(ku91) inode,xx(i),yy(i),zz(i)
        endif
   50 continue


      return
      end
c==
c**********   subroutine location   *********************************

c     purpose : search for the element include additional node
c     last modified : 2006.01.01

      subroutine location(kte,xi,yi,zi,x,y,z,nelm,mtj,jac,jtet)
      implicit real*8(a-h,o-z)
      parameter(err=1.0d-12)

      dimension mtj(kte,4),jac(kte,4),x(*),y(*),z(*)


c  -- jtet; the element include additional node

   10 continue

      do 20 n=1,4
        i = mtj(jtet,mod(n,4)+1)
        j = mtj(jtet,4-(n-1)/2*2)
        k = mtj(jtet,3-mod(n/2,2)*2)
        a = y(i)*z(j)+y(j)*z(k)+y(k)*z(i)-y(i)*z(k)-y(j)*z(i)-y(k)*z(j)
        b = z(i)*x(j)+z(j)*x(k)+z(k)*x(i)-z(i)*x(k)-z(j)*x(i)-z(k)*x(j)
        c = x(i)*y(j)+x(j)*y(k)+x(k)*y(i)-x(i)*y(k)-x(j)*y(i)-x(k)*y(j)
        d = -a*x(i)-b*y(i)-c*z(i)
        if ( a*xi+b*yi+c*zi+d.lt.-err ) then
          jtet = jac(jtet,n)
          if ( jtet.eq.0 ) goto 100
          goto 10
        endif
   20 continue

  100 continue


      return
      end
c==
c**********   subroutine data1   *************************************

c     purpose : hexhedra for octree data output
c     last modified : 2006.01.01

      subroutine data1(ktp,ktc,ncube,num,xx,yy,zz)
      implicit real*8(a-h,o-z)

      dimension xx(*),yy(*),zz(*)

      dimension mthex(8)
      character fname*30

      common /unit/ku90,ku91


      rewind(ku90)

      write(*,500)
      read(*,*) judge

      if ( judge.eq.0 ) then
        write(*,510)
        read(*,110) fname
        open(9,file=fname)
        write(9,600) ncube,num
        do 10 i=1,ncube
          read(ku90) id,(mthex(j),j=1,8),nai
          write(9,610) id,(mthex(j),j=1,8),nai
   10   continue
        do 20 i=1,num
          write(9,620) i,xx(i),yy(i),zz(i)
   20   continue
        close(9)
      end if

      close(ku90)


      return
  500 format(/,' Save cubes for octree? (0-yes,1-no) ')
  510 format(' File name?')
  110 format(a)
  600 format(i7,','i7)
  610 format(i7,9(',',i7))
  620 format(i7,',',1pe15.6,',',e15.6,',',e15.6)
      end
c==
c**********   subroutine data2   ************************************

c     purpose : original and new nodes output
c     last modified : 2006.01.01

      subroutine data2(ktj,ktp,node,x,y,z,inode)
      implicit real*8(a-h,o-z)

      dimension x(*),y(*),z(*)
      character fname*30

      common /unit/ku90,ku91


      rewind(ku91)

      write(*,500)
      read(*,*)judge

      if ( judge.eq.0 ) then
        write(*,510)
        read(*,110) fname
        open(9,file=fname)
        write(9,600) node+inode
        do 10 i=1,node
          write(9,620) x(i),y(i),z(i)
   10   continue
        do 20 i=1,inode
          read(ku91) id,xa,ya,za
          write(9,620) xa,ya,za
   20   continue
        close(9)
      else
        write(*,700)
        write(*,510)
        read(*,110) fname
        open(9,file=fname)
        write(9,600) inode
        do 30 i=1,inode
          read(ku91) id,xa,ya,za
          write(9,620) xa,ya,za
   30   continue
        close(9)
      endif

      close(ku91)


      return
  500 format(/,' Save original and new nodes? (0-yes,1-no) ')
  510 format(' File nameÅH')
  110 format(a)
  600 format(i7)
  620 format(1pe15.6,',',e15.6,',',e15.6)
  700 format(/,' Save only new nodes. ')
      end
