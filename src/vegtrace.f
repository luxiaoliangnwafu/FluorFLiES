!*******************************************************
!     simulate the optical thickness in the canopy
!     along the photon trajectory
!
!     Written by H. Kobayashi
!     Last modified 08/04/02
!*******************************************************

      subroutine vegtrace(tau, x, y, z, uxr, uyr, uzr, sflag)

      implicit none
      include 'common.inc'
      include 'math.inc'

!     sflg: stem flag=0 stem collision, 1= no stem
      integer i, ith, inobj, idiv, l
      integer flag, sflag, ivox, io, tio, io2, io12, unobj
      integer oface, pface

      real x, y, z
      real x0, y0, z0, x1, y1, z1
      real xb, yb, zb, th, rio, rio12, d, cf, cf12
      real zlim(2), uxr, uyr, uzr
      real so, sp, ts, so2, so12, conv, intv(3), mgn, rb12
      real tau, tauc, tauc12, taub, tobj(5),tobjb(5), tobj12(5)
      real facos

!     merginal value
      mgn = 1.d-2

      inobj = 1

!     wieght minimum limit
      conv = 1.d-8

!     initial distance from object
      so = 1.d5

      sflag = 1

!     optical thickness
      tau = 0.0
      tauc12 = 0.0
      tauc = 0.0
      taub = 0.0

      x0 = x
      y0 = y
      z0 = max(z,mgn)


      rb12 = 1.0
c      rb12 = 0.8


      do i = 1, 3
         intv(i) = 50. / res
      end do

!     z boudary for trace flag=1 (zlim=0.0), flag=2 (zlim=zmax)
      zlim(1) = 0.0
      zlim(2) = zmax
      flag = 1.5 + sign(0.5, uzr)

!     do while photon reaches the terminal point
      do

!     determinatin of first input voxel
         x1 = aint(x0 / intv(1))
         y1 = aint(y0 / intv(2))
         z1 = aint(z0 / intv(3))
         ivox = ixmax * iymax * int(z1)
         ivox = ivox + int(y1) * iymax
         ivox = ivox + int(x1) + 1

         x1 = x1 * intv(1)
         y1 = y1 * intv(2)
         z1 = z1 * intv(3)

         io = 1
         !write(*,*) "running_mode = ", running_mode
         !write(*,*) "ndivs(41) =",ndivs(41)
         !write(*,*) "ndivs(50) =",ndivs(50)

!     check the photon intersection with big-voxel walls
         call planes(sp, x0, y0, z0,
     &        uxr, uyr, uzr, x1, y1, z1, pface, intv)
          !write(*,*) "x0,y0,z0=", x0,y0,z0
          !write(*,*) "x1,y1,z1=", x1,y1,z1
          !write(*,*) "uxr, uyr, uzr=", uxr, uyr, uzr
          !write(*,*) "sp=", sp
          !write(*,*) "ivox=", ivox
          !write(*,*) "ndivs(ivox)=", ndivs(ivox)
!     check the photon intersection with objects
         if(ndivs(ivox) .ne. 0)then

            so = 1.d5
            do idiv = 1, ndivs(ivox)

!     selected object number
               i = divs(ivox, idiv)
               !write(*,*) "idiv =",idiv
               !write(*,*) "i =",i
               do l = 1, 5
                  tobj(l) = obj(i,l)
               end do

               if(sobj(i) .eq. 1)then
                  call cones(ts, x0, y0, z0,
     &                 uxr, uyr, uzr, tobj, oface, tio)
               else if((sobj(i) .eq. 2) .or. (sobj(i) .eq. 4))then
                  call cyls(ts, x0, y0, z0,
     &                 uxr, uyr, uzr, tobj, oface, tio)
               else if(sobj(i) .eq. 3) then
                  call elpss(ts ,x0, y0 ,z0 ,
     &                 uxr ,uyr ,uzr ,tobj, tio)
               else if(sobj(i) .eq. 5) then
                  call helpss(ts ,x0, y0 ,z0 ,
     &                 uxr ,uyr ,uzr ,tobj,oface, tio)
               end if

               if(ts .lt. so) then
                  so = ts
                  inobj = i
                  io = tio
               end if

               if(io .eq. 0) exit

            end do

!     if stem collision, return
            if(sobj(inobj) .eq. 4 .and. so .lt. 1.d5) then
               sflag = 0
               return
            end  if
         end if

!     increment of the optical path
         if(io .eq. 0) then
!     check branch optical thickness
            tobjb(1) = tobj(1)
            tobjb(2) = tobj(2)
            tobjb(3) = tobj(3) - tobj(4) * (1. - rb)
     &           * min(1, abs(sobj(inobj) - 5))
            tobjb(4) = tobj(4) * rb
            tobjb(5) = tobj(5) * rb

            tobj12(1) = tobj(1)
            tobj12(2) = tobj(2)
            tobj12(3) = tobj(3) - tobj(4) * (1. - rb12)
     &           * min(1, abs(sobj(inobj) - 5))
            tobj12(4) = tobj(4) * rb12
            tobj12(5) = tobj(5) * rb12

            if(sobj(i) .eq. 1) then
               call cones(so2, x0, y0, z0,
     &              uxr, uyr, uzr, tobjb, oface, io2)
               call cones(so12, x0, y0, z0,
     &              uxr, uyr, uzr, tobj12, oface, io12)
            else if(sobj(i) .eq. 2) then
               call cyls(so2, x0, y0, z0,
     &              uxr, uyr, uzr, tobjb, oface, io2)
               call cyls(so12, x0, y0, z0,
     &              uxr, uyr, uzr, tobj12, oface, io12)
            else if(sobj(i) .eq. 3) then
               call elpss(so2, x0, y0, z0,
     &              uxr ,uyr ,uzr ,tobjb, io2)
               call elpss(so12, x0, y0, z0,
     &              uxr ,uyr ,uzr ,tobj12, io12)
            else if(sobj(i) .eq. 5) then
               call helpss(so2, x0, y0, z0,
     &              uxr ,uyr ,uzr ,tobjb, oface, io2)
               call helpss(so12, x0, y0, z0,
     &              uxr ,uyr ,uzr ,tobj12, oface, io12)
            end if

!     if outside of branch (io12=1) go into the branch
            if(io12 .eq. 1) then
               xb = x0 + (so12 + mgn) * uxr
               yb = y0 + (so12 + mgn) * uyr
               zb = z0 + (so12 + mgn) * uzr

               if(sobj(i) .eq. 1)then
                  call cones(so12, xb, yb, zb,
     &                 uxr, uyr, uzr, tobj12, oface, io12)
              else if(sobj(i) .eq. 2)then
                  call cyls(so12, xb, yb, zb,
     &                 uxr, uyr, uzr, tobj12, oface, io12)
               else if(sobj(i) .eq. 3) then
                  call elpss(so12, xb ,yb ,zb,
     &                 uxr ,uyr ,uzr ,tobj12, io12)
               else if(sobj(i) .eq. 5) then
                  call helpss(so12, xb ,yb ,zb,
     &                 uxr ,uyr ,uzr ,tobj12, oface,io12)
               end if
            end if

!     if outside of branch (io2=1) go into the branch
            if(io2 .eq. 1) then
               xb = x0 + (so2 + mgn) * uxr
               yb = y0 + (so2 + mgn) * uyr
               zb = z0 + (so2 + mgn) * uzr

               if(sobj(i) .eq. 1)then
                  call cones(so2, xb, yb, zb,
     &                 uxr, uyr, uzr, tobjb, oface, io2)
               else if(sobj(i) .eq. 2)then
                  call cyls(so2, xb, yb, zb,
     &                 uxr, uyr, uzr, tobjb, oface, io2)
               else if(sobj(i) .eq. 3) then
                  call elpss(so2, xb ,yb ,zb,
     &                 uxr ,uyr ,uzr ,tobjb, io2)
               else if(sobj(i) .eq. 5) then
                  call helpss(so2, xb ,yb ,zb,
     &                 uxr ,uyr ,uzr ,tobjb, oface,io2)
               end if
            end if

!     calculation of the optical path
            th = facos(uzr)
            ith = int(th * 180. / pi)
            rio = 1. - real(io2)
            rio12 = 1. - real(io12)
            cf = sbar(iobj(i))
            cf12 = sbar(iobj(i))
c            cf12 = 0.25

            taub = rio * (so2 + mgn) * bad(iobj(i)) * gtblb(ith) * bp2
            taub = taub
     &           + rio * (so2 + mgn) * u(iobj(i)) * gtblc(ith)
     &           * 4. * cf * (1. - bp2)

            tauc12 = ((so12 + mgn) - (so2 + mgn) * rio) * rio12
            tauc12 = tauc12
     &           * (u(iobj(i)) * gtblc(ith) * 4. * cf12 * (1. - bp1)
     &           + (bad(iobj(i)) * gtblb(ith)) * bp1)

            tauc = ((so + mgn) - (so12 + mgn) * rio12)
            tauc = tauc
     &           * (u(iobj(i)) * gtblc(ith) * 4. * cf * (1. - bp1)
     &           + (bad(iobj(i)) * gtblb(ith)) * bp1)

            tau = tau + tauc + tauc12 + taub

            tauc = 0.0
            taub = 0.0
            tauc12 = 0.0

         end if

!     refresh  photon posiiton
         d = so * (1. - real(io)) + min(so, sp) * real(io)
         x0 = x0 + (d + mgn) * uxr
         y0 = y0 + (d + mgn) * uyr
         z0 = z0 + (d + mgn) * uzr

!     check upper or bottom boudary condition
         if(sign(1.0, uzr) * (z0 - zlim(flag)) .ge. 0.0) return

         x0 = x0 - (aint(x0 / xmax) - 0.5 + sign(0.5, x0)) * xmax
         y0 = y0 - (aint(y0 / ymax) - 0.5 + sign(0.5, y0)) * ymax
      end do
      return
      end
