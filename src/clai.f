!*****************************************
! LAI calculation
!
!*****************************************

      subroutine clai(lai, plai, spn, cv)

      implicit none
     
      include 'common.inc'
      include 'math.inc'

      integer i, j, k, l, m, n, its
      integer ivox, idiv, fc, io1, io2, sflag
      real x, y, z, x1, y1, z1
      real ux, uy, uz, s, spn
      real lai, plai(100), tlai, intv
      real tobj(5), tobjb(5), cv, tau

      lai = 0.0
      cv = 0.0
      sflag = 1 
      tau = 0.0
c      spn = 0.1
      intv = 50. / res
      ux = 0.0 
      uy = 0.0
      uz = 1.0
      io1 = 1
      io2 = 1

      m = 1
      its = 0

      do i = 1, 100
         plai(i) = 0.0
      end do

      do k = 1, int(zmax / spn)
         z = (real(k) -0.5) * spn
         if(mod(k,int(1. / spn)) .eq. 0) m = m + 1
c         write(*,*) k

         do i = 1, int(xmax / spn)
            x = (real(i) -0.5) * spn

            do j = 1, int(ymax / spn)
               y = (real(j) -0.5) * spn
                
               x1 = aint(x / intv) 
               y1 = aint(y / intv)
               z1 = aint(z / intv)  
               ivox = ixmax * iymax * int(z1)
               ivox = ivox + int(y1) * iymax
               ivox = ivox + int(x1) + 1
               
               if(ndivs(ivox) .ne. 0)then
                  
                  s = 1.d5
                  do idiv = 1, ndivs(ivox)
                   
!     selected object number
                     n = divs(ivox, idiv)
                     
                     do l = 1, 5
                        tobj(l) = obj(n,l)
                     end do
               
!     define the branch dominant region
                     tobjb(1) = tobj(1)
                     tobjb(2) = tobj(2)
                     tobjb(3) = tobj(3) - tobj(4) * rb
                     tobjb(4) = tobj(4) * rb
                     tobjb(5) = tobj(5) * rb
                     
                     if(sobj(n) .eq. 1)then
                        call cones(s, x, y, z, ux, uy, uz, tobj, fc,io1)      
                        call cones(s, x, y, z, ux, uy, uz, tobjb,fc,io2)   
                     else if(sobj(n) .eq. 2)then
                        call cyls(s, x, y, z, ux, uy, uz, tobj, fc, io1)
                        call cyls(s, x, y, z, ux, uy, uz, tobjb,fc, io2)
                     else if(sobj(n) .eq. 3) then
                        call elpss(s ,x ,y ,z ,ux ,uy ,uz ,tobj, io1)
                        call elpss(s ,x ,y ,z ,ux ,uy ,uz,tobjb, io2)
                     else if(sobj(n) .eq. 5) then
                        call helpss(s ,x ,y ,z ,ux ,uy ,uz ,tobj,fc,io1)
                        call helpss(s ,x ,y ,z ,ux ,uy ,uz,tobjb,fc,io2)
                     end if
                     
                     if(io2. eq. 0) then
                        tlai = tlai + u(iobj(n)) * (1. - bp2)
                        its = its + 1
                     else if(io1 .eq. 0)then
                        tlai = tlai + u(iobj(n)) * (1. - bp1)
                        its = its + 1
                     end if

                  end do
                  
                  if (its .ne. 0) then
                     plai(m) = plai(m) + tlai / real(its)
                  end if
                  tlai = 0.0
                  its = 0
               end if
               
            end do
         end do
      end do
      
      do i = 1, 100
         plai(i) = plai(i) * spn * spn * spn / (xmax * ymax)
         lai = lai + plai(i)
      end do

!     crown cover calcualtion

      write(*,*) "Crown cover calculation ..."

      z = 0.01

!     dummy clumping factor and lea area density
      do i = 1, 5
         sbar(i) = 0.25
         u(i) = 1.0
      end do

      do i = 1, int(xmax / spn)
         x = (real(i) -0.5) * spn
         
         do j = 1, int(ymax / spn)
            y = (real(j) -0.5) * spn

            call  vegtrace(tau, x, y, z, ux, uy, uz, sflag)
            if(tau .ne. 0.0) cv = cv + 1.0

         end do
      end do

      cv = cv / real(int(xmax / spn) * int(ymax / spn))

      return
      end
