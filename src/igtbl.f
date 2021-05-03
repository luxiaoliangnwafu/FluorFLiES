!     ****************************
!     Ross-Nilson's G-function
!     2008/01/16
!     Hideki kobayashi
!     ***************************

      subroutine igtbl

      implicit none
      include 'common.inc'
      include 'math.inc'
 
      integer i,j
      real th,fgl,fpsi
      real dx,xm,xr,w(5),x(5),max,min,sc,sb,sf
      save w,x
      data w/0.2955242247, 0.2692667193, 0.2190863625,
     &     0.1494513491, 0.0666713443/
      data x/0.1488743389, 0.4333953941, 0.6794095682, 
     &     0.8650633666, 0.9739065285/

      real fsin,fcos,facos


!     Integration of the G-function
      max = pi / 2.0
      min = 0.0

      xm = 0.5 * (max + min)
      xr = 0.5 * (max - min)

      do i = 0, 180
         th = real(i) * rad
         sc = 0.0
         sb = 0.0
         sf = 0.0

         do j = 1, 5
            dx = xr * x(j)
!     for canopy            
            sc = sc + w(j) * (fgl(xm + dx, mc) * fpsi(th, xm + dx)
     &           + fgl(xm - dx, mc) * fpsi(th, xm - dx))
!     for branch area
            sb = sb + w(j) * (fgl(xm + dx, mb) * fpsi(th, xm + dx)
     &           + fgl(xm - dx, mb) * fpsi(th, xm - dx))
!     for forest floor
            sf = sf + w(j) * (fgl(xm + dx, mf) * fpsi(th, xm + dx)
     &           + fgl(xm - dx, mf) * fpsi(th, xm - dx))      
         end do
         gtblc(i) = sc * xr
         gtblb(i) = sb * xr
         gtblf(i) = sf * xr

      end do

      end

!     ********************************
!     psi funciton defined in Shultis and Myneni (1989)
!     this function has been validated 
!     2008/01/16
!     ********************************
      
      real function fpsi(th, thl)
      
      implicit none
      include 'math.inc'

      real absa,th,thl,pht
      real fsin,fcos,facos

      pht = -fcos(th) * fcos(thl)

      if(fsin(th) * fsin(thl) .le. 1.0d-5) then

         pht = pht / 1.0d-5

      else

         pht = pht / (fsin(th) * fsin(thl))

      end if

      absa = abs(pht)
      
      if(absa .gt. 1.0) then

         fpsi = abs(fcos(th) * fcos(thl))

      else

         pht = facos(pht)
         fpsi = fcos(th) * fcos(thl) * (2.0 * pht / pi - 1.0)
     &        + (2.0 / pi) * sin(th) * sin(thl) * sin(pht)

      end if     
      end


!     ********************************
!     leaf angle distribution function
!     (gl x sin(th))
!     using Bunnik (1978) definition
!     this function has been validated 
!     2008/01/16
!     ********************************

      real function fgl(thl, i)

      implicit none
      include 'math.inc'

      integer i
      real thl,ag(3),bg(3),cg(3),dg(3)
      real fsin,fcos

!     Preparation of the coeffifient
!     i=1: spherical, i=2 planophile i=3 erectrophile
      ag(1) = 0.0 
      bg(1) = 0.0
      cg(1) = 0.0
      dg(1) = 1.0

      ag(2) = 1.0
      bg(2) = 1.0
      cg(2) = 1.0
      dg(2) = 0.0

      ag(3) = 1.0 
      bg(3) = -1.0
      cg(3) = 1.0 
      dg(3) = 0.0

      fgl = (2.0 / pi) * 
     &     (ag(i) + bg(i) * fcos(2.0 * cg(i) * thl)) + dg(i) * fsin(thl)
      return
      end
