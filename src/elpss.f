!************************************************
! This routine calculates the distance between the (x,y,x)
! and rotational elliptical sphere external boudary
!                             written H. Kobayashi
!                             last modified 05/1/13
!***************************************************

!***************************************************
      subroutine elpss(s ,x ,y ,z ,ux ,uy ,uz ,tobj, io)
!************************************************

      implicit none 

      include 'common.inc'

!     io : the positon in indide object=0 outside=1
      integer io
      
!     photon position (x,y,z)
!     photon direction (ux,uy,uz)
!     D =quadratic judgements
!     radius for z-axis : r1
!     radius for x - y plane: r2
      real x,y,z
      real t1,t2,t,ux,uy,uz,D
      real tobj(5)
      real s
      real a,b,c,elps
      real cx,cy,cz,r1,r2
      
!     crown related parameters
      cx = tobj(1)
      cy = tobj(2)
      cz = tobj(3)
      r1 = tobj(4)
      r2 = tobj(5)

      t = 0.0
      t1 = 1.d5
      t2 = 1.d5
      elps = 0.0
      s = 0.0

! calculate the quadratic parameters 
      a = r1 * r1 * ux * ux
      a = a + r1 * r1 * uy * uy
      a = a + r2 * r2 * uz * uz
      a=max(a,0.0001)

      b = x * ux * r1 * r1 
      b = b -cx * ux * r1 * r1
      b = b +y * uy * r1 * r1
      b = b - cy * uy * r1 * r1
      b = b +z * uz * r2 * r2 
      b = b -cz * uz * r2 * r2      
      b = 2.0 * b

      c = r1 * r1 * (x-cx)**2
      c = c + r1 * r1 * (y-cy)**2 
      c = c + r2 * r2 * (z-cz)**2
      c = c - r2 * r2 * r1 * r1

! calculate D
      D = b**2 - 4.0 * a * c

      if(D .ge. 0)then
         
         t1 = (-b - sqrt(D)) / (2.0 * a)
         t2 = (-b + sqrt(D)) / (2.0 * a)

         if(t1 .lt. 0) t1 = 1.d5
         if(t2 .lt. 0) t2 = 1.d5
         
      end if

! determination of the s
      s = min(t1,t2)

! inout check 

      elps = (x - cx) * (x - cx) / r2**2
      elps = elps + (y - cy) * (y - cy) / r2**2
      elps = elps + (z - cz) * (z - cz) / r1**2

      if(elps .lt. 1.0 .and. s. lt. 1.d5)then
         io = 0
      else
         io = 1
      end if
     
      return
      end

