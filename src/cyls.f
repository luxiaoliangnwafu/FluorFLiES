!************************************************
!     This routine calculates the distance between the (x,y,x)
!     and cylinder external boudary
!     written H. Kobayashi
!     last modified 01/08/14
!***************************************************
      
!***************************************************
      subroutine cyls(s, x, y, z, ux, uy, uz, tobj, face, io)
!************************************************
     
      implicit none 
      
      include 'common.inc'
      
!     parameter definition
!     photon position (x,y,z)
!     photon direction (ux, uy, uz)
!     cone apex(cx,cy,cz)
!     cone height h and radius bottom circle r
!     D =quadratic judgements
!     t1,t2,t3=solution
!     face: the face of cylinder 1=side;2=bottom circle;3=upper circle     
!     io : the positon in indide object=0 outside=1 

      integer face,io

      real x,y,z,ux,uy,uz
      real cx,cy,cz,h,r
      real circle
      real t1,t2,t3,t4,t,D
      real tobj(5)
      real s
     
!     we have to solve following eq.
!     a*t**2+b*t+c=0
!     for cone boundary
      real a,b,c
      
!     set initial value
      t = 0.0
      t1 = 1.e5
      t2 = 1.e5
      t3 = 1.e5
      t4 = 1.e5
      
!     crown related parameters
      cx = tobj(1)
      cy = tobj(2)
      cz = tobj(3)
      h = tobj(4)
      r = tobj(5)
      
!     calculate the quadratic parameters 
!     a*t^2+b*t+c=0
      a = ux**2 + uy**2
      a = max(a, 1.e-4)
      b = 2 * (ux * (x - cx) + uy * (y - cy))
      c = (x - cx)**2 + (y - cy)**2 - r**2
      
!     calculate D
      D = b**2 - 4.0 * a * c
      
      if(D .ge. 0) then
         t1 = (-b - sqrt(D)) / (2.0 * a)
         t2 = (-b + sqrt(D)) / (2.0 * a)
         
!     if t do not meet the following conditions, 
!     penalties are added in the solution.     
         if(((cz - h - 1.e-4) .le. (z + t1 * uz)).and.
     &        ((z + t1 * uz) .le. (cz + 1.e-4))) then
            if (t1 .lt. 0)then
               t1 = 1.e5               
            end if
         else
            t1 = 1.e5
         end if
         
         if(((cz - h - 1.e-4) .le. (z + t2 * uz)).and.
     &        ((z + t2 * uz) .le. (cz + 1.e-4)))then
            
            if (t2 .lt. 0)then
               t2 = 1.e5
            end if
         else 
            t2 = 1.e5   
         end if
         
!     end if of D judge      
!     else 
!     t1=100000
      end if
      
!     if ray-line cross the bottom circle of the cylinder crown. 
      uz = sign(max(abs(uz), 1.e-4), uz)     
      t3 = (cz - h - z) / uz
      circle = (x + t3 * ux - cx)**2 + (y + t3 * uy - cy)**2
      
      if((t3 .lt. 0) .or. (circle .gt. (r*r))) then
         t3 = 1.e5           
      end if
      
!     if ray-line cross the upper circle of the cylinder crown.
      uz = sign(max(abs(uz), 1.e-4), uz)
      t4 = (cz - z) / uz
      circle = (x + t4 * ux - cx)**2 + (y + t4 * uy - cy)**2
      
      if((t4 .lt. 0) .or. (circle .gt. (r*r))) then
         t4 = 1.e5           
      end if      

!     determination of the s      
      t = t1
      face = 1
      if(t .ge. t2)then
         t = t2 
      end if       
      if (t .ge. t3)then
         t = t3
         face = 2
      end if
      if (t .ge. t4)then
         t = t4
         face = 3
      end if
      
! there are no solution, if t=1.e5      
      s = t
      if(s .gt. 0.9d5) face = -1

! inout check

      if(((cz - h .lt. z) .and. (z .lt. cz))
     & .and. (((x - cx)**2 + (y - cy)**2) .lt. r**2))then
         io = 0
      else
         io = 1
      end if

      if(s .ge. 1.e5) io = 1

      return 
      end
