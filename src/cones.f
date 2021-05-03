!************************************************
! This routine calculates the distance between the (x,y,x)
! and cone external boudary
!                             written H. Kobayashi
!                             last modified 04/11/15
!***************************************************

!***************************************************
      subroutine cones(s, x, y, z, ux, uy, uz, tobj, face, io)
!************************************************

      implicit none 

      include 'common.inc'

! parameter definition
! inout : the positon in indide object=0 outside=1
      integer face,io
      
! photon position (x,y,z)
! photon direction (ux, uy, uz)
! cone apex(cx,cy,cz)
! cone height h and radius bottom circle r
! rp=(r/h)
!     face: the face of cylinder 1=side;2=bottom circle
! D =quadratic judgements
! t1,t2,t3=solution

      real x,y,z,ux,uy,uz
      real cx,cy,cz,h,r,rp
      real circle
      real t1,t2,t3,t,D
      real tobj(5)
      real s

! we have to solve following eq.
! a*t**2+b*t+c=0
! for cone boundary
      real a,b,c

!     set initial value
      t = 0.0
      t1 = 1.e5
      t2 = 1.e5
      t3 = 1.e5
      
!     crown related parameters
      cx = tobj(1)
      cy = tobj(2)
      cz = tobj(3)
      h = tobj(4)
      r = tobj(5)
      rp = r / h
      
      
!     calculate the quadratic parameters 
!     a*t^2+b*t+c=0
      a = ux**2 + uy**2 - (rp * uz)**2
      a = sign(max(abs(a), 0.0001), a)
      
      b = 2. * (ux * (x - cx) + uy * (y - cy) - uz * (z - cz) * rp**2)
      c = (x - cx)**2 + (y - cy)**2 - (rp * (z - cz))**2
      
!     calculate D
      D = b**2 - 4.0 * a * c
      
      if(D .ge. 0)then
         
         t1 = (-b - sqrt(D)) / (2.0 * a)
         t2 = (-b + sqrt(D)) / (2.0 * a)
         
!     if t do not meet the following conditions, 
!     penalties are added in the solution.     
         if(((cz - h - 1.e-4) .le. (z + t1 * uz)).and.
     &        ((z + t1 * uz) .le. cz + 1.e-4)) then
            
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

      end if
      
!     if ray-line cross the bottom circle of the cone crown. 
      uz = sign(max(abs(uz), 1.e-4), uz)     
      t3 = (cz - h - z) / uz
      
      circle = (x + t3 * ux - cx)**2 + (y + t3 * uy - cy)**2
      
      if((t3 .lt. 0) .or. (circle .gt. (r * r))) then
         t3 = 1.e5  
      end if
      
!     determination of the sb      
      t = t1
      face = 1
      if(t .ge. t2)then
         t = t2 
      end if  
     
      if (t .ge. t3)then
         t = t3
         face = 2
      end if
      
!     there are no solution, if t=1.e5    
      s = t
      if(s .gt. 0.9d5) face = -1
      
!     inout check     
      if(((cz - h .lt. z) .and. (z .lt. cz)))then
         if((x - cx)**2 + (y - cy)**2 .lt. (rp * (cz - z))**2)then
            io = 0
         else
            io = 1
         end if
      else
         io = 1
      end if

! to prevent numerical error
      if(s. ge. 1.e5) io = 1
 
      return 
      end
      


      
