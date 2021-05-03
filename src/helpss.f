!************************************************
! This routine calculates the distance between the (x,y,x)
! and rotational half of elliptical sphere external boudary
!                             written H. Kobayashi
!                             last modified 05/1/13
!***************************************************

!***************************************************
      subroutine helpss(s ,x ,y ,z ,ux ,uy ,uz ,tobj, face, io)
!************************************************

      implicit none 

      include 'common.inc'

!     io : the positon in indide object=0 outside=1
      integer io, face
      
!     photon position (x,y,z)
!     photon direction (ux,uy,uz)
!     D =quadratic judgements
!     radius for z-axis : r1
!     radius for x - y plane: r2
!     face: the face of half elpsd 1=upper elpsd;2=bottom circle

      real x,y,z,z1,z2
      real t1,t2,t3,t,ux,uy,uz,D
      real tobj(5)
      real s, circle
      real a,b,c,elps
      real cx,cy,cz,r1,r2
      
!     crown related parameters
      cx = tobj(1)
      cy = tobj(2)
      cz = tobj(3)
      r1 = tobj(4)
      r2 = tobj(5)

      t = 0.0
      t1 = 1.e5
      t2 = 1.e5
      elps = 0.0
      s = 0.0

! calculate the quadratic parameters 
      a = r1 * r1 * ux * ux
      a = a + r1 * r1 * uy * uy
      a = a + r2 * r2 * uz * uz
      a = max(a, 0.0001)

      b = x * ux * r1 * r1 
      b = b -cx * ux * r1 * r1
      b = b +y * uy * r1 * r1
      b = b - cy * uy * r1 * r1
      b = b +z * uz * r2 * r2 
      b = b -cz * uz * r2 * r2      
      b = 2.0 * b

      c = r1 * r1 * (x - cx)**2
      c = c + r1 * r1 * (y - cy)**2 
      c = c + r2 * r2 * (z - cz)**2
      c = c - r2 * r2 * r1 * r1

! calculate D
      D = b**2 - 4.0 * a * c

      if(D .ge. 0)then

         t1 = (-b - sqrt(D)) / (2.0 * a)
         t2 = (-b + sqrt(D)) / (2.0 * a)

         z1 = z + t1 * uz
         z2 = z + t2 * uz

         if((t1 .lt. 0.) .or. (z1 .lt. cz)) then
            t1 = 1.e5
         end if
         if((t2 .lt. 0.) .or. (cz .gt. z2)) then
            t2 = 1.e5
         end if
      end if

!     if ray-line cross the bottom circle of the cone crown. 
      uz = sign(max(abs(uz), 1.e-4), uz)     
      t3 = (cz - z) / uz
    
      circle = (x + t3 * ux - cx)**2 + (y + t3 * uy - cy)**2
      
      if((t3 .lt. 0) .or. (circle .gt. (r2 * r2))) then
         t3 = 1.e5  
      end if

! determination of the s
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

! inout check 

      elps = (x - cx) * (x - cx) / r2**2
      elps = elps + (y - cy) * (y - cy) / r2**2
      elps = elps + (z - cz) * (z - cz) / r1**2

      if((elps .lt. 1.0) .and. (z .gt. cz))then
         io = 0
      else
         io = 1
      end if

      if(s .ge. 1.e5) io = 1

      return
      end

