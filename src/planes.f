!************************************************
! This routine calculates the distance between the (x,y,z)
! and cube
!                             written H. Kobayashi
!                             last modified 08/04/01
!***************************************************

!***************************************************
      subroutine planes(s, x, y, z, ux, uy, uz, x1, y1, z1, face, intv)
!************************************************

      implicit none 

      include 'common.inc'
     
! face 1:x=x1, 2:x=x1+intv, 3:y=y1, 
!      4:y=y1+intv, 5:z=z1 (ground), 6:z=z1+intv (sky)
!     (x1, y1, z1) are  minumun potiion of cubic apex
!     intv(1):x, intv(2):y, intv(3):z

      integer face

      real x,y,z,x1,y1,z1,xt,yt,zt
      real ux,uy,uz,intv(3)
      real s,sb,su
      real conv
      real OUTPUT_WL
      OUTPUT_WL = 3

!     parallel condition
      conv = 1.d-6

!     for algorithm see Frontier Techinical Report NO7 p61-62

!     searh intersection
!     for y - z plane (x = x1,x = x1 + intv)

      if(abs(ux) .lt. conv) goto 100 ! for parallel condition exit

      sb = (x1 - x) / ux
      su = (x1 + intv(1) - x) / ux 
      if(sb .gt. su) then
         s = sb
         face = 1
      else
         s = su
         face = 2
      end if
      if (CURRENT_WL .ge. OUTPUT_WL) then
            write(*,*) "*sb, su, s = ", sb,su,s
      end if
      yt = y + s * uy
      if(( yt .ge. y1 ) .and. ( yt .le. y1 + intv(2) )) then
         zt = z + s * uz
         if(( zt .ge. z1 ) .and. ( zt .le. z1 + intv(3) )) then 
c            write(*,*)"y - z plane",sb,su,s
            return
         end if
      end if

!     for z - x plane (y = y1, y = y1 + intv)

 100  if(abs(uy) .lt. conv) goto 200 ! for parallel condition exit
      sb = (y1 - y) / uy
      su = (y1 + intv(2) - y) / uy 
      if(sb .gt. su) then
         s = sb
         face = 3
      else
         s = su
         face = 4
      end if
      if (CURRENT_WL .ge. OUTPUT_WL) then
            write(*,*) "**sb, su, s = ", sb,su,s
      end if

      zt = z + s * uz
      if(( zt .ge. z1 ) .and. ( zt .le. z1 + intv(3) )) then
         xt = x + s * ux
         if(( xt .ge. x1 ) .and. ( xt .le. x1 + intv(1) ))then
c            write(*,*)"z - x plane",sb,su,s
            return
         end if
      end if
      
!     for x - y plane (z = z1, z = z1 + intv)
 200  sb = (z1 - z) / uz
      su = (z1 + intv(3) - z) / uz
      if(sb .gt. su) then
         s = sb
         face = 5
      else
         s = su
         face = 6
      end if
      if (CURRENT_WL .ge. OUTPUT_WL) then
            write(*,*) "***sb, su, s = ", sb,su,s
      end if

      xt = x + s * ux
      if(( xt .ge. x1 ) .and. ( xt .le. x1 + intv(1) )) then
         yt = y + s * uy
         if(( yt .ge. y1 ) .and. ( yt .le. y1 + intv(2) )) then
c            write(*,*)"x - y plane",sb,su,s
            return
         end if

      end if


!     if cannot find the intersection due to the limiation of numerical calculation
!     move the photon position slightly then back to the top of this code 

      x = x + 0.1 * sign(1.0, (x1 + intv(1) * 0.5) - x)
      y = y + 0.1 * sign(1.0, (y1 + intv(2) * 0.5) - y)
      z = z + 0.1 * sign(1.0, (z1 + intv(3) * 0.5) - z)

c      write(*,*) "can't find cube intersection",x,y,z,ux,uy,uz
c      stop

      end

