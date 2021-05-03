!     preparation of the look up table for the math functions
!      Hideki Kobayashi

      subroutine imath
      implicit none
      include 'math.inc'

      integer i, j
 
!     sin & cos
      do i = -65000, 65000
         tsin(i) = 0.0
         tcos(i) = 1.0
      end do
      do i = -62832, 62832
         tsin(i) = sin(real(i) * 0.0001)
         tcos(i) = cos(real(i) * 0.0001)
      end do

!     arccos
      do i = -10000, 10000
         tacos(i) = acos(real(i) * 0.0001)
      end do
    
!     exp (0.0 < x <80.0)
      do i = 0, 10000
         texp(i) = exp(-1.0 * 0.008 * real(i))
      end do

!     identity matrix (delta function)
      do i = 1, 6
         do j = 1, 6
            if(i .eq. j) then
               dlt(i, j) = 1.
            else
               dlt(i, j) = 0.
            end if
         end do
      end do
      

      return
      end
