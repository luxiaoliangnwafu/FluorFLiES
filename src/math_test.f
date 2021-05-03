
      include "math.inc"

      integer i
      real x,y,z
      real fsin,fcos,facos,fexp


      call imath

      do i=-10000,10000
c         write(*,*)i, fsin(real(i)*0.0001),fcos(real(i)*0.0001)
c         write(*,*) i, fsin(real(i)*0.0001),sin(real(i)*0.0001)
c         write(*,*) i,fexp(real(i)),exp(-1.0*real(i))
         write(*,*) i,facos(real(i)*0.0001),acos(real(i)*0.0001)
      end do

      write(*,*)  fsin(0.785),fcos(0.785)
      stop
      end
