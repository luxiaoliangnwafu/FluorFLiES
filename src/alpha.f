      subroutine alpha(ti,pi,to,po,cosa,a)
      
      implicit none

      real pi

      parameter(pi=3.14159265358979)

      cosa=sin(ti)*sin(to)*cos(pi-po)+cos(ti)*cos(to)

      if(abs(cosa).lt.0.99)then
         a=acos(cosa)
      else

      end if




      end
