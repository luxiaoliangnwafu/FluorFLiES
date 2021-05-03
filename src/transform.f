c****************************************
c     coordinate transformation
c
c      by Hideki Kobayashi
c****************************************

      subroutine trans(uxl,uyl,uzl,a,b,uxo,uyo,uzo)
      implicit none

      include 'common.inc'
      include 'math.inc'

      real uxl,uyl,uzl,uxo,uyo,uzo
      real a,b,c
      real fsin,fcos
      real sina,cosa,sinb,cosb,sint,sinp,cosp

      sina = fsin(a)
      cosa = fcos(a)
      sinb = fsin(b)
      cosb = fcos(b)
      sint = sqrt(uxl * uxl + uyl * uyl)
      sinp = uyl / sint
      cosp = uxl / sint
      
      if((uxl * uxl + uyl * uyl) .gt. 1.0d-15)then
         uxo = cosa * uxl
     &        + sina * (cosb * uzl * cosp - sinb * sinp)
         uyo = cosa * uyl
     &        + sina * (cosb * uzl * sinp + sinb * cosp)
         uzo = cosa * uzl - sina * cosb * sint
      else
         uxo = sina * cosb * sign(1., uzl)
         uyo = sina * sinb * sign(1., uzl)
         uzo = cosa * sign(1., uzl)
      end if

!     conversion to unit vector

      c = uxo * uxo + uyo * uyo + uzo * uzo
      c = sqrt(c)
      c =1. / c
      uxo = uxo * c
      uyo = uyo * c
      uzo = uzo * c
      
      end
