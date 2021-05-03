!     ****************************************************
!     Scattering direction
!     Hideki Kobayashi
!     We use the method of Frontier Technical Report No7
!     p90-91, based on Rejection method
!     ****************************************************
      subroutine scatvec(lr, lt, uxl, uyl, uzl, uxi, uyi, uzi,
     &                   uxo, uyo, uzo, m)

      implicit none
      include 'common.inc'
      include 'math.inc'

      integer m                 ! m leaf angle distribution
      real lr,lt,thi,phi,tho,pho
      real uxi,uyi,uzi,uxo,uyo,uzo
      real uxl,uyl,uzl
      real a,b
      real thl,phl,cosa,acosa
      real fsin,fcos,facos,fgl,fglm(3)
      real*8 frnd
      real rnd,ref

      fglm(1) = 1.0
      fglm(2) = 4.0 / pi
      fglm(3) = 4.0 / pi

!     step 1: determination of the leaf normal vector

 100  phl = 2. * pi * real(frnd())          ! phi direction
 101  thl = 0.5 * pi * real(frnd())         ! th direction
      rnd = real(frnd())
      ref = fgl(thl, m) / fglm(m)

      if(rnd .gt. ref) goto 101

      uxl = fsin(thl) * fcos(phl)
      uyl = fsin(thl) * fsin(phl)
      uzl = fcos(thl)

!     step 2: Adjustment to follow the |omega*omegaL|=cosa

      cosa = uxi * uxl + uyi * uyl + uzi * uzl
      rnd = real(frnd())

      if(rnd .gt. abs(cosa)) goto 100

!     step 3: Determination of the scattering direction on the leaf
!     Reflection or transmittion

      b = 2.0 * pi * real(frnd())
      a = sqrt(real(frnd()))
      a = facos(a)

      if(real(frnd()) .le. (lr / (lr + lt))) then ! reflection
         !write(*,*) "scatter: reflectance"
         if(cosa .ge. 0.0) then
            b = b + pi
            if(b .gt. (2. * pi)) b = b - 2. * pi
            a = pi - a
         end if
      else                      ! Transmission
         !write(*,*) "scatter: transmittance"
         if(cosa .lt. 0.0) then
            b = b + pi
            if(b .gt. (2. * pi)) b = b - 2. * pi
            a = pi - a
         end if
      end if

!     For tho and pho, coodinate transformation
      call trans(uxl, uyl, uzl, a, b, uxo, uyo, uzo)

      end


!     ****************************************************
!     get random leaf normal vector
!
!                             written Sicong Gao
!                             last modified 20/02/18
!     ****************************************************
      subroutine getleafnormvec(m, uxl, uyl, uzl)

      implicit none
      include 'common.inc'
      include 'math.inc'

      integer m                 ! m leaf angle distribution
      real uxl,uyl,uzl
      real thl,phl
      real fsin,fcos,facos,fgl,fglm(3)
      real*8 frnd
      real rnd,ref

      fglm(1) = 1.0
      fglm(2) = 4.0 / pi
      fglm(3) = 4.0 / pi

      phl = 2.0 * pi * real(frnd())

      do
        thl = 0.5 * pi * real(frnd())
        rnd = real(frnd())
        ref = fgl(thl, m) / fglm(m)

        if (rnd .le. ref) exit
      end do

      uxl = fsin(thl) * fcos(phl)
      uyl = fsin(thl) * fsin(phl)
      uzl = fcos(thl)

      end

!     ****************************************************
!     get scatter direction based on leaf normal vectors
!
!                             written Sicong Gao
!                             last modified 12/04/18
!     ****************************************************
      subroutine getscattervec(uxi, uyi, uzi,
     &                         uxl, uyl, uzl, uxo, uyo, uzo,
     &                         scatType)

      implicit none
      include 'common.inc'
      include 'math.inc'

      ! input
      real uxi, uyi, uzi, uxl, uyl, uzl, uxo, uyo, uzo
      integer scatType

      ! work
      real*8 frnd
      real a, b, cosa, facos

      cosa = uxi * uxl + uyi * uyl + uzi * uzl

      b = 2.0 * pi * real(frnd())
      a = sqrt(real(frnd()))
      a = facos(a)

      if (scatType .eq. 0) then ! reflection
         !write(*,*) "scatter: reflectance"
         if(cosa .ge. 0.0) then
            b = b + pi
            if(b .gt. (2. * pi)) b = b - 2. * pi
            a = pi - a
         end if
      else                      ! Transmission
         !write(*,*) "scatter: transmittance"
         if(cosa .lt. 0.0) then
            b = b + pi
            if(b .gt. (2. * pi)) b = b - 2. * pi
            a = pi - a
         end if
      end if

!     For tho and pho, coodinate transformation
      call trans(uxl, uyl, uzl, a, b, uxo, uyo, uzo)

      end
