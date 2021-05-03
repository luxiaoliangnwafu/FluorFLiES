!*************************************************************************
! init fluorescence weight
! wf (fluorescence weight), absorbed,
! leaf normal vector (vlx, vly, vlz), incident light vector(vix, viy, viz)
! incident light intensity (ilight)
!
!                             written Sicong Gao
!                             last modified 19/02/18
!**************************************************************************
      subroutine initfluorweight(wf, par, ilight, wl,
     &                           vlx, vly, vlz, vix, viy, viz)

      implicit none
      include "common.inc"
      ! input
      real vlx, vly, vlz, vix, viy, viz
      real wf, absorbed, ilight
      real wl, f_efficient

      ! work
      real vec, conv, par

      wf = 0
      if ((wl .lt. WL_PAR_START) .or. (wl .gt. WL_PAR_END)) then
          return
      end if

      conv = 1.d-8

      ! par is week, it cannot active fluorescence
      if (par .lt. conv) return

      !call getFluoyield(par, f_efficient)

      vec = abs(vlx * vix + vly * viy + vlz * viz)
      ! wf is decided by PAR (should include non-fluorescence and fluorescence)
      wf = vec * par * f_efficient

      !write(*,*) "wf = ", wf

      end
