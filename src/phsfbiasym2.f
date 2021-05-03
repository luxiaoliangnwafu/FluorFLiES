!--License
!    Copyright (C) 2006 Hironobu Iwabuchi
!   
!    This file is part of MCARaTS.
!   
!    MCARaTS is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published
!    by the Free Software Foundation; either version 2 of the License,
!    or (at your option) any later version.
!   
!    MCARaTS is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!   
!    You should have received a copy of the GNU General Public License
!    along with this program; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
!    MA 02110-1301, USA.
!--End
      subroutine phsfbiasym2(ang, phs, nang, angs, iangs,
     &     g0fwd, g0bwd, g1, g1fwd, g1bwd, g2, g2fwd, g2bwd)
!*********************************************************************
! Get asymmetry factors of scattering phase function for each 
!  forward and backward regions.
!
! Note: ang() should be an increasing function.
!*****

      parameter (pi = 3.1415926535898, rad = pi / 180.0)
      real ang(*), phs(*)
      real*8 sum0, sum1, sum2

! Separation point & P interpolation
      iangs = ifunc_vctrbinsrch(ang, nang, 0, angs)
      iangs = max(1, min(nang - 1, iangs))
      dang = ang(iangs + 1) - ang(iangs)
      if (dang .gt. 1.0e-35) then
         phss = phs(iangs)
      else
         rat = (angs - ang(iangs)) / dang
         phss = phs(iangs) * (1.0 - rat) + phs(iangs + 1) * rat
      end if

! Forward region
      sum0 = 0.0
      sum1 = 0.0
      sum2 = 0.0
      a0 = rad * ang(1)
      sina0 = sin(a0)
      cosa0 = cos(a0)
      do iang = 1, iangs - 1
         a1 = rad * ang(iang + 1)
         sina1 = sin(a1)
         cosa1 = cos(a1)
         w = (a1 - a0) * (phs(iang + 1) * sina1 + phs(iang) * sina0)
         sum0 = sum0 + w
         sum1 = sum1 + w * 0.5 * (cosa0 + cosa1)
         sum2 = sum2 + w * 0.5 * (cosa0 * cosa0 + cosa1 * cosa1)
         a0 = a1
         sina0 = sina1
         cosa0 = cosa1
      end do
      a1 = rad * angs
      sina1 = sin(a1)
      cosa1 = cos(a1)
      w = (a1 - a0) * (phss * sina1 + phs(iangs) * sina0)
      sum0 = sum0 + w
      sum1 = sum1 + w * 0.5 * (cosa0 + cosa1)
      sum2 = sum2 + w * 0.5 * (cosa0 * cosa0 + cosa1 * cosa1)
      g0fwd = sum0
      if (sum0 .gt. 1.0d-35) then
         g1fwd = sum1 / sum0
         g2fwd = sum2 / sum0
      else
         g1fwd = 0.0
         g2fwd = 0.0
      end if

! Backward region
      sum0 = 0.0
      sum1 = 0.0
      sum2 = 0.0
      a1 = rad * ang(nang)
      sina1 = sin(a1)
      cosa1 = cos(a1)
      do iang = nang - 1, iangs + 1, -1
         a0 = rad * ang(iang)
         sina0 = sin(a0)
         cosa0 = cos(a0)
         w = (a1 - a0) * (phs(iang + 1) * sina1 + phs(iang) * sina0)
         sum0 = sum0 + w
         sum1 = sum1 + w * 0.5 * (cosa0 + cosa1)
         sum2 = sum2 + w * 0.5 * (cosa0 * cosa0 + cosa1 * cosa1)
         a1 = a0
         sina1 = sina0
         cosa1 = cosa0
      end do
      a0 = rad * angs
      sina0 = sin(a0)
      cosa0 = cos(a0)
      w = (a1 - a0) * (phs(iangs + 1) * sina1 + phss * sina0)
      sum0 = sum0 + w
      sum1 = sum1 + w * 0.5 * (cosa0 + cosa1)
      sum2 = sum2 + w * 0.5 * (cosa0 * cosa0 + cosa1 * cosa1)
      g0bwd = sum0
      if (sum0 .gt. 1.0d-35) then
         g1bwd = sum1 / sum0
         g2bwd = sum2 / sum0
      else
         g1bwd = 0.0
         g2bwd = 0.0
      end if

! Moments
      sum = g0fwd + g0bwd
      if (sum .gt. 1.0e-35) then
         g0fwd = g0fwd / sum
         g0bwd = g0bwd / sum
      else
         g0fwd = 0.0
         g0bwd = 0.0
      end if
      g1 = g0fwd * g1fwd + g0bwd * g1bwd
      g2 = g0fwd * g2fwd + g0bwd * g2bwd

      end
