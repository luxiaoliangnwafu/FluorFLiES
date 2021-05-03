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
      subroutine nrmlzphsf(nang, ang, pf, sumnorm)
!*******************************************************************
! Normalize phase function to -sumnorm-.
!*****

      parameter (pi = 3.1415926535898, rad = pi / 180.0)
      real ang(*), pf(*)
      real*8 sum

      sum = 0.0
      do iang = nang - 1, 1, -1
         a0 = rad * ang(iang)
         a1 = rad * ang(iang + 1)
         sum = sum + (a0 - a1) * 0.5
     &        * (pf(iang + 1) * sin(a1) + pf(iang) * sin(a0))
      end do
      sum = abs(sum)

      f = 2.0 * sumnorm / sum
      do iang = 1, nang
         pf(iang) = pf(iang) * f
      end do

      end
