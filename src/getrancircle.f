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
      subroutine getrancircle(r2min, r2, r, sinf, cosf)
!*************************************************************************
! Get a random point withing a unit circle
!****

      real*8 frnd
c      real rand

! Rejection method
      do i = 1, 200
         w1 = 1.0 - 2.0 * real(frnd())
         w2 = 1.0 - 2.0 * real(frnd())

         r2 = w1 * w1 + w2 * w2

         if (r2 .gt. r2min .and. r2 .le. 1.0) goto 1
      end do
      w1 = 1.0                  ! for unexpected case
      w2 = 0.0
      r2 = 1.0
 1    continue

! Azimuth
      r = sqrt(r2)
      cosf = w1 / r
      sinf = w2 / r

      end
