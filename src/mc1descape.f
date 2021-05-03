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
      subroutine mc1descape(x, y, z, uxr, uyr, uzr, iz, ichi, ikd,
     &     xr, yr, tau)
!**********************************************************************
! Traces a trajectory in plane-parallel vertically-inhomogeneous
! atmosphere
!*****

      implicit none
      include 'common.inc'
      include 'math.inc'
! In
      real x, y, z, uxr, uyr, uzr
      integer iz, ichi, ikd
! Out
      real xr, yr, tau
! Work
      integer izr
      real zr, fpath,temp

! Ininial location
      xr = x
      yr = y
      zr = z

! TAU integration
      tau = 0.0

      if(ikd.eq.0) return

      if(uzr.ge.0.0)then
         do izr = iz, nz
            tau = tau + (zgrd(izr) - zr)
     &           * (absg1d(izr, ikd) + extt1d(izr, ichi))
            zr = zgrd(izr)
         end do
         tau = tau / uzr
      else
         do izr = iz,1,-1
            temp = (zr - zgrd(izr-1))
     &           * (absg1d(izr, ikd) + extt1d(izr, ichi))
            if (temp .ge. 1.0) then
              !write(*,*) "tau wrong"
            else
              tau = tau + temp
              zr = zgrd(izr-1)
            end if
         end do
         tau=tau/(-1.*uzr)
      end if

! Escape location
      fpath = (zr - z) / abs(uzr)
      call arthshift(xr, uxr, fpath, xmax)
      call arthshift(yr, uyr, fpath, ymax)

      end
