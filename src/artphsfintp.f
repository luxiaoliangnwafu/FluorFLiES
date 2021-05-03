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
      subroutine artphsfintp(rawang, rawphs, wrkc1, wrkc2, wrkc3,
     &     wrkang, wrkphs, nraw, nwrk)
!**********************************************************************
! Interpolation of phase function by natural cubic spline
!****

      parameter (eps = 1.0e-30)
      real rawang(*), rawphs(*), wrkc1(*), wrkc2(*), wrkc3(*)
      real wrkang(*), wrkphs(*)

! Natural cubic spline coefficients
      wrkc1(1) = (rawphs(2) - rawphs(1)) / (rawang(2) - rawang(1))
      wrkc1(2) = (rawphs(3) - rawphs(2)) / (rawang(3) - rawang(2))
      wrkc2(2) = 2.0 * (rawang(3) - rawang(1))
      wrkc3(2) = wrkc1(2) - wrkc1(1)
      do i = 3, nraw - 1
         h1 = rawang(i + 1) - rawang(i)
         h0 = rawang(i) - rawang(i - 1)
         tmp = h0 / wrkc2(i - 1)
         wrkc1(i) = (rawphs(i + 1) - rawphs(i)) / h1
         wrkc2(i) = 2.0 * (h0 + h1) - h0 * tmp
         wrkc3(i) = wrkc1(i) - wrkc1(i - 1) - wrkc3(i - 1) * tmp
      end do
      sig0 = wrkc3(nraw - 1) / wrkc2(nraw - 1)
      h = rawang(nraw) - rawang(nraw - 1)
      wrkc1(nraw - 1) = wrkc1(nraw - 1) - h * 2.0 * sig0
      wrkc2(nraw - 1) = 3.0 * sig0
      wrkc3(nraw - 1) = -sig0 / h
      sig1 = sig0
      do i = nraw - 2, 2, -1
         h = rawang(i + 1) - rawang(i)
         sig0 = (wrkc3(i) - h * sig1) / wrkc2(i)
         wrkc1(i) = wrkc1(i) - h * (sig1 + 2.0 * sig0)
         wrkc2(i) = 3.0 * sig0
         wrkc3(i) = (sig1 - sig0) / h
         sig1 = sig0
      end do
      h1 = rawang(2) - rawang(1)
      wrkc1(1) = wrkc1(1) - h1 * sig1
      wrkc2(1) = 0.0
      wrkc3(1) = sig1 / h1

! Interpolation
      iraw = 1
      do iwrk = 1, nwrk
         x = wrkang(iwrk)
         iraw = i_rvctrssrch1(rawang, x, iraw, iraw, nraw - 1)
         dx = x - rawang(iraw)
         y = rawphs(iraw) + (wrkc1(iraw) + (wrkc2(iraw) + wrkc3(iraw)
     &        * dx) * dx) * dx
         wrkphs(iwrk) = max(0.0, y)
c         if (y .lt. 1.0e-8 .or. y .gt. 1.0e10)
c     &        write (*,*) 'artphsfintp: Phase function is invalid,', y
      end do
      
      end
