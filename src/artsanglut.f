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
      subroutine artsanglut(wrkang, wrkphs, wrkcum, wrkpdf, nwrk, nlut,
     &     wrksca)
!************************************************************************
! Make LUTs of scattering angle for uniform distributions and others.
!
! Notes
!   wrkang: interpolated angles (degree)
!   wrkphs: interpolated phase functions
!   wrkpdf: normalized PDFs for spherical integration
!   wrkcum: normalized cumulative PDFs by backward integration (pi to 0)
!   wrksca: scattering angle LUT for uniform distribution
!*****

      parameter (pi = 3.1415926535898, rad = pi / 180.0)
      real wrkang(*), wrkphs(*), wrkcum(*), wrkpdf(*), wrksca(0:*)
      real*8 sum

! Integrate P (backward)
      sum = 0.0
      wrkpdf(nwrk) = 0.0
      wrkcum(nwrk) = 0.0
      a1 = rad * wrkang(nwrk)
      sina1 = sin(a1)
      do iwrk = nwrk - 1, 1, -1
         a0 = rad * wrkang(iwrk)
         sina0 = sin(a0)
         p = (a1 - a0) * (wrkphs(iwrk + 1)*sina1 + wrkphs(iwrk)*sina0)
         sum = sum + p
         wrkpdf(iwrk) = p
         wrkcum(iwrk) = sum
         a1 = a0
         sina1 = sina0
      end do

! Normalizations
!    (2*pi)/(4*pi)*sum{Q=0,pi;(P(Q)*sinQ*dQ)} = 1
      asum = 1.0 / sum
      f = 4.0 / sum
      do iwrk = 1, nwrk
         wrkphs(iwrk) = wrkphs(iwrk) * f
         wrkcum(iwrk) = wrkcum(iwrk) * asum
         wrkpdf(iwrk) = wrkpdf(iwrk) * asum
      end do
      wrkcum(1) = 1.0
      wrkcum(nwrk) = 0.0

! Scattering angle
      pdlt = 1.0 / real(nlut)
      iwrk = 1
      do ilut = 1, nlut - 1
         p = pdlt * real(nlut - ilut)
         iwrk = i_rvctrssrch1(wrkcum, p, iwrk, iwrk, nwrk - 1)
         if (wrkpdf(iwrk) .gt. 1.0e-6) then
            rat = (p - wrkcum(iwrk + 1)) / wrkpdf(iwrk)
            wrksca(ilut) = rat * wrkang(iwrk)
     &           + (1.0 - rat) * wrkang(iwrk + 1)
         else
            wrksca(ilut) = wrkang(iwrk)
         end if
      end do
      wrksca(0) = 0.0
      wrksca(nlut) = 180.0

      end
