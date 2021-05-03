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
      subroutine mc1doptics(ext, omg, phs, ang, knmix, nmix,iwl)
!*********************************************************************
! Make LUTs for phase function & scattering angle.
!*****

! Input
      integer nmix, knmix
      real ang(*), ext(knmix, *), omg(*), phs(knmix, *)

      include 'common.inc'
      include 'math.inc'

! Work variables
      parameter (knraw = 5000, knwrk = nlut1)
      real rawang(knraw), rawphs(knraw)
      real wrkc1(knraw), wrkc2(knraw), wrkc3(knraw)
      real wrkphs(knwrk), wrkang(knwrk), wrkcum(knwrk), wrkpdf(knwrk)
      real wrkcos(knwrk), wrksca(0:nlut)
      real delt

      integer iwl

! Initializations
      fsang = real(nlut) / pi
      nchi = 6
      chihi = 0.9
      chilo = 0.4
      fmax = 0.8
      chibin = (chihi - chilo) / real(nchi - 2)
      do ichi = 1, nchi - 1
         chigrd(ichi) = chihi - chibin * real(ichi - 1)
      end do
      chigrd(nchi) = -1.0

! Angle & cosine
      delt = 180.0 / real(nlut)
      nwrk = nlut + 1
      do iwrk = 1, nwrk
         a = delt * (iwrk - 1)
         wrkang(iwrk) = a
         wrkcos(iwrk) = cos(rad * a)
      end do

! Normalization of the phase functions
      do imix = 1, nmix
         do iang = 1, nang
            rawphs(iang) = phs(imix, iang)
         end do
         call nrmlzphsf(nang, ang, rawphs, 1.0)
         do iang = 1, nang
            phs(imix, iang) = rawphs(iang)
         end do
      end do

! Loop for layers
      do iz = 1, nz

!     Mix optical properties
         suma = 0.0
         sume = 0.0
         sums = 0.0
         do iang = 1, nang
            rawphs(iang) = 0.0
         end do
         do imix = 1, nmix
            e = ext(imix, iz)
            o = omg(imix)
            s = e * o
            suma = suma + e * (1.0 - o)
            sume = sume + e
            sums = sums + s

            do iang = 1, nang
               rawphs(iang) = rawphs(iang) + s * phs(imix, iang)
            end do
         end do

         extt1d(iz, 1) = sume
         abst1d(iz)    = suma

         f = 1.0 / max(1.0e-35, sums)
         do iang = 1, nang
            rawang(iang) = ang(iang)
            rawphs(iang) = rawphs(iang) * f
         end do
         nraw = nang

!     Make LUTs
         call artphsfintp(rawang, rawphs, wrkc1, wrkc2, wrkc3, wrkang,
     &        wrkphs, nraw, nwrk)
         call artsanglut(wrkang, wrkphs, wrkcum, wrkpdf, nwrk, nlut,
     &        wrksca)

         do ilut = 0, nlut
            pflut(ilut, iz) = wrkphs(ilut + 1)
            salut(ilut, iz) = wrksca(ilut) * rad
         end do

         pflut(nlut1, iz) = pflut(nlut, iz)
         salut(nlut1, iz) = salut(nlut, iz)

!     Truncation
         call phsfbiasym2(wrkang, wrkphs, nwrk, 90.0, iangs, g0fh,
     &        g0bh, g1a, g1fh, g1bh, g2a, g2fh, g2bh)
         cftmax = fmax * g0fh * g1fh**4 / real(nchi - 1)
         e = extt1d(iz, 1)
         s = e - abst1d(iz)
         do ichi = 1, nchi
            fd = cftmax * real(ichi - 1)
            if (fd .le. 0.0) then
               trulut(1, ichi, iz) = 1.0
               trulut(2, ichi, iz) = 0.0
               trulut(3, ichi, iz) = 1.0
               trulut(4, ichi, iz) = -1.0e+35
               trulut(5, ichi, iz) =  1.0e+35
               trulut(6, ichi, iz) = abs(g1a)
            else
               call phsftrunc5(wrkcos, wrkcum, wrkpdf, wrkang, nwrk,
     &              g1a, g2a, fd, ftf, ftb, g1t, g2t, angtf, angtb,
     &              iangtf, iangtb)
               ftt = 1.0 - (ftf + ftb)
               trulut(1, ichi, iz) = 1.0 - fd
               trulut(2, ichi, iz) = ftf
               trulut(3, ichi, iz) = ftt
               trulut(4, ichi, iz) = rad * angtf
               trulut(5, ichi, iz) = rad * angtb
               trulut(6, ichi, iz) = abs(g1t)
            end if
            extt1d(iz, ichi) = e - s * fd
         end do
      end do

      end
