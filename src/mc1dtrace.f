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
      subroutine mc1dtrace(w, wq,x, y, z, ux, uy, uz, ftau, chi, iz,
     &     ikd, nscat, ichi,iwl)
!*********************************************************************
! Traces a trajectory in plane-parallel vertically-inhomogeneous
! atmosphere
!*****

      implicit none

      include 'common.inc'
      include 'math.inc'

!     Arguments (I/O)
      real w, wq, x, y, z, ux, uy, uz, ftau, chi
      real sinq, s, rlut, rilut, cosf, sinf, r, r2, p, tau
      real xr, yr, adf, pf, rat, q, uxr, uyr, uzr, cosq, utb
      real utf, ftt, ftf, fpath, path, absg, extm, absm
      real th,ph
      integer iz, ikd, nscat, ichi, irdc
      integer ixr, iyr, ilut
      integer i,j,iwl

! Work
      real*8 frnd

! Loop for layers
      do while(iz .ge. 1 .and. iz .le. nz)

         absg = absg1d(iz, ikd)
         extm = absg + extt1d(iz, ichi)
         absm = absg + abst1d(iz)
         if (uz .lt. 0.0) then
            path = (zgrd(iz - 1) - z) / uz
         else
            path = (zgrd(iz)     - z) / uz
         end if

!     Scattering events
         if (extm .gt. 1.0e-20) then
            fpath = ftau / extm
            do while(fpath .le. path)
                                ! motion to the collision point
               z = z + fpath * uz
!               if(z .lt. zmin) write(*,*) "z=",z
               call arthshift(x, ux, fpath, xmax)
               call arthshift(y, uy, fpath, ymax)
                                ! weight scaling
               w = w * (1.0 - absm / extm)

               call artrroulette(w, wrr)
               if (w .le. 0.0) then
                  return
               end if
                                ! properties
               ftf = trulut(2, ichi, iz)
               ftt = trulut(3, ichi, iz)
               utf = trulut(4, ichi, iz)
               utb = trulut(5, ichi, iz)
               chi = chi * trulut(6, ichi, iz)
               nscat = nscat + 1
               do while(chi .lt. chigrd(ichi))
                  ichi = ichi + 1
               end do
                                ! local estimates
               do irdc = 1, nrdc
                  uxr = uxrtab(irdc)
                  uyr = uyrtab(irdc)
                  uzr = uzrtab(irdc)
                  call twovecang(uxr, uyr, uzr, ux, uy, uz, cosq, q)

                  if (q .lt. utf .or. q .gt. utb) goto 1
                  rilut = fsang * q
                  ilut = int(rilut)
                  rat = rilut - real(ilut)
                  pf = (1.0 - rat) * pflut(ilut, iz)
     &                 + rat * pflut(ilut + 1, iz)
                  adf = pf / (4.0 * pi * abs(uzr) * ftt)
                  call mc1descape(x, y, z, uxr, uyr, uzr, iz, ichi,
     &                 ikd, xr, yr, tau)
                  ixr = xr / xmax * knxr + 1 ! give here the # of pixels
                  iyr = yr / ymax * knyr + 1
                  p = w * adf * exp(-tau)
                  ixr = min(ixr,knxr)
                  iyr = min(iyr,knyr)
                  prdcF(ixr, iyr, irdc) = prdcF(ixr, iyr, irdc) + p
                  prdcQ(ixr, iyr, irdc) = prdcQ(ixr, iyr, irdc) + p * wq

 1             end do

               do i = 1, nangc

                  uxr = uxrc(i)
                  uyr = uyrc(i)
                  uzr = uzrc(i)

                  call twovecang(uxr, uyr, uzr, ux, uy, uz, cosq, q)

                  if (q .lt. utf .or. q .gt. utb) goto 2
                  rilut = fsang * q
                  ilut = int(rilut)
                  rat = rilut - real(ilut)
                  pf = (1.0 - rat) * pflut(ilut, iz)
     &                 + rat * pflut(ilut + 1, iz)
                  adf = pf / (4.0 * pi * abs(uzr) * ftt)
                  call mc1descape(x, y, z, uxr, uyr, uzr, iz, ichi,
     &                 ikd, xr, yr, tau)
                  ixr = xr / xmax * knxr + 1 ! give here the # of pixels
                  iyr = yr / ymax * knyr + 1
                  p = w * adf * exp(-tau)
                  brf(2, i) = brf(2, i) + p

 2             end do

                                ! scattering
               call getrancircle(1.0e-12, r2, r, sinf, cosf)

               rlut = nlut * (ftf + ftt * r2)
               ilut = int(rlut)
               rat = rlut - real(ilut)

               s = (1.0 - rat) * salut(ilut, iz)
     &              + rat * salut(ilut + 1, iz)
               cosq = cos(s)
               sinq = sin(s)
               call artscat(ux, uy, uz, sinq, cosq, sinf, cosf)
               if (abs(uz) .lt. 1.0e-17) uz = 1.0e-17
                                ! new path
               ftau = -log(max(1.0e-35, real(frnd())))
               if (uz .lt. 0.0) then
                  path = (zgrd(iz - 1) - z) / uz
               else
                  path = (zgrd(iz)     - z) / uz
               end if
               extm = absg + extt1d(iz, ichi)
               fpath = ftau / extm
            end do
            ftau = max(0.0, ftau - extm * path)
         end if

!     Bounds
         if (uz .lt. 0.0) then  ! downward
            iz = iz - 1
            z = zgrd(iz)
         else                   ! upward
            z = zgrd(iz)
            iz = iz + 1
         end if
         call arthshift(x, ux, path, xmax)
         call arthshift(y, uy, path, ymax)

      end do

      end
