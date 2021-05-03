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
      subroutine phsftrunc5(cosa, cump, pdf, ang, nang, g1, g2,
     &     fd, ftf, ftb, g1t, g2t, angtf, angtb, iangtf, iangtb)
!***********************************************************************
! Dual-end truncation approximation (DTA):
!   Truncation of forward peak of phase function with a correction
!   for conservativation of 1st and 2nd moments of cosQ.
!
! Note:
!  1, cump() is in the reverse order (1 to 0)
!  2, The renormalized phase function is
!       P = phsf(ang) / (1 - ftf - ftb) if angtf < ang < angtb
!           0                           else
!****

      parameter (pi = 3.1415926535898, rad = pi / 180.0)

      real cosa(*), cump(*), pdf(*), ang(*)
      real*8 sum0, sum1, sum2

! No truncation
      if (fd .le. 0.0001) then
         ftf = 0.0
         ftb = 0.0
         g1t = g1
         g2t = g2
         angtf = 0.0
         angtb = 180.0
         iangtf = 1
         iangtb = nang
         return
      end if

! Initializations
      fdc = 1.0 - fd
      iangtf = ifunc_vctrbinsrch(cump, nang, 1, fdc)
      iangtb = nang
      g1t = (g1 - fd) / (1.0 - fd)
      g2t = (g2 - fd) / (1.0 - fd)
      g1t = max(-1.0, min(1.0, g1t))
      g2t = max( 0.0, min(1.0, g2t))
      sum0 = 0.0d0
      sum1 = 0.0d0
      sum2 = 0.0d0
      g1tlo = g1
      g1thi = g1
      g2tlo = g2
      g2thi = g2

! Integration over the delta part
      do iang = 1, iangtf
         w = pdf(iang)
         c0 = cosa(iang)
         c1 = cosa(iang + 1)
         sum0 = sum0 + w
         sum1 = sum1 + w * 0.5 * (c0 + c1)
         sum2 = sum2 + w * 0.5 * (c0 * c0 + c1 * c1)
         if (iang .ge. iangtf - 1) then
            g1tlo = g1thi
            g1thi = (g1 - sum1) / (1.0 - sum0)
         end if
      end do

! Loop
      iexit = -1
      do while(iexit .ne. 1)

! Truncation of backward region
         if (iexit .eq. 0) then
            if (iangtb .lt. nang .and. g2tlo .le. g2t) then
               iexit = 1
               if (abs(g2thi - g2tlo) .lt. 1.0e-6) then
                  ratb = 0.0
               else 
                  ratb = (g2t - g2tlo) / (g2thi - g2tlo)
               end if
               angtb = ang(iangtb) * (1.0 - ratb)
     &              + ang(iangtb + 1) * ratb
               ftb = cump(iangtb + 1) + pdf(iangtb) * (1.0 - ratb)
               w = -pdf(iangtb) * ratb
               c0 = cos(angtb * rad)
               c1 = cosa(iangtb)
            else
               g2thi = g2tlo
               iangtb = iangtb - 1
               angtb = ang(iangtb)
               ftb = cump(iangtb)
               w = pdf(iangtb)
               c0 = cosa(iangtb)
               c1 = cosa(iangtb + 1)
            end if
            sum0 = sum0 + w
            sum1 = sum1 + w * 0.5 * (c0 + c1)
            sum2 = sum2 + w * 0.5 * (c0 * c0 + c1 * c1)
            g1thi = (g1 - sum1) / (1.0 - sum0)
            if (sum0 .gt. 0.99) iexit = 1
            w = pdf(iangtf)
            c0 = cosa(iangtf)
            c1 = cosa(iangtf + 1)
            sum0a = sum0 - w
            sum1a = sum1 - w * 0.5 * (c0 + c1)
            g1tlo = (g1 - sum1a) / (1.0 - sum0a)
         else
            iexit = 0
         end if

! Truncation of forward region
         do while(g1t .lt. g1thi)
            iangtf = iangtf + 1
            g1tlo = g1thi
            w = pdf(iangtf)
            c0 = cosa(iangtf)
            c1 = cosa(iangtf + 1)
            sum0 = sum0 + w
            sum1 = sum1 + w * 0.5 * (c0 + c1)
            sum2 = sum2 + w * 0.5 * (c0 * c0 + c1 * c1)
            g1thi = (g1 - sum1) / (1.0 - sum0)
            if (sum0 .gt. 0.99) then
               angtf = ang(iangtf)
               ftf = 1.0 - cump(iangtf)
               goto 1
            end if
         end do
         do while(g1t .gt. g1tlo)
            iangtf = iangtf - 1
            g1thi = g1tlo
            w = -pdf(iangtf)
            c0 = cosa(iangtf)
            c1 = cosa(iangtf + 1)
            sum0 = sum0 + w
            sum1 = sum1 + w * 0.5 * (c0 + c1)
            sum2 = sum2 + w * 0.5 * (c0 * c0 + c1 * c1)
            g1tlo = (g1 - sum1) / (1.0 - sum0)
            if (iangtf .le. 1) then
               angtf = 0.0
               ftf = 0.0
               goto 1
            end if
         end do
         if (abs(g1thi - g1tlo) .lt. 1.0e-6) then
            ratf = 0.0
         else 
            ratf = (g1t - g1tlo) / (g1thi - g1tlo)
         end if
         angtf = ang(iangtf) * (1.0 - ratf) + ang(iangtf + 1) * ratf
         ftf = 1.0 - (cump(iangtf + 1) + pdf(iangtf) * (1.0 - ratf))
         w = pdf(iangtf) * (1.0 - ratf)
         c0 = cos(angtf * rad)
         c1 = cosa(iangtf + 1)
         sum0a = sum0 - w
         sum2a = sum2 - w * 0.5 * (c0 * c0 + c1 * c1)
         g2tlo = (g2 - sum2a) / (1.0 - sum0a)
      end do

 1    return

      end
