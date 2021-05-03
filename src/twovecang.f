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
      subroutine twovecang(ux1, uy1, uz1, ux2, uy2, uz2, cosa, a)
!********************************************************************
! Computes angle and cosine between two unit vectors.
!****

      parameter (pi = 3.1415926535898)

      cosa = ux1 * ux2 + uy1 * uy2 + uz1 * uz2

      if (cosa .gt. -0.99 .and. cosa .lt. 0.99) then
         a = r_acos(cosa)
      else
         if (cosa .lt. 0.0) then
            dux = ux1 + ux2
            duy = uy1 + uy2
            duz = uz1 + uz2
         else
            dux = ux1 - ux2
            duy = uy1 - uy2
            duz = uz1 - uz2
         end if
         aa = dux * dux + duy * duy + duz * duz
         sinaa = aa * (1.0 - 0.25 * aa)
         a = sqrt(sinaa) * ((6.0 - 2.0 * sinaa) / (6.0 - 3.0 * sinaa))
                                ! approximation using Newtonian method
         if (cosa .lt. 0.0) a = pi - a

      end if

      end
