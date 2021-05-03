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
      real function r_acos(x)
!********************************************************************
! Get single-precision acos(x) using LUT, where x should be in the
!  range from -1 to 1.
!****

      parameter (nlut = 1000000, zmax = 0.99)
      parameter (zbin = zmax / nlut, fz = nlut / zmax)
      parameter (pi = 3.1415926535898)
      real actab(0:nlut)
      save init, actab
      data init/0/

! Z in the range from 0 to 1
      if (x .ge. 0.0) then
         z = x
      else
         z = -x
      end if

! Small Z: use LUT
      if (z .lt. zmax) then
         if (init .eq. 0) then
            do ilut = 0, nlut
               actab(ilut) = acos(zbin * (ilut + 0.5))
            end do
            init = 1
         end if
         y = actab(int(z * fz))

! Large Z: 1st-order Newtonian approximation
      else if (z .lt. 1.0) then
         dz = 1.0 - z
         tmp = 30.0 + dz * dz
         y = sqrt(2.0 * dz) * ((tmp - 7.5 * dz) / (tmp - 10.0 * dz))

! Too large Z (invalid input!)
      else
         y = 0.0
      end if

! Result
      if (x .ge. 0.0) then
         r_acos = y
      else
         r_acos = pi - y
      end if

      end
