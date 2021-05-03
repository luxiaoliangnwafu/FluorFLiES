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
      integer function ifunc_vctrbinsrch(grd, n, irvrs, dat)
!**********************************************************************
! Search grid number by binary search method.
!  grd(i) <= dat < grd(i+1) or grd(i) > dat >= grd(i+1).
!
! Note :
!  irvrs = 0    : grd(i) is increasing with i
!          1    : grd(i) is decreasing with i
!*****

      real grd(*)

      i0 = 0
      i1 = n + 1

      if (irvrs .le. 0) then
         do while(i1 .gt. i0 + 1)
            i = (i0 + i1) / 2
            if (dat .ge. grd(i)) then
               i0 = i
            else
               i1 = i
            endif
         end do
         ifunc_vctrbinsrch = i0
      else
         do while(i1 .gt. i0 + 1)
            i = (i0 + i1) / 2
            if (dat .lt. grd(i)) then
               i0 = i
            else
               i1 = i
            endif
         end do
         ifunc_vctrbinsrch = i0
      endif

      end
