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
      integer function i_rvctrssrch1(xgrd, x, ix, ixmin, ixmax)
!**********************************************************************
! Find a grid number i, where 
!       xgrd(i) <= x  < xgrd(i+1) when xgrd is increasing vector, or
!       xgrd(i)  > x >= xgrd(i+1) when xgrd is decreasing vector,
!   by sequential search method with a given initial estimate, ix. 
! Output i will be in the range [ixmin, ixmax]. x should follow
!       xgrd(ixmin) <= x  < xgrd(ixmax+1), or
!       xgrd(ixmin)  > x >= xgrd(ixmax+1).
!*****

      real xgrd(*)

! Increasing vector
      if (xgrd(ixmin) .lt. xgrd(ixmax + 1)) then
         if (x .lt. xgrd(ix)) then
            do i = ix - 1, ixmin, -1
               if (x .ge. xgrd(i)) exit
            end do
            i = max(ixmin, i)
         else
            do i = ix + 1, ixmax + 1
               if (x .lt. xgrd(i)) exit
            end do
            i = min(ixmax, i - 1)
         end if

! Decreasing vector
      else
         if (x .ge. xgrd(ix)) then
            do i = ix - 1, ixmin, -1
               if (x .lt. xgrd(i)) exit
            end do
            i = max(ixmin, i)
         else
            do i = ix + 1, ixmax + 1
               if (x .ge. xgrd(i)) exit
            end do
            i = min(ixmax, i - 1)
         end if
      end if

      i_rvctrssrch1 = i

      end

