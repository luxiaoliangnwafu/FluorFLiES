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

      subroutine artrroulette(w, wrr)
!*******************************************************************
! Russian roulette: survive or killed?
!     Modifed by H. Kobayashi Sep. 5, 2006
!****
      
      real*8 frnd

      if (w .lt. wrr * 0.5) then
         r = frnd()
         if (w .gt. wrr * r) then
            w = wrr             ! survive
         else
            w = 0.0             ! killed
         end if
      end if

      end
