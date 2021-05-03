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
      subroutine nrmlzuvctr(ux, uy, uz)
!******************************************************************
! Renormalize a unit vector by using Newton method.
!   A given unit vector must satisfy that 
!    ux*ux + uy*uy + uz*uz is nearly equal to 1.
!*****

      fone = 0.5 * (3.0 - (ux * ux + uy * uy + uz * uz))
                       ! This is approximation by 1st-order Newton method

      ux = ux * fone
      uy = uy * fone
      uz = uz * fone

      end
