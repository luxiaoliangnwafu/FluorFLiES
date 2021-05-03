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
      subroutine artscat(ux, uy, uz, sinq, cosq, sinf, cosf)
!*********************************************************************
! Scattering (renewal of the direction). 
!*****

      parameter (eps = 1.0e-15)

      sinq0 = sqrt(ux * ux + uy * uy)

      if (sinq0 .gt. eps) then
         tmp = 1.0 / sinq0
         cosf0 = ux * tmp
         sinf0 = uy * tmp
         ux = cosq * ux + sinq * (cosf * uz * cosf0 - sinf * sinf0)
         uy = cosq * uy + sinq * (cosf * uz * sinf0 + sinf * cosf0)
         uz = cosq * uz - sinq * cosf * sinq0
      else if (uz .gt. 0.0) then
         ux = sinq * cosf
         uy = sinq * sinf
         uz = cosq
      else
         ux = -sinq * cosf
         uy = -sinq * sinf
         uz = -cosq
      end if

      call nrmlzuvctr(ux, uy, uz)

      end
