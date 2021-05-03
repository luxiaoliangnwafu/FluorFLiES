!**********************************************
!     3-Dimensional monte carlo 
!     calculation of the radiation field
!     in fish eye field
!     by H. Kobayashi
!     modified 04/12/08
!     last modified 08/04/07  
!**********************************************

      subroutine vegfeye(w, tgx, tgy, ichi, ikd)


      implicit none
      
      include 'common.inc'
      include 'math.inc'

      integer i, j
      integer ichi, ikd, ith, iph, sflag

      real x, y, z, tgx, tgy, tgz
      real ux, uy, uz, ph, th
      real w
      real tau, Id, rr, c, kzy
      real fcos, fsin, facos, fexp

      Id = 0.0

      do ith = 0, 89
         do iph = 0, 360
            ux = fsin(real(ith) * rad) * fcos(real(iph) * rad)
            uy = fsin(real(ith) * rad) * fsin(real(iph) * rad)
            uz = fcos(real(ith) * rad)

            tgz = 0.0

            call  vegtrace(tau, tgx,tgy, tgz, ux, uy, uz, sflag)           

            Id = w * fexp(-tau)
            Id = Id * real(sflag)
            feye(ith, iph) = feye(ith, iph) + Id
            
            Id = w
            rfeye(ith, iph) = rfeye(ith, iph) + Id
            
         end do
      end do

      return 
      end

      
