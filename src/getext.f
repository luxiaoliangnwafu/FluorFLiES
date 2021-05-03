! *********************************************************************
      subroutine getext(nmix,cflg, ext, Qext, Qext_ref, d, cbnz, ctnz,
     &     taur, ctaur, rat)
! ********************************************************************
      implicit none
      
      include 'common.inc'
      include 'math.inc'

      integer knmix,knang,knzext
      parameter (knmix=10, knang=2000, knzext=200)
      integer i,imix, nmix, cflg, cbnz, ctnz
      real ext(knmix, knzext), zmd(100), sfc(10), zsum
      real Qext(*), Qext_ref(*), d, rat(*)
      real taur,ctaur

      real fzpf

               
!     Calcualte the extinction coef. from total tau    
         do i=1, nz
            zmd(i) = 0.5 * (fzpf(zgrd(i),d) + fzpf(zgrd(i-1), d))
            zmd(i) = -d * log(zmd(i)) 
         end do 
        
!     Calculate the scale factor (sfc)         
         zsum = 0.0
         do i = 1, nz
            zsum = zsum + fzpf(zmd(i), d)
         end do
         
         do imix = 2, nmix
            sfc(imix - 1) = rat(imix - 1) * taur / zsum
         end do

         do imix = 2, nmix
            do i = 1, nz
               ext(imix,i) =
     &              (sfc(imix - 1) * fzpf(zmd(i), d))
     &              / (zgrd(i) - zgrd(i - 1))
!     Convert to ext from ext_ref
               ext(imix, i) = ext(imix, i) * (Qext(imix - 1) 
     &              / Qext_ref(imix - 1))
            end do    
         end do

         if(cflg .eq. 1)then
            do i = cbnz + 1, ctnz
            ext(nmix + cflg, i) = ctaur / (zgrd(ctnz) - zgrd(cbnz))
            ext(nmix + cflg, i) =
     &           ext(nmix + cflg, i) * (Qext(nmix) / Qext_ref(nmix))
            end do
         end if

         return
         end
      

c******************************************
      real function fzpf(z,d)
c******************************************
      implicit none

      real z,d

      fzpf = exp(-z / d)

      end
