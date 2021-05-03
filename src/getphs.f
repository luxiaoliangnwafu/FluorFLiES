      subroutine getatm(wl0, nkd, wkd, ext, rfname)
      
      implicit none

      include 'common.inc'
      include 'math.inc'

      integer i, j, iz, nkd, ikd, k, nl
      integer knmix,knzext
      parameter(knmix = 10, knzext = 200)
      integer inz,iwave
      real wl0,wl1,wl2,wkd(knkd)
      real ext(knmix, knzext)
      real ext_back(knmix, knzext)
      real absg1d_back(knz,knkd)
      real zmed

      character*81 rfname

      open(12,file = rfname)
      do i = 1, 740
         read(12, *) 
         read(12, *) iwave, wl1, wl2
         read(12, *) 
         
         do iz = 1, nz
            read(12,*) inz, ext(1,iz), (absg1d(iz,j), j = 1, nkd)
     &           ,(wkd(k), k = 1, nkd)
!            write(*,*) iz,inz, ext(1,iz), (absg1d(iz,j), j = 1, nkd)
!     &           ,(wkd(k), k = 1, nkd)
         end do

!     the ext(1,iz) and absg1d values are adjusted according to the 
!     actual height level
!     caution!! this is very rough mapping, I didn't do a extrapolation
!     I just used the values from closest layers.
         do iz = 1, nz
            ext_back(1, iz) = ext(1, iz)
!            ext(1, iz) = ext(1, klayer(iz))
            nl = klayer(iz)
            zmed = 0.5 * (zgrd(iz) + zgrd(iz - 1)) 
            ext(1, iz) = (zgrdm(nl + 1) - zmed) * ext(1, nl)
            ext(1, iz) 
     &           = ext(1, iz) 
     &           + (zmed - zgrdm(nl)) * ext(1, nl + 1)
            ext(1, iz) = ext(1, iz) / (zgrdm(nl + 1) - zgrdm(nl))

            if(iz .eq. nz)   ext(1, iz) = ext(1, klayer(iz)) 

            do j = 1, 3
               absg1d_back(iz, j) = absg1d(iz, j)
!               absg1d(iz, j) = absg1d(klayer(iz), j)
               absg1d(iz, j) = (zgrdm(nl + 1) - zmed) * absg1d(nl, j)
               absg1d(iz, j) 
     &              = absg1d(iz, j) 
     &              + (zmed - zgrdm(nl)) * absg1d(nl + 1, j)
               absg1d(iz, j) = absg1d(iz, j)/(zgrdm(nl + 1) - zgrdm(nl))

               if(iz .eq. nz)  absg1d(iz, j) = absg1d(klayer(iz), j) 
            end do

!            write(*,*) iz,klayer(iz),ext(1, iz) ,ext_back(1, iz),
!     absg1d(iz, 1), absg1d_back(iz, 1),
!     &           zmed, zgrdm(nl + 1),zgrdm(nl)
!     &           absg1d(iz, 3), absg1d_back(iz, 3)

         end do

      if(wl2 .ge. wl0 .and. wl0 .ge. wl1) then
         wkd(nkd) = 1.000001    ! because of security         
         exit
      end if 
        
      end do
      close(12)
      return
      end

c **********************************************************************
      subroutine getphs
     &   (wl0, span, nmix, re, Qext_ref, Qext, omg,
     &     G, Qabs, ang, phs, fname, imode)
c *********************************************************************
      implicit none

      include 'common.inc'
      include 'math.inc'

      integer i, j, iang
      integer imix, nmix
      integer knmix, imode
      parameter (knmix = 10)
      real fac, span
      real re(knmix),Qext_ref(knmix),Qext(knmix),G(knmix),Qabs(knmix)
      real wref,crref(knmix),ciref(knmix),cr(knmix),ci(knmix)
      real omg(knmix),phs(knmix,*),ang(*)
      real wl0, wl1, wl2
      real dQext, dG, dQabs, dphs(1000)
      real dcr,dci, domg
      real RvDry, ReDry, rhoDry, RvWet, ReWet, rhoWet
      real dum
      
      character*81 fname(knmix)
      
      do imix = 2, nmix
         open(13,file = fname(imix - 1))
         if(imode .eq. 1)then
            do i = 1, 740
               read(13,*)
               read(13,*) wl1, dum
               read(13,*)
               read(13,*)
               read(13,*) RvDry, ReDry, rhoDry, RvWet, ReWet, rhoWet
               read(13,*) 
               read(13,*) Qext_ref(imix - 1), Qext(imix - 1), 
     &              Qabs(imix - 1), omg(imix), G(imix - 1)
               read(13,*) 
               read(13,*) nang
               read(13,*)
               
               do j = 1, nang
                  read(13,*) ang(j), phs(imix, j)
               end do
               
               if(wl1 + 0.005 .gt. wl0) then
                  
                  read(13,*)
                  read(13,*) wl2, dum
                  read(13,*)
                  read(13,*)
                  read(13,*) RvDry, ReDry, rhoDry, RvWet, ReWet, rhoWet 
                  read(13,*) 
                  read(13,*) Qext_ref(imix - 1), dQext, dQabs, domg, dG
                  read(13,*) 
                  read(13,*) nang
                  read(13,*)
                  
                  do j = 1, nang
                     read(13,*) ang(j), dphs(j)
                  end do
                  
                  Qext(imix - 1) = ((wl0 - wl1) * dQext + 
     &                 (wl2 - wl0) * Qext(imix - 1)) * (1. / 0.005)
                  Qabs(imix - 1) = ((wl0 - wl1) * dQabs + 
     &                 (wl2 - wl0) * Qabs(imix - 1)) * (1. / 0.005)
                  omg(imix) = ((wl0 - wl1) * domg +
     &                 (wl2 - wl0) * omg(imix)) *  (1. / 0.005)
                  G(imix - 1) = ((wl0 - wl1) * dG + 
     &                 (wl2 - wl0) * G(imix - 1)) * (1. / 0.005)
                  
                  do j = 1, nang
                     phs(imix, j) = ((wl0 - wl1) * dphs(j) +
     &                    (wl2 - wl0) * phs(imix, j)) * (1. / 0.005)
                  end do
                  exit
               end if
            end do
         else
            do i=1, 53
               read(13,*)
               read(13,*) wl1, wl2
               read(13,*)
               read(13,*)
               read(13,*) RvDry, ReDry, rhoDry, RvWet, ReWet, rhoWet
               read(13,*) 
               read(13,*) Qext_ref(imix - 1),Qext(imix - 1),
     &              Qabs(imix - 1), omg(imix), G(imix - 1)
               read(13,*) 
               read(13,*) nang
               read(13,*)
               
               do j = 1, nang
                  read(13,*) ang(j), phs(imix, j)
               end do
            
               if(wl0 .gt. wl1 .and. wl0 .lt. wl2)then
                  exit
               end if
            end do   
         end if
         close(13)         
      
      end do

c     make rayleigh scattering albedo
      omg(1) = 1.0              ! Rayleigh scatteing albedo
      fac = 3.0 / (16.0 * pi)
      do iang = 1, nang         ! Rayleigh phase function
         phs(1, iang) = fac * (1.0 + cos(ang(iang) * rad)**2)
      end do
      return
      end
