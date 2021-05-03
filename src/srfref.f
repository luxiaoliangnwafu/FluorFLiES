! ******************************************
!     scattering in the surface boudary
! ******************************************
      subroutine srfref(ux, uy, uz, sor, nscat, mode)

      implicit none
      include 'common.inc'
      include 'math.inc'

      integer nscat, mode
      real ux, uy, uz, w, sor
      real th, ph, facos, fsin, fcos
      real*8 frnd

c      write(*,*) ux,',',uy,',',uz,',',nscat,',',w,"srfref 1"

! lambertian reflection
      if(mode .eq. 1) then
         th = 0.5 * facos(1. - 2. * real(frnd()))
         ph = 2.0 * pi * real(frnd())
         ux = fsin(th) * fcos(ph)
         uy = fsin(th) * fsin(ph)
         uz = fcos(th)
         if(abs(uz) .lt. 0.0174524) uz = sign(0.0174524, uz)

!         call vegrroulette(w, epsi)
         nscat = nscat + 1
      end if

c      write(*,*) ux,',',uy,',',uz,',',nscat,',',w,"srfref 3"

 1001 format(F10.5,A2,F10.5,A2,F10.5,A2,I5,A2,A10)

      return
      end
