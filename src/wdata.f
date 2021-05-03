      subroutine wdata(np, nwl, stype, imode, cmode, amode,
     &     span, RF, RQ, cosq0, wl0, wls, tflx, bflx, dflx,
     &     tpfd, bpfd, dpfd, scmpf, scmpp,rflx, rbflx, rdflx,
     &     fmode)

      implicit none
      include 'common.inc'
      include 'math.inc'

      integer i, j, k, ith, iph
      integer np, nwl, stype, cmode, amode, imode
      real tm, pixnp, pixnpa

      real RF, RQ, wls, wl0, cosq0, th, ph
      real*8 tflx, bflx, dflx, tpfd, bpfd, dpfd, rflx, rbflx, rdflx
      real*8 scmpf(3,1000),scmpp(3,1000)
      real span(1000)
      real app(100), ffsum, sfsum
      real facos, eye(180, 180),T(5),W(5), S(5),aaT(0:90), ave, elai

      real aveF, aveQ, sdvF, sdvQ, sumF, sumQ, sum2F, sum2Q, pF, pQ
      real pFmax, pQmax
      real phi, uxc

      ! SIF
      integer fmode

      if(np .eq. -4) goto 200

      do k = 1, 100
         app(k) = 0.0
      end do

!     summary of the radiative flux

!     summary calculation of flux
      pixnp = real(np) / real(size * size)
      tm = real(tflx) / real(np)
c      write(*,*) tpfd / real(np)

      do i = 1, nwl
         scmpf(1, i) = 100. * scmpf(1, i) / tflx
         scmpp(1, i) = 100. * scmpp(1, i) / tpfd
         scmpf(2, i) = 100. * scmpf(2, i) / bflx
         scmpp(2, i) = 100. * scmpp(2, i) / bpfd
         scmpf(3, i) = 100. * scmpf(3, i) / dflx
         scmpp(3, i) = 100. * scmpp(3, i) / dpfd
      end do

      tflx = RF * abs(cosq0) * tflx / real(np)
      bflx = RF * abs(cosq0) * bflx / real(np)
      dflx = RF * abs(cosq0) * dflx / real(np)

      rflx = RF * abs(cosq0) * rflx / real(np)
      rbflx = RF * abs(cosq0) * rbflx / real(np)
      rdflx = RF * abs(cosq0) * rdflx / real(np)

      tpfd = RQ * abs(cosq0) * tpfd / real(np)
      bpfd = RQ * abs(cosq0) * bpfd / real(np)
      dpfd = RQ * abs(cosq0) * dpfd / real(np)

      cfpr =  RQ * abs(cosq0) * cfpr / real(np)
      bfpr =  RQ * abs(cosq0) * bfpr / real(np)
      ffpr =  RQ * abs(cosq0) * ffpr / real(np)
      sfpr =  RQ * abs(cosq0) * sfpr / real(np)
      tfpr = cfpr + bfpr + ffpr

      tfpr = tfpr / real(tpfd)
      cfpr = cfpr / real(tpfd)
      bfpr = bfpr / real(tpfd)
      ffpr = ffpr / real(tpfd)
      sfpr = sfpr / real(tpfd)

!     summary of 3D far etc

      th = facos(cosq0)
      ith = int(th * 180. / pi)

      do k = int(zmax) + 1, 1, -1
         do i = 1, size
            do j = 1, size
               ap(i ,j, k) = ap(i, j, k) * RQ * abs(cosq0) / pixnp
               apd(i, j, k) =  apd(i, j, k) * RQ * abs(cosq0) / pixnp
               apb(i, j, k) = apb(i, j, k) * abs(cosq0) / pixnp
               apb(i, j, k) = apb(i, j, k) * (1. / gtblc(ith))
               apb(i, j, k) = min(apb(i, j, k), 1.)
               app(k) = app(k) + ap(i, j, k)
            end do
         end do
         app(k) = app(k) / (real(tpfd) * real(size * size))
         apnp(k) = apnp(k) * RQ * abs(cosq0) / real(np) / real(tpfd)
      end do

!      write(*,*) tflx,tpfd,tfpr
!      stop

!     summary of forest floor far etc
      do i = 1, size
         do j = 1, size
            apf(i ,j) = apf(i, j) * RQ * abs(cosq0) / pixnp
            apfd(i, j) =  apfd(i, j) * RQ * abs(cosq0) / pixnp
         end do
      end do

!     summary of forest floor and surface downward flux
      ffsum = 0.0
      sfsum = 0.0
      do i = 1, size
         do j = 1, size
            sfdir(i ,j) = sfdir(i, j) * RQ * abs(cosq0) / pixnp
            sfdif(i, j) = sfdif(i, j) * RQ * abs(cosq0) / pixnp
            ffdir(i ,j) = ffdir(i, j) * RQ * abs(cosq0) / pixnp
            ffdif(i, j) = ffdif(i, j) * RQ * abs(cosq0) / pixnp
            ffsum = ffsum + ffdir(i, j) + ffdif(i, j)
            sfsum = sfsum + sfdir(i, j) + sfdif(i, j)
         end do
      end do
      ffsum = ffsum / real(size * size)
      sfsum = sfsum / real(size * size)

!     writing
      open(12, file = "flxsum_"//trim(FILE_O_LABEL)//".txt")

      write(12,*) "-- summary of radiative quantity --"
      write(12,*)

      write(12,*) "-- simulation mode --"
      if(amode .eq. 1) write(12,'(A25,A20)')
     &     " Atmospheric module:", "Yes"
      if(amode .eq. 2) write(12,'(A25,A20)')
     &     " Atmopsheric module:", "No"

      if(amode .eq. 1) then
         if(imode .eq. 1) write(12,'(A25,A20,F10.5)')
     &        "Spectral domain:", "single wavelength",wl0
         if(imode .eq. 2) write(12,'(A25,A20)') "Spectral domain:","PAR"
         if(imode .eq. 3) write(12,'(A25,A20)')
     &        "Spectral domain:","Shortwave total"
      end if

      if(stype .eq. 1) write(12,'(A25,A20)')
     &     "Surface type:", "lambertian"
      if(stype .eq. 2) write(12,'(A25,A20)')
     &     "Surface type:", "3 D canopy"

      write(12,*)

      write(12,*) "-- Solar zenith angle, Elv(m), atm. transmittance --"
      write(12,'(F12.2,A4,F12.5,A4,F10.5)')
     &     (pi - acos(cosq0)) / rad ," ",zmin," ", tm

      write(12,*)
      write(12,*) "-- Downward Irradiance at TOC --"

      if(amode .eq. 2) then
         write(12,*) "warning ! irradiance & PPFD are "
         write(12,*) " normalized by 1000.*cos(th) !!"
         write(12,*)
      end if

      write(12,*) "Energy unit (W/m2)"
      write(12,'(4(A12))') "Total","Beam","Diffuse","Diff.Ratio"
      write(12,'(4(F12.5))') tflx, bflx, dflx, dflx / tflx

      write(12,*)
      write(12,*) "Photon Unit (umol/m2/s)"
      write(12,'(4(A12))') "Total","Beam","Diffuse","Diff.Ratio"
      write(12,'(4(F12.5))') tpfd, bpfd, dpfd, dpfd / tpfd

      write(12,*)
      write(12,*) "Albedo (A)"
      write(12,'(3(A12))') "Actual","black","white"
      write(12,'(3(A12))') " ","(beam)","(diffuse)"
      write(12,'(3(F12.5))') rflx / tflx, rbflx / bflx, rdflx / dflx

      if(stype .eq. 2) then

         write(12,*)
         write(12,*) "-- fraction of downward flux (T) -- "
         write(12,'(2(A15))') "Forest Floor","Surface floor"
         write(12,'(2(F15.8))') ffsum / tpfd, sfsum / tpfd

         write(12,*)
         write(12,*) "-- fraction of absorbed radiation (far) --"
         write(12,'(4(A10))') "total","leaf","non-leaf","floor"
         write(12,'(4(F10.7))') tfpr, cfpr, bfpr, ffpr

         write(12,*) " Vertical profile of far"
         write(12,'(A6,3(A10))') "H(m)","total","leaf","non-leaf"
         do k = int(zmax) + 1, 1, -1
            write(12,'(I6,3(F10.7))') k,app(k) + apnp(k),app(k),apnp(k)
         end do

      end if

      if(imode .ne. 1) then
         write(12,*)
         write(12,*) "-- spctrl fraction. energy & photon flux(%)"
         write(12,'(7(A10))') "wl(um)","Etot","Ebeam","Ediffuse",
     &        "Ptot","Pbeam","Pdiffuse"

         do i = 1, nwl
            if(i .le. 20) then
               write(12,'(7(F10.3))')
     &              wls + span(i) * real(i) - span(i) / 2.,
     &              (scmpf(j,i), j = 1, 3), (scmpp(j,i), j = 1,3)
            else
               write(12,'(7(F10.5))')
     &              0.7 + span(i) * real(i - 20) - span(i) / 2.,
     &              (scmpf(j,i), j = 1, 3), (scmpp(j,i), j = 1,3)
            end if
         end do
      end if

      write(12,*)
      write(12,*) "-- Downward radiance at the top of canopy --"
      write (12,*)
     &     "# idrc theta phi radiance_F stdev_F radiance_Q stdev_Q"

      pixnpa = real(np) / real(knxr * knyr)

      do k = 1, nrdc
         sumF  = 0.0
         sum2F = 0.0
         pFmax = 0.0
         sumQ  = 0.0
         sum2Q = 0.0
         pQmax = 0.0

         do j = 1, knyr
            do i = 1, knxr
               pF = real(prdcF(i, j, k)) / pixnpa
               pQ = real(prdcQ(i, j, k)) / pixnpa

               sumF  = sumF  + pF
               sum2F = sum2F + pF * pF
               pFmax = max(pFmax, pF)

               sumQ  = sumQ  + pQ
               sum2Q = sum2Q + pQ * pQ
               pQmax = max(pQmax, pQ)

            end do
         end do

         aveF = sumF / real(knxr * knyr)
         sdvF = sqrt(sum2F / real(knxr * knyr) - aveF * aveF)
         aveQ = sumQ / real(knxr*knyr)
         sdvQ = sqrt(sum2Q / real(knxr * knyr) - aveQ * aveQ)

         uxc = sign(max(1.e-5, abs(uxrtab(k))), uxrtab(k))
         phi = sqrt(uxrtab(k)**2 + uyrtab(k)**2)
         phi = max(1.e-5,phi)
         phi = acos(uxc / phi) / rad

         write (12,'(I6,2(F8.3),4(F12.5))')
     &        k,acos(uzrtab(k)) / rad,phi,aveF,sdvF,aveQ,sdvQ

      end do

      write(12,*)

      write(12,*) "-- PAR for tree, grass, stem and floor --"
      write(12,'(4(F15.3))') PPFD_T, PPFD_G, PPFD_S, PPFD_F

      close(12)

      if(stype .eq. 1) stop

!     BRDF results output
      if(cmode .eq. 1) then
         do i = 1, nangc
            brfc(1, i) = brfc(1, i) / brf(1, i)
            brfs(1, i) = brfs(1, i) / brf(1, i)
            brff(1, i) = brff(1, i) / brf(1, i)

            if(amode .eq. 1)then
               brf(1, i) = pi * brf(1, i) / (real(np) * tm)
            else
               brf(1, i) = pi * brf(1, i) / real(np)
            end if

            brfc(2, i) = brfc(2, i) / brf(2, i)
            brfs(2, i) = brfs(2, i) / brf(2, i)
            brff(2, i) = brff(2, i) / brf(2, i)
            brf(2, i) = pi * brf(2, i) / real(np)
         end do

         open(22,file="brfsum_"//trim(FILE_O_LABEL)//".txt")

         if(amode .eq. 1) then
            write(22,'(2(A8),5(A10))')
     &           " Theta", "Phi","BRF(TOC)","BRF(TOA)",
     &           "rover", "rbark", "rfloor"
            write(*,'(2(A8),5(A10))')
     &           " Theta", "Phi","BRF(TOC)","BRF(TOA)",
     &           "rover", "rbark", "rfloor"
         else
            write(22,'(2(A8),4(A10))')
     &           " Theta", "Phi","BRF(TOC)",
     &           "rover", "rbark", "rfloor"
            write(*,'(2(A8),4(A10))')
     &           " Theta", "Phi","BRF(TOC)",
     &           "rover", "rbark", "rfloor"
         end if

         k = 1
         do j=1, nph
            do i=1, nth

               if(amode .eq. 1)then
                  write(22,'(2(F8.2), 5(F10.6))')
     &                 angt(i), angp(j), brf(1, k),brf(2,k),
     &                 brfc(1, k),brfs(1, k),brff(1, k)
                  write(*,'(2(F8.2), 5(F10.6))')
     &                 angt(i), angp(j), brf(1, k),brf(2,k),
     &                 brfc(1, k),brfs(1, k),brff(1, k)
               else
                  write(22,'(2(F8.2), 4(F10.6))')
     &                 angt(i), angp(j), brf(1, k),
     &                 brfc(1, k),brfs(1, k),brff(1, k)
                  write(*,'(2(F8.2), 4(F10.6))')
     &                 angt(i), angp(j), brf(1, k),
     &                 brfc(1, k),brfs(1, k),brff(1, k)
               end if


               k = k + 1

            end do
         end do
         close(22)
      end if

!     Nadir view image

      if(cmode .ne. 1 .and. fmode .ne. 1)then
         open(22, file = "TOCIMG_"//trim(FILE_O_LABEL)//".txt")
         do j = size, 1, -1

            if(amode .eq. 1) then
               write(22, '(300(F10.5))')
     &              (pi*refl(1, i, j) / (pixnp * tm), i = 1, size)
            else
               write(22, '(300(F10.5))')
     &              (pi * refl(1, i, j) / pixnp, i = 1, size)
            end if

         end do
         close(22)
         
         open(23,  file = "nTOCIMG_"//trim(FILE_O_LABEL)//".txt")
         do j = size, 1, -1
            write(23, '(300(F10.5))')
     &           (real(irefl(1, i, j)) / pixnp, i = 1, size)
         end do
         close(23)
      end if

      if(cmode .eq. 3) then
         open(22, file = "apar.txt")
         open(23, file = "aparb.txt")
         open(24, file = "apard.txt")

         do k = int(zmax) + 1, 1, -1
            do j = int(size), 1, -1
               write(22,'(300(F12.5))')(ap(i, j, k), i = 1, size)
               write(23,'(300(F12.5))')(apb(i, j, k), i = 1, size)
               write(24,'(300(F12.5))')(apd(i, j, k), i = 1, size)
            end do
         end do
         close(22)
         close(23)
         close(24)

         open(22, file = "aparf.txt")
         open(24, file = "aparfd.txt")
         do j = int(size), 1, -1
            write(22,'(300(F12.5))') (apf(i, j), i = 1, size)
            write(24,'(300(F12.5))') (apfd(i, j), i = 1, size)
         end do
         close(22)
         close(24)

      end if

      ! SIF_results
      if (fmode .eq. 1) then
         !pixnpa =  pixnpa * np
         write(*,*) "cosq0 = ", cosq0
         write(*,*) "abs(cosq0) * RQ = ", RQ_COSQ0
         write(*,*) "abs(cosq0) * RF = ", RF_COSQ0
         write(*,*) "np = ", np
         write(*,*) "pixnpa = ", pixnpa
         write(*,*) "pixnp = ", pixnp
         write(*,*) "tm = ", tm
         call writeSIF(pixnpa, cosq0, np, tm)
      end if

      return

!     fish eye image

 200  open(24,file="feye_pl.txt")
      do i = 1, 90
         write(24,'(360(F12.5))') (feye(i,j), j = 1, 360)
      end do
      close(24)

      open(25,file="feye.txt")
      do i = 1, 180
         do j = 1, 180
            th = (real(i) - 90.)**2 + (real(j) - 90.)**2
            th = sqrt(th)
            ith = int(anint(th))
            if(ith .gt. 85) then
               eye(i, j) = 0.0
            else
               ph = (real(j) - 90.) / max(th,0.01)
               write(*,*) "i",i,j,th,ph
               ph = acos(ph) * 180. / pi
               iph = int(anint(ph))
               if(i .ge. 91) iph = 360 - iph
!               write(*,*) i,j,ith,iph,ph
               eye(i, j) = feye(ith, iph) / rfeye(ith, iph)
c       write(*,*) i,j,ith,iph,ph,(real(i) - 90.) /th
c     eye(i,j) = ph
            end if
         end do
         write(25,'(180(F12.5))') (eye(i,j), j = 1, 180)
      end do
      close(25)

!     calculate effective LAI by LAI-2000 algorithm
      W(1) = 0.034
      W(2) = 0.104
      W(3) = 0.160
      W(4) = 0.218
      W(5) = 0.484

!     S = 1 / cos(theta)

      S(1) = 1.008
      S(2) = 1.087
      S(3) = 1.270
      S(4) = 1.662
      S(5) = 2.670

      do i = 1,5
         T(i) = 0.0
      end do

!     make average over aziumth angle
      do i = 0, 90
         aaT(i) = 0.0
      end do
      do i = 1, 90
         do j = 0, 360
            aaT(i) = aaT(i) + feye(i, j) / rfeye(i, j)
         end do
         aaT(i) = aaT(i) / 361.
      end do

!     first ring
      ave = 0.0
      do i = 0, 12
         T(1) = T(1) + aat(i) * sin(real(i) * rad)
         ave = ave + sin(real(i) * rad)
      end do
      T(1) = T(1) / ave

!     second ring
      ave = 0.0
      do i = 17, 29
         T(2) = T(2) + aat(i) * sin(real(i) * rad)
         ave = ave + sin(real(i) * rad)
      end do
      T(2) = T(2) / ave

!     third ring
      ave = 0.0
      do i = 32, 43
         T(3) = T(3) + aat(i) * sin(real(i) * rad)
         ave = ave + sin(real(i) * rad)
      end do
      T(3) = T(3) / ave

!     forth ring
      ave = 0.0
      do i = 47, 58
         T(4) = T(4) + aat(i) * sin(real(i) * rad)
         ave = ave + sin(real(i) * rad)
      end do
      T(4) = T(4) / ave

!     fifth ring
      ave = 0.0
      do i = 62, 74
         T(5) = T(5) + aat(i) * sin(real(i) * rad)
         ave = ave + sin(real(i) * rad)
      end do
      T(5) = T(5) / ave

      elai = 0.0
      do i = 1,5
         write(*,*) "T= ", T(i)
         elai = elai + 2.0 * (W(i) / S(i)) * (-log(T(i)))
      end do
      write(*,*) "Effective LAI = ",elai

      return
      end

      ! write SIF results
      subroutine writeSIF(pixnpa, cosq0, np, tm)
      implicit none
      include "common.inc"
      include "math.inc"

      ! Input
      real pixnpa, RQ, cosq0, tm
      integer np

      ! work
      integer i, j, k, nm
      real cosBRF, const(100), const_radiation
      real const_observe
      character(5) sza_angle(30),tt

      ! vza string
      sza_angle(1) = '1 '
      sza_angle(2) = '2 '
      sza_angle(3) = '3 '
      sza_angle(4) = '4 '
      sza_angle(5) = '5 '
      sza_angle(6) = '6 '
      sza_angle(7) = '7 '
      sza_angle(8) = '8 '
      sza_angle(9) = '9 '
      sza_angle(10) = '10 '
      sza_angle(11) = '11 '
      sza_angle(12) = '12 '
      sza_angle(13) = '13 '
      sza_angle(14) = '14 '
      sza_angle(15) = '15 '
      sza_angle(16) = '16 '
      sza_angle(17) = '17 '
      sza_angle(18) = '18 '
      sza_angle(19) = '19 '
      sza_angle(20) = '20 '
      sza_angle(21) = '21 '
      sza_angle(22) = '22 '
      sza_angle(23) = '23 '
      sza_angle(24) = '24 '
      sza_angle(25) = '25 '
      sza_angle(26) = '26 '
      sza_angle(27) = '27 '
      sza_angle(28) = '28 '
      sza_angle(29) = '29 '
      sza_angle(30) = '30 '

      ! BRDF view
      cosBRF = 1
      areaTOC = real(size * size)

      ! debug
      !RF_COSQ0 = RQ_COSQ0
      
      write(*,*) "nangc=", nangc

      const_radiation = RQ_COSQ0 / areaTOC / pixnpa  * 1000.0
      const_observe =   RQ_COSQ0 / areaTOC / pixnpa  * 1000.0 

      write(*,*) "Start writing SIF results..."
      open(26, file = "SIF_results_"//trim(FILE_O_LABEL)//".txt")
      
      write(26, *) "View zenith angle: ", (angt(i), i = 1, nth)
      write(26, *) "View azimuth angle: ", (angp(i), i = 1, nph)

      write(26,'(1(A5), 120(A15))')
     &           "WL  ", 
     &           ("SIF.VZA."//trim(sza_angle(i)), i = 1, nangc),
     &           ("SIF.E.VZA."//trim(sza_angle(i)), i = 1, nangc),
     &           ("SIF.S.VZA."//trim(sza_angle(i)), i = 1, nangc),
     &           "SIF.ALL"

      ! i in [640, 850]
      do i = 1, 211

         write(26, *) i + 640 - 1,
     &   (Qfs(i,j) * const_observe, j = 1, nangc),
     &   (Qfs_e(i,j) * const_observe, j = 1, nangc),
     &   (Qfs_s(i,j) * const_observe, j = 1, nangc),
     &   Qfs_all(i)* const_radiation

      end do
      close(26)
      
      open(27, file = "SIF_PLACE_"//trim(FILE_O_LABEL)//".txt")
      ! i in [400, 750]
      do i = 1, 211
          write(27, *) i + 640 - 1, SIFplace(i,1) * const_observe,
     &                 SIFplace(i,2) * const_observe
      end do
      close(27)

      open(28, file = "SIF_PLACE_emit_"//trim(FILE_O_LABEL)//".txt")
      ! i in [400, 750]
      do i = 1, 211
          write(28, *) i + 640 - 1, SIFplace_e(i,1) * const_observe,
     &                 SIFplace_e(i,2) * const_observe
      end do
      close(28)

      nm = 760
      write(*,*) "SIF ANGLE: ", nm, "nm" 
      do i = 1, nangc
         write(*,*) "sza: ", angt(i)
         write(*,*) Qfs(nm - 640 + 1, i) * const_observe, 
     &              Qfs_all(nm - 640 + 1) * const_radiation
      end do

      !open(29, file = "Fluor_PAR_LEAF_"//trim(FILE_LABEL)//".txt")
      !do k = int(zmax) + 1, 1, -1
      !   do j = int(size), 1, -1
      !      write(29,'(300(F12.5))')(PAR_LEAF(i, j, k), i = 1, size)
      !   end do
      !end do
      !close(29)

      !open(30, file = "SIF_AREA_"//trim(FILE_LABEL)//".txt")
      !do j = 300, 1, -1
      !   write(30,'(300(F12.5))') (SIFarea(i, j), i = 1, 300)
      !end do
      !close(30)

      open(31, file = "SIF.emit.1st."//trim(FILE_O_LABEL)//".txt")
      do k = int(zmax) + 1, 1, -1
         do j = int(size), 1, -1
         write(31, *)(SIF_e_1st(i, j, k), i = 1, size)
         end do
      end do
      close(31)

      !open(32, file = "SIF.obs.loc."//trim(FILE_O_LABEL)//".txt")
      !do k = int(zmax) + 1, 1, -1
      !   do j = int(size), 1, -1
      !   write(32, *)(SIF_o_pos(i, j, k), i = 1, size)
      !   end do
      !end do
      !close(32)

  999 continue 
         
      write(*,*) "Photon = ", par0 + np
      write(*,*) "reflectance = ", par0
      write(*,*) "transmittance = ", par1
      write(*,*) "middle = ", par2
      write(*,*) "areaTOC = ", areaTOC
      write(*,*) "cosq0 = ", cosq0

      write(*,*) "C SIF760 = ", SIF_place_sca_num(760 - 640 + 1,1)
      write(*,*) "F SIF760 = ", SIF_place_sca_num(760 - 640 + 1,2)
      write(*,*) "sca C SIF760 = ", SIFplace(760 - 640 + 1,1) 
     &   * const_observe
      write(*,*) "sca F SIF760 = ", SIFplace(760 - 640 + 1,2)
     &   * const_observe
      write(*,*) "-- PAR for tree, grass, stem and floor --"
      write(*,'(4(F15.3))') PPFD_T, PPFD_G, PPFD_S, PPFD_F
      write(*,*) "-- W for tree, grass, stem and floor --"
      write(*,'(4(F15.3))') W_T, W_G, W_S, W_F
      write(*,*) "-- WF for tree, grass, stem and floor --"
      write(*,'(4(F15.3))') WF_T, WF_G, WF_S, WF_F
      end
