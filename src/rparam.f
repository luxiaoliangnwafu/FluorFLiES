!     Reading parameters
!**********************************************

      subroutine rparam(
     &     np,amode,cmode,imode,fmode,nwl,npl,nmix,
     &     atype,rtype,stype,cbnz,ctnz,cflg,
     &     dif,th0,ph0,phi,th,ph,
     &     sinf0,cosf0,cosq0,sinq0,
     &     wq,span,wl0,wls,RF,RQ,
     &     d,spcf,spcq,taur,
     &     ctaur,rfname,fname,alb,
     &     lr,lt,ulr,ult,str,sor,
     &     npproc, Nprod, num)

      implicit none
      include 'common.inc'
      include 'math.inc'

      integer knmix,knang,knzext
      parameter (knmix=10, knang=2000, knzext=200)

      integer*8 temp_np

      integer np,npproc
      integer Nprod,amode,imode,cmode,stype,nwl,iwl,fmode
      integer atype,rtype
      integer npl(1000),nmix
      integer cbnz,ctnz,cflg, num
      real dif,wq,span(1000),wl0,wls
      real th0,ph0,phi,th,ph
      real sinf0,cosf0,cosq0,sinq0
      real RF,RQ,d,spcf(1000),spcq(1000)
      real taur,ctaur,alb(100)

      real lr(5,1000), lt(5,1000), ulr(1000), ult(1000), str(5,1000),sor(1000)

      character*81 fname(knmix),rfname

!     total number of photon
      write(*,*) "np: Input # of photon"
      read(num,*) np
      !np = 1000

!     fish eye simulation mode
      if(np .eq. -4) then
         if(Nprod .ne. 1) then
            write(*,*) "fish eye mode cannot work with multi-processors"
            write(*,*) "please exucute under single processor"
            stop
         endif

         write(*,*) "fish eye mode - selected"
         call vegparam(np, nwl,cmode,fmode,wls,lr,lt,ulr,ult,
     &        str,sor,span(1),num)
         stype = 2
         return
      end if

!     LAI calculation mode
      if(np .eq. -5) then
         if(Nprod .ne. 1) then
            write(*, *) "LAI calcualtion mode cannot work",
     &           "with multi-prpcessor"
            write(*,*) "please exucute under single processor"
            stop
         endif

         write(*,*) "LAI calculation mode - selected"
         call vegparam(np, nwl,cmode,fmode,wls,lr,lt,ulr,ult,
     &        str,sor,span(1),num)
         stype = 2
         return
      end if

      npproc = int(real(np) / real(Nprod))

!     Atmospheric mode
      write(*,*) "amode: 1:atmospheric mode, 2: without atmospheric"
      read(num,*) amode
      !amode = 1

      if(amode.eq.2)then
         write(*,*) "dif: Input frac. of diffuse radiation (0-1)"
         read(num,*) dif
         if(dif.lt.0.0.or.dif.gt.1.0) then
            write(*,*) "fraction of diffuse should be 0.0-1.0... exit "
            stop
         end if
      end if

!     Get geometrical parameters
      call atmgeo
     &     (th0, ph0, sinq0, cosq0, sinf0, cosf0, amode, num)

!     get atmospheric parameters
      call atmparam
     &     (amode, imode, wq, span(1), wl0, wls,
     &     nwl, RF, RQ, npl, np, npproc, Nprod, nmix, atype,
     &     rtype, rfname, fname, d, spcf, spcq,
     &     taur, ctaur, cbnz, ctnz, cflg, num, cosq0)


!     vegetation parameters read/ initialize

      write(*,*) "stye: Surface mode 1: Lambertian, 2: 3-D Vegetation"
      read(num,*) stype
      !stype = 2
      if(stype .eq. 1) then
         do iwl = 1, nwl
            if(imode .eq. 1) then
               write(*,*) "Input albedo"
            else
               if(iwl .le. 20) then
                  write(*, "('Input albedo ',F5.3, ' and ' ,F5.3)")
     &                 wls + span(iwl) * real(iwl - 1),
     &                 wls + span(iwl) * real(iwl)
               else
                  write(*, "('Input albedo ',F5.3, ' and ' ,F5.3)")
     &                 0.7 + span(iwl) * real(iwl - 21),
     &                 0.7 + span(iwl) * real(iwl - 20)
               end if
            end if
            read(num,*) alb(iwl)
         end do
      elseif(stype .eq. 2)then

         call vegparam(np, nwl,cmode,fmode,wls,lr,lt,ulr,ult
     &        ,str,sor,span(1),num)
c         call paraminit(nwl,lr,lt,nts) <= modifir!!!!
      else
         write(*, *) "Bad surface type selection exit"
      end if

      end

!*************************************************
      subroutine atmgeo
     &     (the0, phi0, sinq0, cosq0, sinf0, cosf0, mode, num)
!*************************************************

!     implicit none

      include 'common.inc'
      include 'math.inc'

      integer i,j,k,ntha,npha,mode, num
      real the0, phi0, sinq0, cosq0, sinf0, cosf0, th(18), ph(36)
      real f0,q0

!     solar zenith angle
      write(*,*) "the0: Solar zenith angle (degree) & Elevation (m)"
      read(num,*) the0, zmin
      !the0 = 30
      !zmin = 0
      
      ! solar azimuth angle
      write(*,*) "the0: Solar azimuth angle (degree)"
      read(num,*) phi0

      q0 = pi - rad * the0
      if(phi0 .le. pi)then
         f0 = rad * phi0 + pi
      else
         f0 = rad * phi0 - pi
      end if

      sinq0 = sin(q0)
      cosq0 = cos(q0)
      sinf0 = sin(f0)
      cosf0 = cos(f0)

      if(mode .eq. 2)return

!     radiance sampling angle
      write(*,*) "Radiance at the bottom atmospheric boundary"
      write(*,*) "ntha, th: # of angle and zenith angle(degree,max 18)"
      write(*,*) "ex. 5 100. 120. 140. 160. 170."
      read(num,*) ntha,(th(i),i = 1, ntha)
      !ntha = 1
      !th(1) = 100
      write(*,*) "npha, ph: # of angle and azimuth angle(degree,max 36)"
      write(*,*) "ex. 3 0. 90. 180."
      read(num,*) npha,(ph(i),i = 1, npha)
      !npha = 1
      !ph(1) = 0

      k = 0
      do i = 1, npha
         do j = 1, ntha
            if(th(j) .lt. 91.) then
               write(*,*) "Zenith angle should be greater than 91"
               write(*,*) th(j)," is ignored !"
               write(*,*) "    "
            else
               k = k + 1
               uxrtab(k) = sin(rad * th(j)) * cos(rad * ph(i))
               uyrtab(k) = sin(rad * th(j)) * sin(rad * ph(i))
               uzrtab(k) = cos(rad * th(j))
            end if
         end do
      end do
      nrdc = k                    ! # of radiances

      return
      end

!**************************************************
      subroutine atmparam
     &     (amode,imode, wq, span, wl0, wls,
     &     nwl, RF, RQ, npl, np, npproc,Nprod, nmix, atype,
     &     rtype, rfname, fname, d, spcf, spcq,
     &     tauref, ctauref, cbnz, ctnz, cflg, num, cosq0)
!*************************************************
      implicit none
      include 'common.inc'
      include 'math.inc'

      integer i,j,k,amode,imode, cflg, nmix, nwl, num
      integer ctype, atype, rtype, ctnz, cbnz
      integer nspc, np, npl(1000), npproc, Nprod

      real d, wq, span, wls, wl0
      real ctauref,tauref, cosq0
      real cbot, ctop,dum, wlsd, wl0d(4000)
      real spcf(1000),spcq(1000),spcdf(4000),spcdq(4000)
      real parF,parQ,swF,swQ,RF,RQ, iRF, iRQ
      parameter(parF=531.2593,parQ=2423.93994,
     &     swF=1351.81531,swQ=10213.1367)

      character*81 fname(10),rfname
      character*2 ch(3)

      ch(1) = "hi"
      ch(2) = "lo"
      ch(3) = "lo"

      nspc = 0
      do i = 1, 100
         spcf(i) = 0.0
         spcq(i) = 0.0
         npl(i) = 0
      end do

!     Integration mode
      if(amode .eq. 2)then
!     write(*,*)"Only monochro wavelength calculation!"
         imode = 1
!     write(*,*) "Irradiance, photon flux at the top of atmosphere"
!     read(num,*) RF,RQ
!     these are the fault and tentative value and no physical meaning
         RF = 1000.
         RQ = 1000.
         wq = 1.0
         nwl = 1
      else
         write(*,*) "imode: Integration mode 1:Monochro 2:PAR ",
     &              "3:SW (SIF PART)"
         read(num,*) imode

         ! input PAR ffqd, then convert to PAR
         write(*,*) "PAR photon flux (PPFD):"
         read(num, *) iRQ
         
         ! 1 W/m2 ≈ 4.6 μmole.m2/s
         ! RF is W/m2
         ! RQ is μmole.m2/s
         iRF = iRQ / 4.57

         write(*,*) "PAR PPFD = ", iRQ
         write(*,*) "PAR irradiance = ", iRF

         !imode = 3
         if(imode .lt. 0 .or. imode .gt. 4)then
            write(*,*) "Mode error!",imode
            stop
         elseif(imode .eq. 1)then
            write(*,*) "wl0:wavelength (micron ex 0.55)"
            read(num,*) wl0
            wls = 0.2
            nwl = int(anint((4.0 - 0.3) / span))
            npl(1) = npproc
            write(*,*) "np",np,npl(1)
         elseif(imode .eq. 2)then
            wls = 0.4
            nwl = int(anint((0.7 - 0.4) / span))
            RF = parF
            RQ = parQ
         elseif(imode .eq. 3)then
            wls = 0.4
            WL_START = 0.4
            WL_END = 0.86
            ! SIF 0.3 - 1.5
            nwl = int((WL_END - WL_START) / span)
            !nwl = nwl + int((4.0 - 0.7) / (5. * span))

            ! because SIF has high coorelation in PAR,
            ! so set par incidiance at TOA, TOC
            RF = iRF
            RQ = iRQ
         end if

!     read solar radiation
         open(12, file = "Data/solar_rad")
         read(12, *)
         do i = 1, 3802
            read(12, *) wl0d(i), spcdf(i), spcdq(i), dum
            nspc = nspc + 1
         end do
         close(12)

!     search a spectral irradiance data in monochromatic calculation
         j = 1
         wlsd = wls + 0.0005

         if(imode .eq. 1)then
            do i = 1, nspc
               if(wlsd .gt. (wl0 - 0.0025)
     &              .and. wlsd .lt. (wl0 + 0.0025))then
                  RF = (wl0d(i + 1) * spcdf(i) + wl0d(i) * spcdf(i + 1))
     &                 / (wl0d(i) + wl0d(i + 1))
                  RQ = (wl0d(i + 1) * spcdq(i) + wl0d(i) * spcdq(i + 1))
     &                 / (wl0d(i) + wl0d(i + 1))
                  spcf(j) = RF
                  spcq(j) = RQ
                  wq = 1.0
                  nwl = 1
                  exit
               end if
               wlsd = wlsd + 0.001
            end do
         end if

         if(imode .eq. 2)then
            do i = 1, nspc
               if(abs(wl0d(i) - wlsd) .lt. 1.0d-04) then
                  do k = 0, 19
                     spcf(j) = spcf(j) + spcdf(i + k)
                     spcq(j) = spcq(j) + spcdq(i + k)
                  end do

                  spcf(j) = spcf(j) / (span * 1000.)
                  spcq(j) = spcq(j) / (span * 1000.)

                  j = j + 1
                  if(j .gt. nwl) exit
                  wlsd = wlsd + span
               end if
            end do
         end if

         if(imode .eq. 3)then
            !do i = 1, nspc
            !   if(abs(wl0d(i) - wlsd) .lt. 1.0d-04) then
            !      if(wlsd .lt. 0.7)then
            !         do k = 0, 19
            !            spcf(j) = spcf(j) + spcdf(i+k)
            !            spcq(j) = spcq(j) + spcdq(i+k)
            !         end do

            !         spcf(j) = spcf(j) / (span * 1000.)
            !         spcq(j) = spcq(j) / (span * 1000.)

            !      else
            !         do k = 0, 99
            !            spcf(j) = spcf(j) + spcdf(i+k)
            !            spcq(j) = spcq(j) + spcdq(i+k)
            !         end do

            !         spcf(j) = spcf(j) / (span * 5000.)
            !         spcq(j) = spcq(j) / (span * 5000.)

            !      end if

            !      j = j + 1
            !      if(j .gt. nwl) exit
            !      if(wlsd .lt. 0.7)then
            !         wlsd = wlsd + span
            !      else
            !         wlsd = wlsd + 5. * span
            !      end if
            !   end if
            !end do

            ! modified for SIF [400 - 1200]
            do i = 1, nspc
              if (wl0d(i) .ge. wlsd) then
                 spcf(j) = spcdf(i)
                 spcq(j) = spcdq(i)

                 j = j + 1
                 if (j .gt. nwl) exit
                 wlsd = wlsd + 0.001
              end if
            end do
            write(*,*) "cosq0 = ", cosq0
            write(*,*) "RF = ", RF
            write(*,*) "RQ = ", RQ
            RF_COSQ0 = abs(cosq0) * RF
            RQ_COSQ0 = abs(cosq0) * RQ
            write(*,*) "abs(cosq0) * RF = ", abs(cosq0) * RF
            write(*,*) "abs(cosq0) * RQ = ", abs(cosq0) * RQ

         end if

!     make input photon wieght for each wavelength
         dum = 0.0
c         write(*,*) nwl

         do i = 1, nwl
            !if(i .le. 20)then
            !   dum = dum + spcf(i)
            !else
            !   dum = dum + 5. * spcf(i)
            !end if
            dum = dum + spcf(i)
         end do
         do i = 1, nwl
            !if(i .le. 20)then
            !   npl(i) = anint(real(npproc) * spcf(i) / dum)
            !else
            !   npl(i) = anint(real(npproc) * 5. * spcf(i) / dum)
            !end if
            npl(i) = anint(real(npproc) * spcf(i) / dum)
         end do

         dum = 0.

         do i = 1, nwl
            dum = dum + npl(i)
         end do

         npproc = int(dum)            ! Modification of np for actual value
         np = npproc * Nprod
         write(*,*) "Actural number of photon is: ", np

      end if                    ! end if for amode

!     with/without atmospheric
      if(amode .eq. 2)then
         rtype = 0
         atype = 0
         tauref = 0.0
         ctype = 0
      else

!     read z profile
         open(14, file="Data/zgrd")
         read(14,*) (zgrd(i), i = 0, nz)
         close(14)

!     rescaling of zgrd according to the elevation
         do i = 0, nz
            zgrd_back(i) = zgrd(i)
            zgrd(i) = (zgrd(nz) - zmin) * (zgrd(i) / zgrd(nz)) + zmin
            klayer(i) = i
!            write(*,*) "zgrd, klayer",zgrd(i), klayer(i)
         end do

!     mapping the iz to the actual height level corrected above
         do i = 0, nz
            do j = 0, nz
               if(zgrd(i) .lt. zgrd_back(j)) then
                  klayer(i) = j - 1
!                  write(*,*) "zgrd zgrd_back klayer",
!     &                 zgrd(i), zgrd_back(i), klayer(i)
                  exit
               end if
            end do
         end do

!     make the middle of the layer height of the original zgrd (zgrd_back)
         do i = 1, nz
            zgrdm(i) = 0.5 * (zgrd_back(i) + zgrd_back(i-1))
            write(*,*) "zgrdm",zgrdm(i)
         end do

         write(*,*) "rtype: Atmospheric profile"
         write(*,*) " 1: Tropical"
         write(*,*) " 2: Mid latitude summer"
         write(*,*) " 3: Mid latitude winter"
         write(*,*) " 4: High latitude summer"
         write(*,*) " 5: High latitude winter"
         write(*,*) " 6: US standard atm."
         read(num,*) rtype
         !rtype = 1

         if(rtype .eq. 1)then
            rfname = "Data/gas_TR_"//ch(imode)
         elseif(rtype .eq. 2)then
            rfname = "Data/gas_MS_"//ch(imode)
         elseif(rtype .eq. 3)then
            rfname = "Data/gas_MW_"//ch(imode)
         elseif(rtype .eq. 4)then
            rfname = "Data/gas_HS_"//ch(imode)
         elseif(rtype .eq. 5)then
            rfname = "Data/gas_HW_"//ch(imode)
         elseif(rtype .eq. 6)then
            rfname = "Data/gas_US_"//ch(imode)
         else
            write(*,*) "bad selection "
            stop
         endif

!     read the aerosol data
         write(*,*) "atype: aerosol type"
         write(*,*) " 1:  Continental clean"
         write(*,*) " 2:  Continental average"
         write(*,*) " 3:  Continental polluted"
         write(*,*) " 4:  Urban"
         write(*,*) " 5:  Desert"
         write(*,*) " 6:  Maritime clean"
         write(*,*) " 7:  Maritime polluted"
         write(*,*) " 8:  Maritime Tropical"
         write(*,*) " 9:  Arctic"
         write(*,*) "10:  Antactic"
         write(*,*) "11:  Smoke"
         read(num,*) atype
         !atype = 1

!     AOT
         write(*,*) "tauref: AOT at 0.550 micron"
         read(num,*) tauref
         !tauref = 0.5
         write(*,*) "                "
         write(*,*)" - this version uses a extinction with RH=70% -"

!     current version uses a mixed extinction coef. by Iwabushi san
         nmix=2

         if(atype .eq. 1)then
            fname(1) = "Data/opt_type1_rh0.70_"//ch(imode)
            d = 8000.0            !     scale height (m)
         elseif(atype .eq. 2)then
            fname(1) = "Data/opt_type2_rh0.70_"//ch(imode)
            d = 8000.0            !     scale height (m)
         elseif(atype .eq. 3)then
            fname(1) = "Data/opt_type3_rh0.70_"//ch(imode)
            d = 8000.0            !     scale height (m)
         elseif(atype .eq. 4)then
            fname(1) = "Data/opt_type4_rh0.70_"//ch(imode)
            d = 8000.0            !     scale height (m)
         elseif(atype .eq. 5)then
            fname(1) = "Data/opt_type5_rh0.70_"//ch(imode)
            d = 2000.0            !     scale height (m)
         elseif(atype .eq. 6 .or. atype .eq. 8)then
            fname(1) = "Data/opt_type6_rh0.70_"//ch(imode)
            d = 1000.0            !     scale height (m)
         elseif(atype .eq. 7)then
            fname(1) = "Data/opt_type7_rh0.70_"//ch(imode)
            d = 1000.0            !     scale height (m)
         elseif(atype .eq. 8)then
            fname(1) = "Data/opt_type8_rh0.70_"//ch(imode)
            d = 1000.0            !     scale height (m)
         elseif(atype .eq. 9)then
            fname(1) = "Data/opt_type9_rh0.70_"//ch(imode)
            d = 99000.0           !     scale height (m)
         elseif(atype .eq. 10)then
            fname(1) = "Data/opt_type10_rh0.70_"//ch(imode)
            d = 99000.0           !     scale height (m)
         elseif(atype .eq. 11)then
            write(*,*) "smoke aerosol under construction !!"
            stop
            d = 8000.0            !     scale height (m) (tentative value)
         else
            write(*,*) "Bad aerosol selection!"
            stop
         end if

!     read cloud data
         write(*,*) "ctype: cloud type"
         write(*,*) " 0:  Cloud-free"
         write(*,*) " 1:  Stratus continental"
         write(*,*) " 2:  Stratus maritime"
         write(*,*) " 3:  Cumulus continental clean"
         write(*,*) " 4:  Culumus continental pulluted"
         write(*,*) " 5:  Culumus maritime"
         write(*,*) " 6:  Fog"
         write(*,*) " 7:  Cirrus 1 (-25degC)"
         write(*,*) " 8:  Cirrus 2 (-50 degC)"
         write(*,*) " 9:  Cirrus 3 (-50 degC + small particles"

         read(num,*) ctype
         !ctype = 0

         if(ctype .ne. 0)then
            write(*,*) "ctauref:COT at 0.55 micron"
            read(num,*) ctauref
            write(*,*) " ctop, cbot: cloud top and bottom height (m)"
            read(num,*) ctop,cbot

            cflg = 1

            if(ctop .lt. cbot) then
               write(*,*) "cloud top should be greater",
     &              "than cloud bottom"
               stop
            endif

            if(ctype .eq. 1)then
               fname(nmix) = "Data/opt_type101_"//ch(imode)
            elseif(ctype .eq. 2)then
               fname(nmix) = "Data/opt_type102_"//ch(imode)
            elseif(ctype .eq. 3)then
               fname(nmix) = "Data/opt_type103_"//ch(imode)
            elseif(ctype .eq. 4)then
               fname(nmix) = "Data/opt_type104_"//ch(imode)
            elseif(ctype .eq. 5)then
               fname(nmix) = "Data/opt_type105_"//ch(imode)
            elseif(ctype .eq. 6)then
               fname(nmix) = "Data/opt_type106_"//ch(imode)
            elseif(ctype .eq. 7)then
               fname(nmix) = "Data/opt_type107_"//ch(imode)
            elseif(ctype .eq. 8)then
               fname(nmix) = "Data/opt_type108_"//ch(imode)
            elseif(ctype .eq. 9)then
               fname(nmix) = "Data/opt_type109_"//ch(imode)
            else
               write(*,*) "bad cloud selection"
               stop
            end if

            do i = 0, nz
               if(zgrd(i) .lt. cbot)then
                  cbnz = i
               end if
               if(zgrd(i) .gt. ctop)then
                  ctnz = i
                  exit
               end if
            end do

            write(*,*) cbnz,ctnz
            write(*,*) "clouds are located between",
     &           zgrd(ctnz), zgrd(cbnz)
         end if                 ! if for amode
      end if

      return
      end

!**********************************************
!     Reading vegetation parameters
!     & checking thier defintion and ranges
!     by H. Kobayashi
!     last modified 05/05/26
!
!**********************************************
      subroutine vegparam(np, nwl,cmode,fmode,wls,lr,lt,ulr,ult,
     &     str,sor,span, num)

      implicit none

      include 'common.inc'
      include 'math.inc'

      integer kobj, l, np, num, ispc

      real w,w1(5),w2,w0,umax
      real tp(1000),ctemp,wls,span
      integer i,j,k,m,Obj_nt,nwl,cmode,fmode
      real lr(5,*),lt(5,*),ulr(*),ult(*),str(5,*),sor(*)
      real ramda,avelr,avelt
      real max, min, p1, p2, a(4), b(4), c(4), cc(4)
      real xr, yr, d, dd, rx, ry, rr
      character*81 crowndata_file
      real i_ulr_p, i_ulr_n, i_ult_p, i_ult_n
      real i_sr_p, i_sr_n, i_str_p, i_str_n
      real crown_scale
      integer leaf_incline

!     read the condition file
!      write(*,*) "Input boundary condition 1:periodic 2:non-periodic"
!      read(num,*) bound
      bound = 1

      if((bound .ne. 1) .and. (bound .ne. 2))then
         write(*,*) "boundary condition should be 1 or 2."
         write(*,*) "  1: Periodic"
         write(*,*) "  2: Non-periodic"
         stop
      end if

      if(np .eq. -4 .or. np .eq. -5) goto 201

!     input output mode
      write(*,*) "cmode: calculation mode"
      write(*,*) "1: BRF only ", "2: BRF Nadir Image ",
     &     "3: 3D APAR"
      read (num,*) cmode
      !cmode = 2

      if(cmode .le. 0 .or. cmode .gt. 4) then
         write(*,*) "Bad mode selection exit"
         stop
      end if

      if(cmode .ne. 1) goto 100

      write(*,*) "nth, angt: # of angle and zenith angle(max 18)for BRF"
      write(*,*) "eg. 5 10. 20. 45. 50. 70."
      read(num,*) nth,(angt(i),i = 1, nth)
      !nth = 1
      !angt(1) = 30.0
      write(*,*)  nth,(angt(i),i = 1, nth)
      write(*,*) "nph, angp:# of angle and azimuth angle(max 36)for BRF"
      write(*,*) "eg. 3 0. 90. 180."
      read(num,*) nph,(angp(i),i = 1, nph)
      !nph = 1
      !angp(i) = 0.0

      if(nth * nph .gt. 100000)then
         write(*,*) "The number of sampling angles are too huge!"
         write(*,*) "should be theta*phi<348"
         stop
      end if

      if(nth .gt. 100000)then
         write(*,*) "The number of sampling theta are too huge!"
         write(*,*) "should be theta*phi<100"
         stop
      end if

      if(nph .gt. 100000)then
         write(*,*) "The number of sampling phi are too huge!"
         write(*,*) "should be theta*phi<100"
         stop
      end if

      ! SIF
      if (nth .ne. nph) then
         write(*,*) "vza and vaa are not equal!"
         write(*,*) nth, nph
         stop
      end if


      k = 0
      !do i = 1, nph
      !   do j = 1, nth
      !      if(angt(j) .gt. 80.) then
      !         write(*,*) "Zenith angle should be less than 80"
      !         write(*,*) angt(j)," is ignored !"
      !         write(*,*) "    "
      !         nth = nth - 1

      !      else
      !         k = k + 1
      !         uxrc(k) = sin(rad * angt(j)) * cos(rad * angp(i))
      !         uyrc(k) = sin(rad * angt(j)) * sin(rad * angp(i))
      !         uzrc(k) = cos(rad * angt(j))
      !         write(*,*) "uxrc, uyrc, uzrc:"
      !         write(*,*) uxrc(k), uyrc(k), uzrc(k)
      !      end if
      !   end do
      !end do

      ! SIF
      do i = 1, nph
         if (angt(i) .gt. 80) then 
            write(*,*) "Zenith angle should be less than 80"
            write(*,*) angt(j)," is ignored !"
         else
            k = k + 1
            uxrc(k) = sin(rad * angt(i)) * cos(rad * angp(i))
            uyrc(k) = sin(rad * angt(i)) * sin(rad * angp(i))
            uzrc(k) = cos(rad * angt(i))
            write(*,*) "uxrc, uyrc, uzrc:"
            write(*,*) uxrc(k), uyrc(k), uzrc(k)
         end if
      end do
      nangc = k                    ! # of radiances

 100  continue
      if(cmode .ge. 2) then
         nangc = 1
         uxrc(1) = 0.0
         uyrc(1) = 0.0
         uzrc(1) = 1.0
      end if


!     get the number of forest species
 201  write(*,*) "nts: # of group of tree species"
      read(num,*) nts
      !nts = 1
      if(np .eq. -4 .or. np .eq. -5) goto 202

      u(6) = real(nts)

!     currently max nts should be less than 5
      if(nts .le. 0 .or. nts .gt. 5)then
         write(*,*) "Error # of tree species"
         write(*,*) "tree species should be 1-5"
         stop
      end if

!     read leaf reflectance/transmittance
      write(*,*) "underground reflectance par and nir"
      read(num,*) i_ulr_p, i_ulr_n

      write(*,*) "underground transmittance par and nir"
      read(num,*) i_ult_p, i_ult_n

      write(*,*) "soil reflectance par and nir"
      read(num,*) i_sr_p, i_sr_n
      
      write(*,*) "stem reflectance par and nir"
      read(num,*) i_str_p, i_str_n

      write(*,*) "Simulation label.:"
      read(*,('(A)'))  FILE_LABEL

      write(*,*) "Simulation output label.:"
      read(*,('(A)'))  FILE_O_LABEL

      ramda = wls
      if(np .eq. -4) nwl = 1


      call readspectrum(lr, lt)
      do i = 1, nwl

         !if(nwl .eq. 1)then
            !write(*,*) "Input surface optical properties"
         !else
           !write(*, *) "Input surface optical properties in"
          !  if(i .le. 2000)then
              !write(*,"(F5.3,' and ',F5.3)")
       !&              ramda, ramda + span
           ! else
              !write(*,"(F5.3,' and ',F5.3)")
        !&              ramda, ramda + 5. * span
            !end if
         !end if

         !write(*,*) "lr1 lr2.. lt1 lt2.. ulr ult str1 str2.. sor"
         !write(*,*) "(lr, lt:leaf refl. & leaf transm."
         !write(*,*) "ulr,ult:understory leaf refl. & transm."
         !write(*,*) "stmr: stem refl., soilr: soil refl.)"

         if(nwl .eq. 1) then
            read(num,*) (lr(j,i), j = 1, nts),(lt(j,i), j = 1, nts),
     &           ulr(i), ult(i),(str(j,i), j = 1, nts), sor(i)
         else
            if(i .eq. 1) then
              !write(*,*) "Input PAR average values"

               !read(num,*) (lr(j,i), j = 1, nts),(lt(j,i), j = 1, nts), ulr(i), ult(i),(str(j,i), j = 1, nts), sor(i)

               !lr(1,1) = 0.1539255
               !lt(1,1) = 0.1469605
               ulr(1) = i_ulr_p
               ult(1) = i_ult_p
               str(1,1) = i_str_p
               sor(1) = i_sr_p
               ispc = i
            !elseif(i .eq. 21) then
            elseif (i + 400 .eq. 801) then
               !write(*,*) "Input NIR average values"
      !        read(num,*) (lr(j,i), j = 1, nts),(lt(j,i), j = 1, nts), ulr(i), ult(i),(str(j,i), j = 1, nts), sor(i)
               !lr(1,i) = 0.4492953
               !lt(1,i) = 0.4881736
               ulr(i) = i_ulr_n
               ult(i) = i_ult_n
               str(1,i) = i_str_n
               sor(i) = i_sr_n

               ispc = i
            else
               do j = 1, nts
                  !lr(j, i) = lr(j, ispc)
                  !lt(j, i) = lt(j, ispc)
                  str(j, i) = str(j, ispc)
               end do
               ulr(i) = ult(ispc)
               ult(i) = ult(ispc)
               sor(i) = sor(ispc)
               !write(*,*)  (lr(j,i), j = 1, nts),(lt(j,i), j = 1, nts),
              !   ulr(i), ult(i),(str(j,i), j = 1, nts), sor(i)
            end if
         end if

!     check the parameter range
         do j = 1, nts
            write(*,*) i, lr(j,i), lt(j,i)
            if((lr(j, i) + lt(j, i)) .gt. 1)then
               write(*,*) "canopy leaf reflectance+transmittance"
               write(*,*) "      is too large, exit!"
               stop
            end if

            if((lr(j, i) + lt(j, i)) .lt. 0.0001)then
               lr(j, i) = 0.0001
               lt(j, i) = 0.0001
            end if

         end do


         if((ulr(i) + ult(i)) .gt. 0.99)then
            write(*,*) "floor leaf reflectance+transmittance"
            write(*,*) "      is too large, exit!"
            stop
         end if

         if((ulr(i) + ult(i)) .lt. 0.0001)then
            ulr(i) = 0.0001
            ult(i) = 0.0001
         end if

         do j = 1, nts
            if(str(j,i) .gt. 1.00)then
               write(*, *) "stem reflectance"
               write(*, *) "      is too large, exit!"
               stop
            end if

            if(str(j, i) .lt. 0.0001)then
               str(j, i) = 0.0001
            end if

         end do

         if(sor(i) .gt. 1.00)then
            write(*,*) "soil reflectance"
            write(*,*) "      is too large, exit!"
            stop
         end if

         !if(nwl .le. 20)then
          !  ramda = ramda + span
         !else
          !  ramda = ramda + 5. * span
         !end if
         ramda = ramda + span
      end do

 202  write(*,*) "u: leaf area density 1,2,3...# tree species"
      read(num,*) (u(i),i = 1, nts)
      !u(1) = 0.5

      if(.not.(np .eq. -4 .or. np .eq. -5)) then
         write(*,*) "gLAI: forest floor LAI"
         read(num,*) gLAI
         !gLAI = 0.5
      end if

      write(*,*) "BAD: branch area density 1,2,3... # of tree species"
      read(num,*) (BAD(i), i = 1, nts)
      !BAD(1) = 0.5

      do i = 1, nts
         if(u(i) .lt. 0.0 .or. u(i) .gt. 8.0)then
            write(*,*) i,"th leaf area density",u(i),
     &           " should be set in the range (0.0-8.0)"
            write(*,*) "   exit!"
            stop
         end if
      end do

      if(gLAI .lt. 0.0 .or. gLAI .gt. 8.0)then
         write(*,*) gLAI, "floor LAI",
     &       " should be set in the range (0.0-8.0)"
         write(*,*) "   exit!"
         stop
      end if

      if(np .ne. -5) then
         write(*,*) "sbar: Spherical ave. shoot silhouette to ",
     &        "totall needle area ratio"
         write(*,*) "1,2,3... # of tree species (0.0-0.25)"
         write(*,*) "For broadleaves, please input 0.25"
         read(num,*) (sbar(i), i = 1, nts)
         !sbar(1) = 0.25

         !SIF
         write(*,*) "fmode: active SIF simulation? 1: on; 0: off"
         read(num,*)  fmode
         !fmode = 1
      end if 

      write(*,*) "Crowndata filename:"
      read(*,('(A)'))  crowndata_file
      
      if(np .ne. -5) then

         write(*,*) "leaf inclination"
         write(*,*) "(1=spherical 2 = planophile 3 = erectrophile):"
         read(num,*) leaf_incline

         write(*,*) "Crown shape scale factor:"
         read(num,*)  crown_scale

         write(*,*) "Grass layer exist (0: non, 1:exist):"
         read(num,*)  GRASS

         if (GRASS .eq. 1) then
            write(*,*) "Grass layer EF-matrix folder name:"
            read(*,('(A)'))  FILE_LABEL_G

            write(*,*) "canopy and grass use "
            write(*,*) "the same spectrum (0: no, 1:yes):"
            read(num,*)  C_G_SPECTRUM
         endif

         

      end if

      mc = leaf_incline
      mb = leaf_incline
      mf = leaf_incline

      !write(*,*) "fluorescence quantum efficiency for PSI, PSII: "
      !read(num,*) FQEI, FQEII
      !FQEI = 0.004
      !FQEII = 0.02

      umax = u(1)
      do i = 2, nts
         if(umax .lt. u(i))umax = u(i)
      end do

      if(umax .gt. 1.0)then
         fe = 1.0
      elseif(umax .le. 0.01)then
         fe = 0.01
      else
         fe = umax
      end if
      !fe = 1
      write(*,*) "fe=", fe

      if (fmode .eq. 1)then
      ! read EF matrices
         call readEFMatrix(MBI, 0, 1, FILE_LABEL)
         call readEFMatrix(MBII, 0, 2, FILE_LABEL)
         call readEFMatrix(MFI, 1, 1, FILE_LABEL)
         call readEFMatrix(MFII, 1, 2, FILE_LABEL)

         if (GRASS .eq. 1) then
             call readEFMatrix(gMBI, 0, 1, FILE_LABEL_G)
             call readEFMatrix(gMBII, 0, 2, FILE_LABEL_G)
             call readEFMatrix(gMFI, 1, 1, FILE_LABEL_G)
             call readEFMatrix(gMFII, 1, 2, FILE_LABEL_G)
         endif
      end if

!     canopy object parameters
!     object id initialization
      do i = 1, 6000
         iobj(i) = 1
      enddo

      areaTOC = 0.0
      
      open(23, file="Data/crowndata/"//trim(crowndata_file))
      read(23, *) nobj

!     input obj_nt
      obj_nt = nobj

!     if obj_num=0, dummy data are added
!     to prevent from loop counter error
      if(nobj .eq. 0)then
         nobj = 1
         sobj(1) = 1
         obj(1,1) = 0.01
         obj(1,2) = 0.01
         obj(1,3) = 0.01
         obj(1,4) = 1.d-5
         obj(1,5) = 1.d-5
      else
         k = 1
         do i = 1, nobj
            read(23,*) sobj(k), obj(k,1), obj(k,2)
     &           , obj(k,3), obj(k,4), obj(k,5)
     &           , iobj(k)
            
            ! scale crown shape
            obj(k,5) = obj(k,5) * crown_scale
            !if (sobj(k) .eq. 3) then
                !obj(k,5) = obj(k,5) * crown_scale
            !end if

!     if the canopy object is less that 0.1 (m), this object is neglected
!     to get the efficient computer time.

            if(sobj(k) .ne. 4)then
               if((obj(k,4) .lt. 0.01) .or. (obj(k,5) .lt. 0.01))then
                  write(*,*) k, "th canopy neglected !"
                  Obj_nt = Obj_nt - 1
               end if
            end if
            k = k + 1
         end do
      end if
      close(23)

      nobj = Obj_nt


!     check the objid range
      do i = 1, nobj
         if(iobj(i) .le. 0 .or. iobj(i) .gt. nts)then
            write(*,*) "species id should be the range betweem 0-"
     &           ,nts
            exit
         end if
      end do

      write(*,*) "Total Object is ",nobj

      do i = 1, nobj
         if(sobj(i) .eq. 3)then
!     change from height to radius
            obj(i,4) = obj(i,4) / 2.
         end if

         ! count TOC area
         if (sobj(i) .ne. 4) then
            areaTOC = areaTOC + obj(i, 5) * obj(i, 5) * pi * 4
         end if
      end do
      write(*,*) "TOC area: ", areaTOC

!     in case periodic boudary
!     add objects that are partially in the outside from the simulation scence
      if(bound .eq. 1) then

!     set distance calculation parameters
         a(1) = 1.0
         a(2) = 1.0
         a(3) = 0.0
         a(4) = 0.0

         b(1) = 0.0
         b(2) = 0.0
         b(3) = 1.0
         b(4) = 1.0

!     preparation of the i-th grid for space divided method
         c(1) = 0.0
         c(2) = xmax
         c(3) = 0.0
         c(4) = ymax

         cc(1) = 0.0
         cc(2) = 1.0
         cc(3) = 0.0
         cc(4) = 1.0

         kobj = nobj

!     definition of rectangular of the i-th object
         do j=1, kobj

               xr = obj(j,1)
               yr = obj(j,2)

!     check the intersection on the x-y plane
            do k = 1, 4
               d = abs(xr * a(k) + yr * b(k) - c(k))
c     d=d/sqrt(a(k)*a(k)+b(k)*b(k))
               if(d .le. obj(j,5))then
                  dd = sqrt(obj(j,5) * obj(j,5) - d * d)
                  p1 = b(k) * xr + a(k) * yr - dd
                  p2 = b(k) * xr + a(k) * yr + dd
                  min = 0.0
                  max = xmax

                  if((min .lt. p1 .and. max .gt. p1) .or.
     &                 (min .lt. p2 .and. max .gt. p2)) then

                     nobj = nobj + 1

                     obj(nobj, 1) = obj(j, 1)
     &                    + a(k) * xmax * (1. - 2. * cc(k))
                     obj(nobj, 2) = obj(j, 2)
     &                    + b(k) * ymax * (1. - 2. * cc(k))
                     obj(nobj, 3) = obj(j, 3)
                     obj(nobj, 4) = obj(j, 4)
                     obj(nobj, 5) = obj(j, 5)

                     sobj(nobj) = sobj(j)
                     iobj(nobj) = iobj(j)

c                     write(*,*) nobj
c                     write(*,*) (obj(nobj,i),i=1,5)
                  end if
               end if
            end do

            do k = 1, 2
               do l = 1, 2
                  rx = xr - c(k)
                  ry = yr - c(l + 2)
                  rr = sqrt(rx * rx + ry * ry)
                  if(rr .le. obj(j,5)) then
                     nobj = nobj + 1
                     obj(nobj, 1) = obj(j, 1) + xmax - 2. * c(k)
                     obj(nobj, 2) = obj(j, 2) + ymax - 2. * c(l + 2)
                     obj(nobj, 3) = obj(j, 3)
                     obj(nobj, 4) = obj(j, 4)
                     obj(nobj, 5) = obj(j, 5)

                     sobj(nobj) = sobj(j)
                     iobj(nobj) = iobj(j)

c                     write(*,*) nobj
c                     write(*,*) (obj(nobj,i),i=1,5)
                  end if

               end do
            end do
         end do

      end if

!     determination of the epsi (epsiron) that is used to perform the
!     Russian Roulette
!     epsi is setted to be c*(leaf_r+leaf_t)

      if(np .ne. -4 .and. np .ne. -5) then
         avelr = 0.0
         avelt = 0.0
         do i = 1, nwl
            do j = 1, nts
               avelr = avelr + lr(j,i)
               avelt = avelt + lt(j,i)
            end do
         end do
         avelr = avelr / real(nwl * nts)
         avelt = avelt / real(nwl * nts)
         epsi = 0.5 * 0.1 * (avelr + avelt)
      end if

      return
      end

      ! read leaf spectrum: reflectance & transmittance
      subroutine readspectrum(lr, lt)

      implicit none
      include "common.inc"

      ! Input
      real lr(5,*),lt(5,*)

      ! work
      integer iwl, wl, index, wlstart, wlend
      real ref, trans
      character*100 path

      wlstart = int(WL_START * 1000)
      wlend = int(WL_END * 1000)
      iwl = wlstart

      path = "Data/leaf_spectrum/spectrum."
      path = trim(path)//trim(FILE_LABEL)//".txt"

      open(24, file=path)

      do while (iwl .le. wlend)

         read(24, *) wl, ref, trans

         if (iwl .eq. wl) then
            index = iwl - wlstart + 1
            lr(1, index) = ref
            lt(1, index) = trans
         end if

         iwl = iwl + 1
      end do

      close(24)

      end

      !***********************************************************************
      ! read EF-matrix
      !
      ! *INPUT VARIABLES:
      !
      ! direction: 0 means backward, 1 means forward
      ! ps: 1 means PSI, 2 means PSII
      !
      !                             written Sicong Gao
      !                             last modified 28/02/18
      !***********************************************************************
      subroutine readEFMatrix(matrix, d, ps, folder)
        implicit none
        include "common.inc"

        integer i, d, ps,ii
        character*100 fname
        character*100 path
        character*25 folder
        real matrix(211, 351)

        path = "Data/EF_MATRIX/"//trim(folder)
        if ((d .eq. 0) .and. (ps .eq. 1)) then
           fname = trim(path)//"/MbI.csv"
        elseif ((d .eq. 0) .and. (ps .eq. 2)) then
           fname = trim(path)//"/MbII.csv"
        elseif ((d .eq. 1) .and. (ps .eq. 1)) then
           fname = trim(path)//"/MfI.csv"
        elseif ((d .eq. 1) .and. (ps .eq. 2)) then
           fname = trim(path)//"/MfII.csv"
        end if
        write(*,*) "folder:", folder
        write(*,*) "path:", path
        write(*,*) "fname:", fname
        open(9, file = fname)
        do i = 1, 211
           read(9,*) matrix(i,:)
        end do

      end

