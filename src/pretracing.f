!*************************************************************************
! pre-tracing method. Calculate incident par and fluorescence (400-750 nm)
! and then reset
!
! *INPUT VARIABLES:
!
!
!                             written Sicong Gao
!                             last modified 26/02/18
!**************************************************************************
      subroutine pretracing(nwl, npl, wl0, wls, span, up,
     &                      sinf0, cosf0, sinq0, cosq0,
     &                      wkd, conv,
     &                      lr, lt, str, ulr, ult, sor,
     &                      cmode, imode,
     &                      RF, RQ, spcq,
     &                      np, npproc, nkd, nmix, cflg,
     &                      ext,
     &                      rfname, fname,
     &                      re, G, Qabs, Qext, Qext_ref,
     &                      d, cbnz, ctnz, taur, ctaur,
     &                      scmpf, scmpp)

      implicit none

      include 'common.inc'
      include 'math.inc'

      ! input
      integer nwl, up, cmode, imode
      integer npl(1000)
      integer np, npproc
      integer nkd, nmix, cflg
      integer cbnz,ctnz, taur, ctaur
      real wl0, wls, span(1000), spcq(1000)
      real conv, RF, RQ, d
      real sinf0, cosf0, sinq0, cosq0
      real wkd(knkd)
      real lr(5,1000), lt(5,1000), ulr(1000), ult(1000)
      real str(5,1000), sor(1000)
      integer knmix,knang,knzext
      parameter (knmix=10, knang=2000, knzext=200)
      real ext(knmix, knzext)
      real re(knmix),Qext_ref(knmix),Qext(knmix),G(knmix),Qabs(knmix)
      character*81 fname(knmix),rfname
      real*8 scmpf(3,1000), scmpp(3,1000)

      ! work
      integer ip, iz, i, j, k, iwl
      real ftau, chi, wl, pixnp
      real*8 frnd, rnd
      !     --- photon data  ---
      integer nscat, nscata, ichi, ikd
      real x, y, z, ux, uy, uz, w, wq
      !     --- SIF data  ---
      real wf
      integer nfluor
      !     --- spectral data  ---
      real tlr(5), tlt(5), tstr(5)
      !     --- collect data  ---
      real*8 rflx, rbflx, rdflx
      real*8 tflx, bflx, dflx, tpfd, bpfd, dpfd
      !     --- random sequence  ---
      integer nrand, lrand, cplrand, cnt, cpcnt
      parameter (nrand=1000)
      dimension lrand(-249:nrand)
      dimension cplrand(-249:nrand)
      common /rnd/lrand,cnt

      ! *init
      rflx =0.0                 ! total reflected irradiaice
      rbflx = 0.0               ! beam reflected irradiance
      rdflx = 0.0               ! diffuse reflected irradiance

      ! Irradiance
      tflx = 0.0                ! Total downward flux at top of canopy
      bflx = 0.0                ! Beam downward flux at TOC
      dflx = 0.0                ! Diffuse downward flux at TOC

      ! Photon flux density
      tpfd = 0.0                ! Total downward PFD at top of canopy
      bpfd = 0.0                ! Beam downward PFD at top of canopy
      dpfd = 0.0                ! Diffuse downward PDF at top of

      ! record list of random and pointer place
      cplrand = lrand
      cpcnt = cnt

      write(*,*) "Start pre-tracing..."

      write(*,*) "# of sampling wavelength"
      do iwl = 1, nwl

         write(*, *) iwl, "  of", nwl, " (", npl(iwl),") "

         ! calculate the wl
         !if (iwl .le. 20) then
         !   wl = wls + span(iwl) * real(iwl) - span(iwl) / 2.0
         !else
         !   wl = 0.7 + span(iwl) * real(iwl - 20) - span(iwl) / 2.0
         !end if
         wl = WL_START + span(iwl) * real(iwl)
         ! record current iwl, use later
         ! CURRENT_WL = wl

         write(*,*) "WL (nm):", int(wl * 1000)

         ! pre-tracing only in 400-700 nm
         if ((wl .lt. 0.4) .or. (wl .gt. 0.7)) then
            write(*,*)  "pre-tracing won't run in current wavelength."
            CYCLE
         end if

         ! preparation of the atmospheric optical parameters
         call  prepatm(imode, wq, RQ, span, spcq, np, npproc,
     &                 wl0, wls, iwl, npl, nkd, wkd, ext, rfname, fname,
     &                 nmix, cflg, re, G, Qabs, Qext, Qext_ref, d,
     &                 cbnz, ctnz, taur, ctaur)

         if(npl(iwl) .eq. 0) then
            write(*, *) "npl error !"
            exit
         end if

         ! loop for single wavelength
         do ip = 1, npl(iwl)

            if(mod(ip, up) .eq. 0)  write(*,*) ip," of ", npl(iwl)

            w = 1.0
            x = xmax * real(frnd())
            y = ymax * real(frnd())
            z = zgrd(nz)
            ux = sinq0 * cosf0
            uy = sinq0 * sinf0
            uz = cosq0

            ! SIF
            wf = 0
            ! fluorescence times
            nfluor = 0

            ftau = -log(max(1.0e-35, real(frnd())))
            chi = 1.0
            ichi = 1
            iz = nz
            nscat = 0
            ! determination of the k-term
            rnd = real(frnd())
            ikd = 1
            do while(rnd .gt. wkd(ikd))
               ikd = ikd + 1
            end do

            do               ! Monte Carlo core loop
              call mc1dtrace(w, wq, x, y, z, ux, uy, uz, ftau,
     &                       chi, iz, ikd, nscat, ichi, iwl)
              ! write(*,*) "nscat = ",nscat
              ! write(*,*) "iz = ",iz

              ! nscata is the scatter times in ATM, nscat is the total scatter times
              nscata = nscat
              if (w .le. 0.0 .or. iz .gt. nz) exit

              ! scmpf is the energy flux
              ! scmpf(1, wl) means the total wl energy flux
              ! scmpf(2, wl) means the beam wl energy flux
              ! scmpf(3, wl) means the diffuse wl energy flux
              ! nscat > 0 means the photo has been collision
              ! write(*,*) "w, wq",w,wq
              scmpf(1, iwl) = scmpf(1, iwl) + w
              scmpf(2, iwl) = scmpf(2, iwl)
     &                      + w * (1. - min(real(nscat), 1.))
              scmpf(3, iwl) = scmpf(3, iwl)
     &                      + w * min(real(nscat), 1. )

              ! scmpp is the photo flux
              ! scmpp(1, wl) means the total wl photon flux
              ! scmpp(2, wl) means the beam wl photon flux
              ! scmpp(3, wl) means the diffuse wl photon flux
              scmpp(1, iwl) = scmpp(1, iwl) + w * wq
              scmpp(2, iwl) = scmpp(2, iwl)
     &                 + w * wq * (1. - min(real(nscat), 1.))
              scmpp(3, iwl) = scmpp(3, iwl)
     &                 + w * wq * min(real(nscat), 1.)

              !     surface reflection
              !     3-D surface
              do i = 1, nts
                tlr(i) = lr(i,iwl)
                tlt(i) = lt(i,iwl)
                tstr(i) = str(i,iwl)
              end do
                
              init_nscat = nscat
              ! call the canopy radiation transfer module
              ! cmode : BRF NADIR IMAGE
              ! tlr, tlt, ulr, ult : leaf / understory leaf reflectance & transmittance
              ! tstr: strm reflectance
              ! sor: soil reflectance
              ! ichi, ikd: ATM sth.
              call canort(x, y, zmax - 1.d-2, ux, uy, uz,
     &                    w, wq, cmode, nscat,
     &                    tlr, tlt, ulr(iwl), ult(iwl), tstr,
     &                    sor(iwl), ichi, ikd,
     &                    wf, nfluor, wl)
              if(w .lt. conv) exit
              iz = 1
              rflx = rflx + w
              rbflx = rbflx + w * (1. - min(real(nscata), 1.))
              rdflx = rdflx + w * min(real(nscata), 1. )
           end do           ! single photon loop
        end do              ! photon loop

      !     summary of spectral flx
        tflx = tflx + scmpf(1,iwl)
        bflx = bflx + scmpf(2,iwl)
        dflx = dflx + scmpf(3,iwl)

        tpfd = tpfd + scmpp(1,iwl)
        bpfd = bpfd + scmpp(2,iwl)
        dpfd = dpfd + scmpp(3,iwl)

      end do                 ! loop for nwl (spectral loop)

      write(*,*) "Finish pre-tracing..."

      ! reset variables
      ! *init
      rflx =0.0                 ! total reflected irradiaice
      rbflx = 0.0               ! beam reflected irradiance
      rdflx = 0.0               ! diffuse reflected irradiance

      ! Irradiance
      tflx = 0.0                ! Total downward flux at top of canopy
      bflx = 0.0                ! Beam downward flux at TOC
      dflx = 0.0                ! Diffuse downward flux at TOC

      ! Photon flux density
      tpfd = 0.0                ! Total downward PFD at top of canopy
      bpfd = 0.0                ! Beam downward PFD at top of canopy
      dpfd = 0.0                ! Diffuse downward PDF at top of

      do i = 1, 3
         do j = 1, 1000
            scmpf(i,j) = 0.0
            scmpp(i,j) = 0.0
         end do
      end do

      lrand = cplrand
      cnt = cpcnt
      !write(*,*) "*cnt = ",cnt

      pixnp = real(np) / real(size * size)


      do i = 1, size
         do j = 1, size
            do k = 1, 100
              !write(*,*) PAR_LEAF(i,j,k)
               PAR_LEAF(i,j,k) = PAR_LEAF(i,j,k) * 1000 /
     &                  (HC * NAVO * pi * LEAFAREA * LEAFAREA) *
     &                  RQ * abs(cosq0)
               if (PAR_LEAF(i,j,k) .ne. 0) then
                  !write(*,*) "No 0:",(PAR_LEAF(i,j,k))
               end if
            end do
         end do
      end do

      !write(*,*) "F_scat"
      !do i = 1, 211
          !write(*,*) F_scat(i)
      !end do

      !write(*,*) "F_excited"
      !do i = 1, 211
          !write(*,*) F_excited(i)
      !end do

      !write(*,*) "QFS"
      !do i = 1, 211
          !write(*,*) QFS(i,1)
      !end do

      end


      subroutine writerand(arr,id)

      implicit none

      include 'common.inc'
      integer arr, nrand
      character*1 id
      parameter (nrand=1000)
      dimension arr(-249:nrand)
      integer i

      open(12,file = "rand"//id)

      do i=-249, 1000
         write(12,*) arr(i)
      end do
      close(12)

      end
