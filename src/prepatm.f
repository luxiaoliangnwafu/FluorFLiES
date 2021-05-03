! **********************************************
!     preparation of the atmospheric parameters
!
!     By Hideki Kobayashi
! ***********************************************

      subroutine prepatm(imode, wq, RQ, span, spcq, np, npproc,
     &     wl0, wls, iwl, npl, nkd, wkd, ext, rfname, fname,
     &     nmix, cflg, re, G, Qabs, Qext, Qext_ref, d,
     &     cbnz, ctnz, taur, ctaur)

      implicit none

      include 'common.inc'
      include 'math.inc'

      integer nkd, iwl
      integer cflg, nmix, imode
      integer np, npl(1000),npproc
      integer knmix, knzext, knang
      parameter(knmix = 10, knzext = 200, knang = 2000)
      real wl0, wls, d, cbnz, ctnz, rat(2)
      real wkd(knkd), anG(knang)
      real re(knmix), Qext_ref(knmix),Qext(knmix),G(knmix),Qabs(knmix)
      real omg(knmix), phs(knmix, knang)
      real ext(knmix, knzext)
      real RQ, RF, span(1000)
      real spcq(1000), spcf(1000), wq
      real ctaur, taur

      character*81 rfname, fname(knmix)

c      write(*,*) "prep 1"

      ! imode = par or sw , RQ is photon flux at the top of atmosphere
      if(imode .ne. 1) then
!     weight of photon for PPFD
        ! ppfd is photosynthetic photon flux density
         !write(*,*) "iwl = ", iwl
         !write(*,*) "spcq = ", spcq(iwl)
         !wq = span(iwl) * spcq(iwl) / RQ
         wq = span(iwl) * spcq(iwl) / 500.0
         wq = wq * real(npproc) / real(npl(iwl)) / 2
         !write(*,*) "npproc = ", npproc
         !write(*,*) "npl(iwl) = ", npl(iwl)
         wl0 = wls + iwl * span(iwl)
         !if(iwl .le. 20) then
         !   wl0 = wls + real(iwl - 1) * span(iwl) + span(iwl) / 2.
         !else
         !   wl0 = 0.7 + real(iwl - 21) * span(iwl) + span(iwl) / 2.
         !end if
      end if

c      write(*,*) "prep 2"

!     get atmospheric parameters
      call getatm(wl0, nkd, wkd, ext, rfname)

c      write(*,*) "prep 3"

      call getphs(wl0, span(1), nmix + cflg, re, Qext_ref, Qext, omg,
     &     G, Qabs, ang, phs, fname, imode)

c      write(*,*) "prep 4"

!     get the ratio of extinction corf. of each aerosol
      rat(1) = 1.
      rat(2) = 2.

c      write(*,*) "prep 5"
      call getext(nmix, cflg, ext, Qext, Qext_ref, d, cbnz, ctnz,
     &     taur, ctaur, rat)

c      write(*,*) "prep 6"

      call mc1doptics(ext, omg, phs, ang, knmix, nmix + cflg, iwl)

c      write(*,*) "prep 7"

      return
      end
