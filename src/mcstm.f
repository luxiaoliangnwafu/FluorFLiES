!*******************************************************
!     This subroutine calculates photon reflection
!     on the stem surface
!
!     Written by H. Kobayashi
!     Last modified 08/04/01
!*******************************************************
!*******************************************************
      subroutine mcstm (w, wq, x, y, z, ux, uy, uz, nscat,
     &     cmode, tobj, face, str, ichi, ikd,
     &     wf, nfluor, wl)
!*******************************************************

      implicit none

!     include global parameters
      include 'common.inc'
      include 'math.inc'

      integer face, cmode, iz
      integer ichi, ikd, nscat, cb

      real x, y, z, ux, uy, uz
      real rx, ry, rr, a, fd
      real nux, nuy, uzm
      real w, wq, tobj(5)
      real*8 frnd
      real*8 rnd
      real str, th, ph
      real fsin, fcos, facos
      real conv, mgn

      ! SIF
      integer nfluor
      real wf, oldw, par, wl
      ! leaf normal vector
      real uxl, uyl, uzl, uxi, uyi, uzi, uxo, uyo, uzo
      real flagfscat

      conv = 1.d-8
      fd = 0.0
      mgn = 1.d-2
      uzm = 0.0174524
      cb = 3

      uxi = ux
      uyi = uy
      uzi = uz
      uxl = 1.0
      uyl = 1.0
      uzl = 1.0

!     reflectance at the side of stem
      if(face .eq. 1) then
!     stem normal vector

         rx = (x - tobj(1))
         ry = (y - tobj(2))
         rr = sqrt(rx *rx + ry *ry)
         nux = rx / rr
         nuy = ry / rr

         th = 0.5 * facos(1. - 2. * real(frnd()))
         ph = 2.0 * pi * real(frnd())

         a = facos(nux)
         cb = 6

         call trans(nux, nuy, 0.0, th, ph, ux, uy, uz)
         if(abs(uz) .lt. uzm) uz = sign(uzm, uz)

!     reflectance at the bottom of stem
      else if(face .eq. 2) then

         th = 0.5 * pi + 0.5 * facos(1. - 2. * real(frnd()))
         ph = 2.0 * pi * real(frnd())
         ux = fsin(th) * fcos(ph)
         uy = fsin(th) * fsin(ph)
         uz = fcos(th)
         z = z - mgn
         if(abs(uz) .lt. uzm) uz = sign(uzm, uz)
         a = 1.

!     reflectance at the top of stem
      else

         th = 0.5 * facos(1. - 2. * real(frnd()))
         ph = 2.0 * pi * real(frnd())
         ux = fsin(th) * fcos(ph)
         uy = fsin(th) * fsin(ph)
         uz = fcos(th)
         z = z + mgn
         if(abs(uz) .lt. uzm) uz = sign(uzm, uz)
         a = 1.
      end if
      
      ! record outgoing vector
      uxo = ux
      uyo = uy
      uzo = uz

!     fpar samping (leave or branch)
      ! bfpr: none leaf absorbed radiation
      ! apnp: vertical radiation at k meter
      bfpr = bfpr + w * wq * (1. - str)
      iz = int(z) + 1
      apnp(iz) = apnp(iz) + w * wq * (1. - str)

      oldw = w
      w = w * str
      !oldw = w

      if (w .gt. 0) then
         nscat = nscat + 1
      end if

!     Russian Roulette
      call vegrroulette(w, epsi)

      if (w .lt. conv) w = 0.0
      ! ************************************************************
      ! SIF
      par = 0.0
      wf = wf * str
      PPFD_S = PPFD_S + (wf + w) * wq
      W_S = W_S + oldw
      WF_S = WF_S + wf
      
      call fluorrad(wf, oldw, wq, par, wl,
     &              x, y, z,
     &              uxi, uyi, uzi, uxl, uyl, uzl,
     &              uxo, uyo, uzo,
     &              str, 1 - str, 1.0,
     &              nfluor, 1, 1,
     &              cb, ikd, ichi, fd)
       call vegrroulette(wf, epsi)
       if (wf .lt. conv) wf = 0.0
      ! ************************************************************

      call vegrad(w, x, y, z, ux, uy, uz,
     &     1.0, 0.0, cb, a, fd, ichi, ikd)

      if (w .lt. conv) w = 0
      if ((w .lt. conv) .and. (wf .lt. conv)) return

      return
      end
