!*******************************************************
!     This subroutine calculates photon trajectry
!     in the forest floor
!
!
!                         Written by H. Kobayashi
!                         Last modified 08/04/03
!*******************************************************

!*******************************************************
      subroutine mcflr(w, wq, x, y, z, ux, uy, uz, nscat,
     &     cmode, ulr, ult, sor, ichi, ikd,
     &     wf, nfluor, wl)
!*******************************************************

      implicit none

!     include global parameters
      include 'common.inc'
      include 'math.inc'

!     local paraeters ---------------
!
!     rnd      : random numbers
!     s        : traveling distance
!     sgm      : extinction
!     w        : weight of the photon
!     z (m)    : Height of the photon position from the uppper grass bounary
!     (-1.0<=z<0)

      integer cmode, i, nscat, mode, m
      integer ichi, ikd, ith, face
      integer ix, iy, iz

      real s, sp
      real th,ph,intv(3)
      real sgm,w,wq,conv,mgn
      real x, y, z, ux, uy, uz, uxo, uyo, uzo
      real zu, zb, uzm, cf, psh
      real ulr, ult, sor, fd
      real*8 rnd
      real*8 frnd
      real facos
      real wl

      ! SIF
      integer nfluor
      real wf, oldw, par
      real uxl, uyl, uzl
      real uxi, uyi, uzi

!     tentative fd
      fd = 0.0

      conv = 1.d-8
      mgn = 1.e-2
      uzm = 0.0174524

!     clumping factor
      cf = sbar(1)

!     mode = 1 lambertian
      mode = 1

!     uppper and bottom of the forest floor layer
      zu = 0.0
      zb = -1.0

      intv(1) = xmax
      intv(2) = ymax
      intv(3) = abs(zu - zb)

!     Monte Carlo loop
      do

         rnd = frnd()
         th = facos(uz)
         ith = int(th * 180. / pi)
         sgm = gtblf(ith) * gLAI
         sgm = max(1.e-5, sgm)
         s = -log(max(rnd,1.e-10)) / sgm

!     check the intersection
         call planes(sp, x, y, z, ux, uy, uz, 0.0, 0.0,-1.0, face, intv)

!     to avoid underflow,
         sp = max(sp, 1.e-10)

         if(s .lt. sp) then
		 ! over-story
		 ! scattering

            x = x + s * ux
            y = y + s * uy
            z = z + s * uz

!     recollision
            do

!     fpar sampling
               par =  w * wq * (1. - ulr - ult)
               ffpr = ffpr + par

               ix = int(x * res) + 1
               iy = int(y * res) + 1
               iz = int(z) + 1

               ix = min(ix, size)
               iy = min(iy, size)

               apf(ix, iy) = apf(ix, iy) + w * wq * (1. - ult - ulr)
               apfd(ix, iy) = apfd(ix, iy) + w * wq * (1. -ult - ulr)
     &              * min(real(nscat), 1.)

               call calParLeaf(ix, iy, iz, wl, (w + wf) * wq)
               oldw = w

               w = w * (ult + ulr)
               !oldw = w

               if (w .gt. 0) then
                  nscat = nscat + 1
               end if

               call vegrroulette(w, epsi)

               if (w .lt. conv) w = 0
               if ((w .lt. conv) .and. (wf .lt. conv)) return

               ! for broad leaf, psh = 0, exit for sure
               psh = (1. - 4. * cf)
               if(frnd() .ge. psh) exit

            end do
            ! save incident vectors
            uxi = ux
            uyi = uy
            uzi = uz

            call vegrad(w, x, y, z, ux, uy, uz,
     &              ulr, ult, 4, 1.0, fd, ichi, ikd)

!     new direction
            call scatvec(ulr, ult, uxl, uyl, uzl, ux, uy, uz,
     &                   uxo, uyo, uzo, mf)
            

            ux = uxo
            uy = uyo
            uz = uzo
            if(abs(uz) .lt. uzm) uz = sign(uzm, uz)

            ! ************************************************************
            ! SIF
            par = (oldw + wf) * wq * (1. - ulr - ult)
            wf = wf * (ult + ulr)
            W_G = W_G + oldw
            WF_G = WF_G + wf
            PPFD_G = PPFD_G + par
            !par = 0.0
            call fluorrad(wf, oldw, wq, par, wl,
     &                    x, y, z,
     &                    uxi, uyi, uzi, uxl, uyl, uzl,
     &                    uxo, uyo, uzo,
     &                    ulr, ult, (ulr + ult),
     &                    nfluor, 0, 2,
     &                    4, ikd, ichi, fd)
            call vegrroulette(wf, epsi)
            if (wf .lt. conv) wf = 0.0
            if ((w .lt. conv) .and. (wf .lt. conv)) return
            ! write(*,*) "** wf = ", wf
            ! ************************************************************

         else
		 ! soil surface
            x = x + (sp + mgn) * ux
            y = y + (sp + mgn) * uy
            z = z + sp * uz
            x = x - (aint(x / xmax) - 0.5 + sign(0.5, x)) * xmax
            y = y - (aint(y / ymax) - 0.5 + sign(0.5, y)) * ymax

            if(z .le. zb)then

!     surface boudary refletance mode 1: Lambertian, 2: RPV model, 3: DSM model

!     fpar sampling
               sfpr = sfpr + w * wq * (1. - sor)
               ix = int(x * res) + 1
               iy = int(y * res) + 1

               ix=min(ix, 300)
               iy=min(iy, 300)
               aps(ix, iy) = aps(ix, iy) + w * wq * (1. - ult - ulr)
!     surface downward flux
               sfdir(ix, iy) = sfdir(ix, iy)
     &              + w * wq * (1. - min(nscat, 1))
               sfdif(ix, iy) = sfdif(ix, iy)
     &              + w * wq * min(nscat, 1)
               oldw = w
               w = w * sor
               !oldw = w

               uxi = ux
               uyi = uy
               uzi = uz
               call vegrad(w, x, y, z, ux, uy, uz,
     &              ulr, ult, 5, 1.0, fd, ichi, ikd)

               call srfref(ux, uy, uz, sor, nscat, mode)
               z = zb + mgn
               uxo = ux
               uyo = uy
               uzo = uz
               uxl = 0.0
               uyl = 0.0
               uzl = -1.0

               ! ************************************************************
               ! SIF
               par = 0.0
               wf = wf * sor
               PPFD_F = PPFD_F + wf + oldw
               W_F = W_F + oldw
               WF_F = WF_F + wf
               call fluorrad(wf, oldw, wq, par, wl,
     &                       x, y, z,
     &                       uxi, uyi, uzi, uxl, uyl, uzl,
     &                       uxo, uyo, uzo,
     &                       sor, 1 - sor, 1.0,
     &                       nfluor, 1, 2,
     &                       5, ikd, ichi, fd)
               call vegrroulette(wf, epsi)
               if (wf .lt. conv) wf = 0.0
               if ((w .lt. conv) .and. (wf .lt. conv)) return
               ! ************************************************************

            else if(z. ge. zu) then
               return
            end if

         end if
      end do

 300  return
      end
