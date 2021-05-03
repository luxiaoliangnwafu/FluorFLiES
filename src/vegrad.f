!**********************************************
!     3-Dimensional monte carlo
!     calculation of the radiation field
!     by H. Kobayashi
!     modified 04/12/08
!     last modified 08/04/07
!**********************************************
      subroutine vegrad(w, x, y, z, ux, uy, uz,
     &     lr, lt, cb, a, fd, ichi, ikd)

      implicit none

      include 'common.inc'
      include 'math.inc'

!     cb:cb = 1 (overstory),cb = 2 (branch)
!     cb = 3 (stem top/bottom),cb = 4 (understory), cb = 5 (soil surface)
!     cb = 6 (stem side)
!     a is only used for lambertian reflection from the stem side
!     other case "a" should be 1.0

      integer i, j, cb
      integer ichi, ikd, ith, ithr
      integer sflag
      integer ix, iy

      real w, x, y, z, ux, uy, uz
      real xt, yt, zt
      real xflag, yflag
      real xr, yr, cosa
      real lr, lt, fd, pf
      real tho, thi, th, thr, phr, rr, ph
      real thrad,phrad
      real ff(6), ua, a, Id
      real taua, tauc, lrds
      real ch, hk, ga, af, nf
      real fexp, fsin, fcos, facos

      real conv

      ! for SIF, if wf excist, w extinct, just pass
      conv = 1.d-8
      if (w .lt. conv) return

      Id = 0.0

      tauc = 0.0
      taua = 0.0

      ff(1) = 1.0
      ff(2) = 1.0
      ff(3) = 0.0
      ff(4) = 1.0
      ff(5) = 0.0
      ff(6) = 0.0

!     if x, and y is outside area
      xt = x - (aint(x / xmax) - 0.5 + sign(0.5, x)) * xmax
      yt = y - (aint(y / ymax) - 0.5 + sign(0.5, y)) * ymax
      zt = z

!     lrds: leaf radius (0.1m)
      lrds = 0.1

      ! nangc: view angles total numbers
      do i = 1, nangc

!     preparation of Haple-type hotspot function
         cosa = ux * uxrc(i)
         cosa = cosa + uy * uyrc(i)
         cosa = cosa + uz * uzrc(i)
         cosa = sign(min(cosa, 0.999999),cosa)
         af = facos(cosa)

         th = facos(uz)
         ith = int(th * 180. / pi)
         thr = facos(uzrc(i))
         ithr = int(thr * 180. / pi)

!     Hapke, types hot spot function
!     af is converted to the opposite angle of scattering angle
!     See Kobayashi and Iwabuchi (2008) Remote Sensing of Environment
         ! scattering occurs on the same interface(canopy, stem, floor)
         af = pi - af
         ga = dlt(cb, 1) * (gtblc(ith) + gtblc(ithr)) * 0.5
         ga = ga + dlt(cb, 2) * (gtblb(ith) + gtblb(ithr)) * 0.5
         ga = ga + dlt(cb, 3) * (gtblc(ith) + gtblc(ithr)) * 0.25
         ga = ga + dlt(cb, 3) * (gtblb(ith) + gtblb(ithr)) * 0.25
         ga = ga + dlt(cb, 4) * (gtblf(ith) + gtblf(ithr)) * 0.5
         ga = ga + dlt(cb, 5) * (gtblf(ith) + gtblf(ithr)) * 0.5
!     for stem side, the hotspot effect will be ignored so hk=1.0
         ga = ga + dlt(cb, 6) * 1.e-6

         ua = 0.0
         do j = 1, nts
            ua = ua + u(j)
         end do
         ua = ua / real(nts)

        ! calculate u (leaf area density)
         ch = gLAI * (dlt(cb, 4) + dlt(cb, 5))
         ch = ch + ua * (dlt(cb, 1) + dlt(cb, 2) + dlt(cb, 3))
         ch = ch + dlt(cb, 6) * 1.e-6

         ch = ch * ga * lrds * 0.5
         ch = 1. / ch
         hk = 1. / (1. + ch * fsin(af * 0.5)/fcos(af * 0.5))
         hk = 1. - hk

         call vegtrace(tauc,xt,yt,zt,uxrc(i),uyrc(i),uzrc(i),sflag)
         call mc1descape(xt, yt, 0.0, uxrc(i), uyrc(i), uzrc(i), 1,
     &        ichi, ikd, xr, yr, taua)

         thi = facos(uz)
         tho = facos(uzrc(i))
         nf = fsin(thi) * fsin(tho)
         nf = max(nf, 1.d-8)
         phr = ux * uxrc(i) + uy * uyrc(i)
         phr = phr / nf
         phr = facos(phr)

         call fpf(pf, thi, tho, phr, lr, lt, cb)

         tauc = tauc + (- zt / uzrc(i))
     &        * gtblf(ithr) *gLAI * (dlt(cb, 4) + dlt(cb, 5))
         tauc = tauc * hk

         if(cb .eq. 6) then
            rr = sqrt(uxrc(i) * uxrc(i) + uyrc(i) * uyrc(i))
            rr = max(1.d-3, rr)
            ph = uxrc(i) / rr
            ph = facos(ph)
            a = abs(fsin(facos(uzrc(i))) * fcos(a - ph))
         end if
         !write(*,*) "veg:tauc=",tauc

         Id = 0.0
         Id = ff(cb) * w * pf * fexp(-tauc) / fcos(thr)
         Id = Id + (1.- ff(cb)) * w * fexp(-tauc) / pi
         Id = Id * abs(a)
         ! floor, strm, branch: fd = 0
         Id = Id * (1. - fd)
         ! sflag = 0 means stem collision
         Id = Id * real(sflag)

         brf(1, i) = brf(1, i) + Id
         brfc(1, i) = brfc(1, i) + Id * dlt(cb, 1)
         brfs(1, i) = brfs(1, i) + Id * (dlt(cb, 2) + dlt(cb, 3))
         brff(1, i) = brff(1, i) + Id * (dlt(cb, 4) + dlt(cb, 5))

!      Nadir image
         ix = int(xt * res) + 1
         iy = int(yt * res) + 1

         ix = min(ix, size)
         iy = min(iy, size)

         refl(1, ix, iy) = refl(1, ix, iy) + Id
         irefl(1, ix, iy) = irefl(1, ix, iy) + 1

         Id = Id * fexp(-taua)
         brf(2, i) = brf(2, i) + Id
         brfc(2, i) = brfc(2, i) + Id * dlt(cb, 1)
         brfs(2, i) = brfs(2, i) + Id * (dlt(cb, 2) + dlt(cb, 3))
         brff(2, i) = brff(2, i) + Id * (dlt(cb, 4) + dlt(cb, 5))

         refl(2, ix, iy) = refl(2, ix, iy) + Id
         irefl(2, ix, iy) = irefl(2, ix, iy) + 1

      end do
      return
      end
