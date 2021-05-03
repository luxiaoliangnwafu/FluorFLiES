!***************************************************
      subroutine fpf(gm, thi, tho, phr, lr, lt, cb)

!     calculate the phase function from LUT
!***************************************************
      implicit none

      include 'common.inc'
      include 'math.inc'
      
!     cb:cb = 1 (overstory),cb = 2 (branch)
!     cb = 3 (stem top/bottom),cb = 4 (understory), cb = 5 (soil surface)
!     cb = 6 (stem side)

      integer i, j, k, l, m, n, cb
      real gm, th1, th2, ph, thi, tho, phr
      real gmr, gmt, lr, lt, gfunc
      real tr(2), tt(2), ttt(2), trr(2)
      real gmrx(2, 2, 2), gmtx(2, 2, 2)

!     change angles rad to degree
      th1 = min(thi * 180. / pi, 179.99999)
      th2 = min(tho * 180. / pi, 179.99999)
      ph = phr * 180. / pi

!     all angles are divided by 10.
      th1 = th1 / 10.
      th2 = th2 / 10.
      ph = ph / 10.

      i = int(aint(th1)) + 1
      j = int(aint(th2)) + 1
      k = int(aint(ph)) + 1
      
c      write(*,*) th1,th2,ph,i,j,k
c      pause

      do l = 0, 1
         do m = 0, 1
            do n = 0, 1
               gmrx(l + 1, m + 1 ,n + 1) =
     &              dlt(cb, 1) * gmrc(i + l ,j + m, k + n) 
     &              + dlt(cb, 2) * gmrb(i + l ,j + m, k + n)
     &              + dlt(cb, 4) * gmrf(i + l ,j + m, k + n)
               gmtx(l + 1, m + 1, n + 1) =
     &              dlt(cb, 1) * gmtc(i + l ,j + m, k + n) 
     &              + dlt(cb, 2) * gmtb(i + l ,j + m, k + n) 
     &              + dlt(cb, 4) * gmtf(i + l ,j + m, k + n)
               
c               write(*,*) i,j,k, gmrx(l + 1, m + 1 ,n + 1), 
c     &              gmtx(l + 1, m + 1, n + 1),cb

            end do
         end do
      end do

!     bi-linear over th1 dimension 
      do n = 1, 2 
         do l = 1, 2
            tr(l) = gmrx(1, l, n) * (real(i) - th1) 
     &           + gmrx(2, l, n) * (th1 - real(i - 1))
            tt(l) = gmtx(1, l, n) * (real(i) - th1) 
     &           + gmtx(2, l, n) * (th1 - real(i - 1))      
c            write(*,*) "flag1"
         end do
         
!     bi-linear over th1, th2 dimension
         
         trr(n) = tr(1) * (real(j) - th2) 
     &        + tr(2) * (th2 - real(j - 1))
         ttt(n) = tt(1) * (real(j) - th2) 
     &        + tt(2) * (th2 - real(j - 1))
         
c         write(*,*) "flag2",trr(1),ttt(1)
      end do

!     bi-linear over (th1+th2) - ph plan
      gmr = trr(1) * (real(k) - ph) + trr(2) * (ph - real(k - 1))
      gmt = ttt(1) * (real(k) - ph) + ttt(2) * (ph - real(k - 1))      

c      write(*,*) gmr, gmt,trr(1),trr(2),tr(1),tr(2)

      gm = lr * gmr + lt * gmt 
      gfunc = dlt(cb, 1) * gtblc(int(th1 * 10.)) +
     &     dlt(cb, 2) * gtblb(int(th1 * 10.)) +
     &     dlt(cb, 4) * gtblf(int(th1 * 10.)) +
     &     (dlt(cb, 3) + dlt(cb, 5) + dlt(cb, 6)) * 1.0
      gm = (1. / gfunc) * gm / ((lr + lt) * pi)

c      gm = gm * 2. / ((lr + lt) * pi)

      return
      end
