!***************************************************
      subroutine ipf
!     preparation of LUT for kernel of phase function
!     by Hideki Kobayashi
!     08/04/10
!****************************************************
      implicit none
     
      include 'common.inc'
      include 'math.inc'

      integer i, j, k
      integer op
      real th1, th2, ph
      real gmr, gmt

      character*81 rfile(3), tfile(3)

      rfile(1) = 'Data/gmr_uni.txt'
      rfile(2) = 'Data/gmr_plano.txt'
      rfile(3) = 'Data/gmr_erect.txt'

      tfile(1) = 'Data/gmt_uni.txt'
      tfile(2) = 'Data/gmt_plano.txt'
      tfile(3) = 'Data/gmt_erect.txt'


!      leaf angle distribution
c      mc = 1
c      mb = 1
c      mf = 1
!     op=1 kernel calculation,else read from txt file
      op = 2

      if(op .eq. 1) then
         
         do i = 1, 19
            th1 = real(i - 1) * 10.* pi / 180.
            
            do j = 1, 19
               th2 = real(j - 1) * 10.* pi / 180.
               
               do k = 1, 37   
                  ph = real(k - 1) * 10.* pi / 180.  
                  
!     canopy
                  call  gmkernel(gmr,gmt,th1, th2, ph, mc)
                  gmrc(i,j,k) = gmr
                  gmtc(i,j,k) = gmt
                  
!     branch
                  call  gmkernel(gmr,gmt,th1, th2, ph, mb)
                  gmrb(i,j,k) = gmr
                  gmtb(i,j,k) = gmt
                  
!     forest floor
                  call  gmkernel(gmr,gmt,th1, th2, ph, mf)
                  gmrf(i,j,k) = gmr
                  gmtf(i,j,k) = gmt

               end do
               write(*,*) i,j, " of 19, 19" 
            end do
         end do
         

c         write(*,*) "after",gmtc(1,1,1)
c         pause

! for reflection 
         open(12, file="gmrc.txt")
         do i = 1, 19
            do j = 1, 19
               write(12, 100) (gmrc(i,j,k), k = 1, 37)
            end do
         end do
         close(12)

         open(13, file="gmrb.txt")
         do i = 1, 19
            do j = 1, 19
               write(13, 100) (gmrb(i,j,k), k = 1, 37)
            end do
         end do
         close(13)

         open(14, file="gmrf.txt")
         do i = 1, 19
            do j = 1, 19
               write(14, 100) (gmrf(i,j,k), k = 1, 37)
            end do
         end do
         close(14)


! for transmission
         open(12, file="gmtc.txt")
         do i = 1, 19
            do j = 1, 19
               write(12, 100) (gmtc(i,j,k), k = 1, 37)
            end do
         end do
         close(12)

         open(13, file="gmtb.txt")
         do i = 1, 19
            do j = 1, 19
               write(13, 100) (gmtb(i,j,k), k = 1, 37)
            end do
         end do
         close(13)

         open(14, file="gmtf.txt")
         do i = 1, 19
            do j = 1, 19
               write(14, 100) (gmtf(i,j,k), k = 1, 37)
            end do
         end do
         close(14)

 100     format(37(F12.8))
         
!      read from txt file
         else

! reflection
            open(12, file=rfile(mc))
            do i = 1, 19
               do j = 1, 19
                  read(12, 100) (gmrc(i,j,k), k = 1, 37)
               end do
            end do
            close(12)
      
            open(13, file=rfile(mb))
            do i = 1, 19
               do j = 1, 19
                  read(13, 100) (gmrb(i,j,k), k = 1, 37)
               end do
            end do
            close(13)

            open(14, file=rfile(mf)) 
            do i = 1, 19
               do j = 1, 19
                  read(14, 100) (gmrf(i,j,k), k = 1, 37)
               end do
            end do
            close(14)
      
! transmission
            open(12, file=tfile(mc))
            do i = 1, 19
               do j = 1, 19
                  read(12, 100) (gmtc(i,j,k), k = 1, 37)
               end do
            end do
            close(12)

            open(13, file=tfile(mb))
            do i = 1, 19
               do j = 1, 19
                  read(13, 100) (gmtb(i,j,k), k = 1, 37)
               end do
            end do
            close(13)

            open(14, file=tfile(mf)) 
            do i = 1, 19
               do j = 1, 19
                  read(14, 100) (gmtf(i,j,k), k = 1, 37)
               end do
            end do
            close(14)

         end if
      return
      end



!****************************************************
      subroutine gmkernel(gmr,gmt,th1, th2, ph, m)
!     Mathmatics is summarized in
!     Shultis and Myneni (1988)
!     J. Q. Spec. Radi. Trans., 39(2), p. 115-129
!******************************************************
      implicit none
     
      include 'common.inc'
      include 'math.inc'

      integer i, j, k, m

      real a, th1, th2, ph, thl1, thl2, phl
      real af1, af2, cosf, gl1, gl2
      real gmr, gmt, gmrp(361), gmtp(361), gmrt(181), gmtt(181)
      real faf, fgl

!     get phase angle 
!     a = faf(th1, th2, ph)

      gmr = 0.0
      gmt = 0.0

      do i = 1, 361
         gmrp(i) = 0.0
         gmtp(i) = 0.0
      end do

      do i = 1, 361

         phl = (real(i) - 1.0) * pi / 180.

         do j = 1, 91
            
            thl1 =(real(j) - 1.0) * pi / 180.
            af1 = faf(th1, thl1, phl)
            af2 = faf(th2, thl1, abs(phl - ph))
            cosf = cos(af1) * cos(af2)

            if(cosf .lt. 0.0)then
               gmrt(j) = abs(cosf)
               gmtt(j) = 0.0
            else
               gmrt(j) = 0.0
               gmtt(j) = cosf
            end if
!            write(*,*) gmrt(j),gmtt(j) 
         end do


!     integration over theta
         do j = 1, 90
            
            thl1 = (real(j) - 1.0) * pi / 180.
            gl1 = fgl(thl1, m) / sin(max(thl1,1.0d-8))
            thl2 = thl1 +  pi / 180.
            gl2 = fgl(thl2, m) / sin(max(thl2,1.0d-8))

            gmrp(i) = gmrp(i) + (gl2 * gmrt(j+1) * sin(thl2) 
     &           + gl1 * gmrt(j) * sin(thl1)) * pi / (2.0 * 2.0 * 90.0)
            gmtp(i) = gmtp(i) + (gl2 * gmtt(j+1) * sin(thl2) 
     &           + gl1 * gmtt(j) * sin(thl1)) * pi / (2.0 * 2.0 * 90.0)

c            write(*,*) gmrp(i),gmtp(i),gl1,gl2

         end do
      end do

!     integration over phi

      do i = 1, 360
         gmr = gmr + (gmrp(i+1) + gmrp(i)) * 2.0 * pi / (2.0 * 360.0)
         gmt = gmt + (gmtp(i+1) + gmtp(i)) * 2.0 * pi / (2.0 * 360.0)
c         write(*,*) i,"flag",gmr,gmt,gmtp(i+1),gmtp(i)
      end do

      gmr = gmr / (2. * pi)
      gmt = gmt / (2. * pi)

      return
      end


!*************************************
!     calculation of alpha
!**************************************

      real function faf(th1, th2, ph)

      implicit none

c set local parameter
      real th1, th2, ph
      
      faf = cos(th1) * cos(th2) + sin(th1) * sin(th2) * cos(ph)
      faf = sign(min(abs(faf), 1. - 1.e-10),faf) 
      faf = acos(faf)

      return

      end
