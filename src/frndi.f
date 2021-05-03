c******************************************************************
c     Initialized function of rondom number sequence
c     Random number function based on the Tausworthe sequence method 
c     Original code and algorithm:
c     Shinomoto (1992), Statistical mechanics: Parity physics cource.,
c     Maruzen press, P27
c*******************************************************************

      integer function frndi(ix)
      implicit none


      integer i,j,k,irep,cnt
      integer ix,nrand,lrand
      real*8 anorm
      common /rnd/lrand,cnt

      parameter (nrand=1000, anorm=1.0/2147483647)
      dimension lrand(-249:nrand)

      save /rnd/

c      ix=715327539
      
c----------------------------------------------------------
c     nrand: the number of random numbers which should be a 
c     multiple of 100, and larger than 300.
c     lrand: integer random number generated
c     anorm:normalization constant
c----------------------------------------------------------

c Here is a initial loop, this is used in only first time. 
      do i=-249,0
         ix=ix*48828125
         if(ix.lt.0) ix=(ix+2147483647)+1
         lrand(i)=ix
      end do

      do i=-102,-72
         lrand(i)=ibset(lrand(i),i+102)
         do j=0,i+101
            lrand(i)=ibclr(lrand(i),j)
         end do
      end do

      do irep=1,10
         do i=1,nrand/100
            k=(i-1)*100
            do j=1+k,100+k
               lrand(j)=ieor(lrand(j-250),lrand(j-103))
            end do
         end do
         do j=-249,0
            lrand(j)=lrand(nrand+j)
         end do
      end do
      
      cnt=1001

      frndi=0

      end
