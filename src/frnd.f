c******************************************************************
c     Random number function based on the Tausworthe sequence method
c     Original code and algorithm:
c     Shinomoto (1992), Statistical mechanics: Parity physics cource.,
c     Maruzen press, P27
c*******************************************************************
      real*8 function frnd()
      implicit none

      integer i,j,k,cnt
      integer nrand,lrand
      real*8 anorm
      common /rnd/lrand,cnt

      parameter (nrand=1000, anorm=1.0/2147483647)
      dimension lrand(-249:nrand)
      save /rnd/
c     to make the sequence many times, following is repeated.
      if(cnt.gt.nrand)then
         do i=1,nrand/100
            k=(i-1)*100
            do j=1+k,100+k
               lrand(j)=ieor(lrand(j-250),lrand(j-103))
            end do
         end do
         do j=-249,0
            lrand(j)=lrand(nrand+j)
         end do
         cnt=1
      end if

      frnd=lrand(cnt)*anorm
      cnt=cnt+1

      end
