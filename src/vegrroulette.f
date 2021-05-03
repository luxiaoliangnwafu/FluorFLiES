c     Russian roulette


      subroutine vegrroulette(w,e)
      implicit none

      real w,e,q
      real*8 frnd
      real rnd, conv

      ! for SIF process, when init wf, wf = 0
      ! or w extinct, just pass
      conv = 1.d-8
      if (w .lt. conv) return

c q is a probability of survival
      q=2.0*e
c      q=0.5

      if(w.lt.e) then
         rnd=real(frnd())

         if(rnd.le.q) then
            w=w/q
         else
            w=1.d-9
         end if
      end if

      return
      end
