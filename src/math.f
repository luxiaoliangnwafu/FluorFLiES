!     *** Math function ***

!     **** sine **** 
      real function fsin(x)
      implicit none
      include 'math.inc'

      real x

      fsin = tsin(int(x * 1.d4))
      
      end

!     **** cos ****
      real function fcos(x)
      implicit none
      include 'math.inc'

      real x

      fcos = tcos(int(x * 1.d4))
      
      end

!     **** acos ****
      real function facos(x)
      implicit none
      include 'math.inc'

      real x
      
      x = min(x, 1.0)
      x = max(x, -1.0)

!      if((x.ge.-1.0).and.(x.le.1.0)) then
         facos = tacos(int(x * 1.d4))
!      else 
!         write(*,*) "Invalid facos input",x
!      end if

      end


!     **** exp ****
!     x shuould be less than 0
      real function fexp(x)
      implicit none
      include 'math.inc'

      integer i
      real x, y
      
      x = max(-80.0, x)
      y = -125. * x
      i = int(y)
      
      fexp = texp(i) * (real(i + 1) - y) 
     &     + texp(1 + i) * (y - real(i))
c      fexp = fexp 

      if(x .gt. 0.0) then
         write(*, *) "Invalid exp(-x)",x
      end if

      end
   
   
