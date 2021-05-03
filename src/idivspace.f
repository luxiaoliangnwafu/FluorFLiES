! ******************************
!     idivspace
!     Hideki Kobayashi
!     2008/01/15
!  ******************************

      subroutine idivspace
      implicit none

      include 'common.inc'
      include 'math.inc'

      integer i, j, k, n, l
      integer idiv, flag
      integer ix, iy, iz

      real x, y, z, intv
      real xr, yr, zu, zb
      real divx, divy, divz
      real a(4), b(4), c(6), cc(4), d, dd
      real rx, ry, rr
      real p1, p2, max, min

!     divided space interval
!     int is a 50-pixel-equivalent length
      intv = 50.0 / res

      flag = 0

!     total gridding

      ixmax = 6
      iymax = 6
      izmax = 20

      ix = 1
      iy = 1
      iz = 1

      idiv = ixmax * iymax * izmax
      write(*,*) "idiv=",idiv

!     set distance calculation parameters
      a(1) = 1.0
      a(2) = 1.0
      a(3) = 0.0
      a(4) = 0.0

      b(1) = 0.0
      b(2) = 0.0
      b(3) = 1.0
      b(4) = 1.0


!     initialize the div array
      do i=1, idiv
         ndivs(i) = 0
         do j=1, 300
            divs(i,j) = 0
         end do
      end do

!     start the idiv loop
      do i=1, idiv
         !write(*,*) "i=",i
         n = 1

!     preparation of the i-th grid for space divided method
         divx = (real(ix) - 1.0) * intv
         divy = (real(iy) - 1.0) * intv
         divz = (real(iz) - 1.0) * intv
         c(1) = divx
         c(2) = divx + intv
         c(3) = divy
         c(4) = divy + intv
         c(5) = divz
         c(6) = divz + intv

!     definition of rectangular of the i-th object
         do j=1, nobj
            !write(*,*) "j=",j
            flag = 0

            if(sobj(j) .eq. 1)then ! cone
               xr = obj(j,1)
               yr = obj(j,2)
               zu = obj(j,3)
               zb = obj(j,3) - obj(j,4)
            elseif(sobj(j) .eq. 2 .or. sobj(j) .eq. 4)then ! cylinder
               xr = obj(j,1)
               yr = obj(j,2)
               zu = obj(j,3)
               zb = obj(j,3) - obj(j,4)
            elseif(sobj(j) .eq. 3)then ! ellipsoid
               xr = obj(j,1)
               yr = obj(j,2)
               zu = obj(j,3) + obj(j,4)
               zb = obj(j,3) - obj(j,4)
            elseif(sobj(j) .eq. 5)then ! half ellipsoid
               xr = obj(j,1)
               yr = obj(j,2)
               zu = obj(j,3) + obj(j,4)
               zb = obj(j,3)
            end if

!     check the intersection on the x-y plane
            do k = 1, 4
               !write(*,*) "k=",k
               d = abs(xr * a(k) + yr * b(k) - c(k))
c     d=d/sqrt(a(k)*a(k)+b(k)*b(k))
               if(d .le. obj(j,5))then
                  dd = sqrt(obj(j,5) * obj(j,5) - d * d)
                  p1 = b(k) * xr + a(k) * yr - dd
                  p2 = b(k) * xr + a(k) * yr + dd
                  min = b(k) * divx + a(k) * divy
                  max = b(k) * divx + a(k) * divy + intv

                  if(min .le. p1 .and. max .ge. p1) flag = 1
                  if(min .le. p2 .and. max .ge. p2) flag = 1

               end if
            end do

            do k = 1, 2
               do l = 1, 2
                  rx = xr - c(k)
                  ry = yr - c(l + 2)
                  rr = sqrt(rx * rx + ry * ry)
                  if(rr .le. obj(j,5)) flag = 1
               end do
            end do

            if(flag .eq. 0)then
               if((xr .ge. c(1) .and. xr .le. c(2)) .and.
     &              (yr .ge. c(3) .and. yr .le. c(4))) then
                  flag = 1
              end if
            end if



!     check the intersection for z-axis
            if(flag .eq. 1)then
               if(zu .gt. c(5) .and. zu .le. c(6))then
                  flag = 2
c     write(*,*) i, zu,zb,c(5),c(6),"1"
             elseif(zb .ge. c(5) .and. zb .lt. c(6))then
                  flag = 2
c     write(*,*) i, zu,zb,c(5),c(6),"2"
               elseif(zu .ge. c(6) .and. zb .le. c(5))then
                  flag = 2
c     write(*,*) i, zu,zb,c(5),c(6),"3"
               end if
            end if


            !write(*,*) "flag=",flag
            !write(*,*) "n=",n
!     input data number for the ndivs & divs
            if(flag .eq. 2)then
               ndivs(i) = ndivs(i) + 1
               divs(i,n) = j
               n = n + 1
               mdiv = i
            end if
         end do

!     count
         ix = ix + 1
         if(ix .gt. ixmax)then
            ix = 1
            iy = iy + 1
            if(iy .gt. iymax)then
               iy = 1
               iz = iz + 1
            end if
         end if
      end do

!     determinatin of zmax at the boudary of the big voxel
      zmax = intv * (1.0 + real((mdiv - 1) / (ixmax * iymax)))

      open(12,file = "div.out")
      write(12,*) zmax
      write(12,*) mdiv
      write(12,*) idiv
      do i=1, idiv
         write(12,*) i,(divs(i,j),j = 1, ndivs(i))
      end do
      close(12)
      end
