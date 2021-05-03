!*************************************************************************
! calculate SIF radiance on leaf level
!
! *INPUT VARIABLES:
!
! photon position: x, y, z
! incident vector: vix, viy, viz
! leaf normal vector: vlx, vly, vlz
! observe vector: vox, voy, voz
! ref, tr means relfectance and transmittance
! bp: 1 means branch, 0 means canopy,
!     when bp = 1, no transmittance occurs
! nfluor means scattering times of fluorescence
! place = 1 means crown, place = 2 means floor
!
!                             written Sicong Gao
!                             last modified 26/02/18
!**************************************************************************
      subroutine fluorrad(wf, w, wq, apar, wl,
     &                    x, y, z,
     &                    vix, viy, viz, vlx, vly, vlz, vox, voy, voz,
     &                    ref, tr, ssa,
     &                    nfluor, bp, place,
     &                    cb, ikd, ichi, fd)

      implicit none
      include "common.inc"

      ! input
      real x, y, z, vix, viy, viz, vlx, vly, vlz, vox, voy, voz
      real wf, w, wq, apar, wl
      real ref, tr, ssa, bp
      integer nfluor, place, cb, ikd, ichi
      real fd

      ! work
      real vec, flagfscat
      real excitLight, calEFMatrix
      integer wlid, sflag, i, ii
      real tauc, xt, yt, zt, uxt, uyt, uzt
      real fwl, rad, conv, ilight
      real tw, vflag
      integer tnfluor, tnscat, conv_t, scatType
      integer ix, iy, iz

      ! write(*,*) "current_iwl = ", CURRENT_IWL
      ! retrieve fwl

      ! CONDITION:
      ! 1. 400 - 750 nm (PAR) can emitte the fluorescence
      ! 2. The fluorescence range is in 640 - 850 nm

	    ! [400, 750] process excitLight
	    ! [640, 850] process scatter fluorescence
      if ((wl .lt. WL_PAR_START) .or. (wl .gt. WL_SIF_END)) then
          return
      end if

      if (running_mode .eq. 0) return

      conv = 1.d-8

      xt = x - (aint(x / xmax) - 0.5 + sign(0.5, x)) * xmax
      yt = y - (aint(y / ymax) - 0.5 + sign(0.5, y)) * ymax
      zt = z

      ix = int(xt * res) + 1
      iy = int(yt * res) + 1
      iz = int(zt) + 1

      ix = min(ix, size)
      iy = min(iy, size)

	    !###############################################
	    ! scattering
	    !###############################################
      if (((wl.ge.WL_SIF_START).and.(wl.le.WL_SIF_END))) then

         wlid = int((wl - WL_SIF_START) * 1000) + 1

         call LEM(wf, x, y, z, vix, viy, viz,
     &            ref, tr, cb, fd, ichi, ikd, 0, wlid, nfluor, place)

      end if

	    !###############################################
	    ! emission
	    !###############################################
      fwl = WL_SIF_START

      do while (fwl .le. WL_SIF_END .and. bp .eq. 0)
		     ! WAITTED: consider par is small enough

         ! scatType = 0 means backward, 1 means forward
         do scatType = 0, 1
		        ! excited fluorescence
            excitLight = calEFMatrix(wl, fwl, w, wf, apar,
     &                            scatType, xt, yt, zt, place)
            vec = abs(vix * vlx + viy * vly + viz * vlz)
            excitLight = excitLight * vec * (1.0 - bp)

            wlid = int((fwl - WL_SIF_START) * 1000) + 1
            Qfs_all(wlid) = Qfs_all(wlid) + excitLight

            ! generate the out direction for excited fluorescence
            call getscattervec(vix, viy, viz, vlx, vly, vlz,
     &                         uxt, uyt, uzt, scatType)
            
            ! single simulate fluorescence, push to stack
            if (excitLight .gt. conv) then
                ! WAIT_PHOTON[x, y, z, ux, uy, uz, wq, wl]
                WP_INDEX = WP_INDEX + 1
                WAIT_PHOTON(WP_INDEX, 1) =  xt
                WAIT_PHOTON(WP_INDEX, 2) =  yt
                WAIT_PHOTON(WP_INDEX, 3) =  zt
                WAIT_PHOTON(WP_INDEX, 4) =  uxt
                WAIT_PHOTON(WP_INDEX, 5) =  uyt
                WAIT_PHOTON(WP_INDEX, 6) =  uzt
                WAIT_PHOTON(WP_INDEX, 7) =  excitLight
                WAIT_PHOTON(WP_INDEX, 8) =  fwl
                WAIT_PHOTON(WP_INDEX, 9) =  real(ichi)
                WAIT_PHOTON(WP_INDEX, 10) =  real(ikd)

                par0 = par0 + 1
            end if

			      ! multi-angle observation
            call LEM(excitLight, x, y, z, vix, viy, viz,
     &               ref, tr, cb, fd, ichi, ikd, 1, wlid, nfluor,place)

          end do  !do scatType = 0, 1

          fwl = fwl + 0.001
        end do ! do while (fwl .le. WL_SIF_END)
        nfluor = nfluor + 1
      end

!***********************************************************************
! calculate EF-matrix with incident light
!
! *INPUT VARIABLES:
!
! apar wavelength: par_wl
! excited fluorescence wavelength: f_wl
! incident non-fluorescence light: q_in
! incident fluorescence light: qf_in
! scatType: 0 means backward, 1 means forward
!
!                             written Sicong Gao
!                             last modified 28/02/18
!***********************************************************************
      real function calEFMatrix(par_wl, f_wl, q_in, qf_in, apar,
     &                          scatType, x, y, z, place)

      implicit none
      include "common.inc"

      ! input
      real par_wl, f_wl, apar, q_in, qf_in, x, y, z
      integer scatType, place

      ! work
      integer i, par_wlid, f_wlid
      real sumvec, psivec, psiivec, f_efficient
      real conv, tMBI, tMBII, tMFI, tMFII

      !write(*,*) "par, q_in, qf_in=", par, q_in, qf_in
      conv = 1.d-8

      if (apar .lt. conv)  then
          calEFMatrix = 0.0
          return
      end if

	    ! f_wl is guranteed by upper level
      if ((par_wl .lt. WL_PAR_START).or.(par_wl .gt. WL_PAR_END)) then
          calEFMatrix = 0.0
        return
      end if

      sumvec = 0.0
      psivec = 0.0
      psiivec = 0.0

      par_wlid = int(par_wl * 1000) - int(WL_PAR_START * 1000) + 1
      f_wlid = int(f_wl * 1000) - int(WL_SIF_START * 1000) + 1
      
      if ((place .eq. 2) .and. (GRASS .eq. 1)) then
          tMBI = gMBI(f_wlid, par_wlid)
          tMBII = gMBII(f_wlid, par_wlid)
          tMFI = gMFI(f_wlid, par_wlid)
          tMFII = gMFII(f_wlid, par_wlid)  
          !write(*,*) "IN grass EF-M"
      else
          tMBI = MBI(f_wlid, par_wlid)
          tMBII = MBII(f_wlid, par_wlid)
          tMFI = MFI(f_wlid, par_wlid)
          tMFII = MFII(f_wlid, par_wlid)             
      endif

      if (scatType .eq. 0) then
        psivec = apar * tMBI
        psiivec = apar * tMBII
      else
        psivec = apar * tMFI
        psiivec = apar * tMFII
      end if
      
      call getFluoyieldPar(apar, f_efficient)
      !f_efficient = 1.0
      sumvec = psivec + psiivec * f_efficient

      calEFMatrix = sumvec

      end

      ! get fluorescence yield depends on the par
      subroutine getFluoyieldPos(x, y, z, f_efficient)
      implicit none
      include "common.inc"

      ! input
      real x, y, z, f_efficient

      ! work
      integer i_par, ix, iy, iz

      ix = int(x * res) + 1
      iy = int(y * res) + 1
      iz = int(z) + 1

      ix = min(ix, size)
      iy = min(iy, size)

      i_par = int(PAR_LEAF(ix, iy, iz))

      !par1 = par1 + 1
      !write(*,*) "G_IRRADIANCE = ", G_IRRADIANCE
      !write(*,*) "par = ", par
      if (i_par .gt. 3000) then
         !write(*,*) "ERROR: par is too large. - getFluoyield ()"
         !write(*,*) "ipar = ", i_par
         f_efficient = Fyeild(3000)
         !par0 = par0 + 1
      else if (i_par .le. 0) then
         !write(*,*) "ERROR: par is too small. - getFluoyield ()"
         !write(*,*) "ipar = ", i_par
         !par0 = par0 + 1
         f_efficient = 0
      else
         f_efficient = Fyeild(i_par)
      end if

      end

      subroutine getFluoyieldPar(par, f_efficient)
      implicit none
      include "common.inc"

      ! input
      real par, f_efficient

      ! work
      real i_par

      i_par =int(par * RQ_COSQ0)
      
      if (i_par .ne. 0) then 
          !write(*,*) "ipar = ", i_par
      end if
      if (i_par .gt. 3000) then
         f_efficient = Fyeild(3000)
      else if (i_par .le. 0) then
         f_efficient = 1
      else
         f_efficient = Fyeild(i_par)
      end if

      end

      subroutine calParLeaf(x, y, z, wl, qincident)
      implicit none
      include "common.inc"

      ! input
      real wl, qincident
      integer x, y, z

      ! work
      integer iwl

      if (running_mode .eq. 1) return

      iwl = int(wl * 1000)
      PAR_LEAF(x, y, z) = PAR_LEAF(x, y, z) + qincident * iwl

      end

      !subroutine collectFluorRad()

      ! check the sunlit sunshade area
      subroutine calSunLeaf(x, y, z, nscat)
      implicit none
      include "common.inc"

      integer x, y, z, nscat, t_nscat
      
      ! LEAFLIT(x,y,z) = 0 : no sun touched leaf or other elements
      ! LEAFLIT(x,y,z) = 1 : sun shade
      ! LEAFLIT(x,y,z) = 2 : sunlit
      if (running_mode .eq. 1) return

      t_nscat = t_nscat - init_nscat
      if (t_nscat .eq. 0) then
        LEAFLIT(x, y, z) = 2 
      else
        ! if the leaf has been sunlit, cannot set the leaf as shade
        if (LEAFLIT(x, y, z) .ne. 2)  LEAFLIT(x, y, z) = 1
      end if

      end


      subroutine LEM(wf, x, y, z, ux, uy, uz,
     &           lr, lt, cb, fd, ichi, ikd, e_s, wlid, nfluor, place)

      implicit none

      include 'common.inc'
      include 'math.inc'

!     cb:cb = 1 (overstory),cb = 2 (branch)
!     cb = 3 (stem top/bottom),cb = 4 (understory), cb = 5 (soil surface)
!     cb = 6 (stem side)
!     a is only used for lambertian reflection from the stem side
!     other case "a" should be 1.0
!     e_s = 1 means emit, 0 means scatter

      integer i, j, cb
      integer ichi, ikd, ith, ithr
      integer sflag, e_s, wlid, nfluor, place
      integer ix, iy, iz

      real wf, x, y, z, ux, uy, uz
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
      if (wf .lt. conv) return

      Id = 0.0
      a = 1.0

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

      ix = int(x * res) + 1
      iy = int(y * res) + 1
      iz = int(z) + 1

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
     &        * gtblf(ithr) * gLAI * (dlt(cb, 4) + dlt(cb, 5))
         tauc = tauc * hk

         if(cb .eq. 6) then
            rr = sqrt(uxrc(i) * uxrc(i) + uyrc(i) * uyrc(i))
            rr = max(1.d-3, rr)
            ph = uxrc(i) / rr
            ph = facos(ph)
            a = abs(fsin(facos(uzrc(i))) * fcos(a - ph))
         end if

         Id = 0.0
         Id = ff(cb) * wf * pf * fexp(-tauc) / fcos(thr)
         Id = Id + (1.- ff(cb)) * wf * fexp(-tauc) / pi
         Id = Id * abs(a)
         ! floor, stem, branch: fd = 0
         Id = Id * (1. - fd)
         ! sflag = 0 means stem collision
         Id = Id * real(sflag)

         
         if (e_s .eq. 1) then
            Qfs(wlid, i) = Qfs(wlid, i) + Id
            Qfs_e(wlid, i) = Qfs_e(wlid, i) + Id
            SIFplace_e(wlid, place) = SIFplace_e(wlid, place) + Id
         else
            Qfs(wlid, i) = Qfs(wlid, i) + Id
            Qfs_s(wlid, i) = Qfs_s(wlid, i) + Id
            !F_scat(wlid) = F_scat(wlid) + Id
         end if
         SIFplace(wlid, place) = SIFplace(wlid, place) + Id
         SIF_place_sca_num(wlid, place) = 
     &          SIF_place_sca_num(wlid, place) + 1

         if (wlid .eq. 121 .and. i .eq. 1) then
             SIF_760_3D(ix,iy,iz) = SIF_760_3D(ix,iy,iz) + Id
            
             if (e_s .eq. 1) then
               SIF_760_3D_e(ix,iy,iz) = SIF_760_3D_e(ix,iy,iz) + Id
            else
               SIF_760_3D_s(ix,iy,iz) = SIF_760_3D_s(ix,iy,iz) + Id
            end if
         end if

         if (wlid .eq. 41 .and. i .eq.1) then
            SIF_680_3D(ix,iy,iz) = SIF_680_3D(ix,iy,iz) + Id
            
            if (e_s .eq. 1) then
               SIF_680_3D_e(ix,iy,iz) = SIF_680_3D_e(ix,iy,iz) + Id
            else
               SIF_680_3D_s(ix,iy,iz) = SIF_680_3D_s(ix,iy,iz) + Id
            end if
         end if

      end do
         
      if (nfluor .eq. 0 .and. e_s .eq. 1) then
          SIF_e_1st(ix,iy,iz) = 1
      end if
        
      if (Id .ne. 0.0) then
          SIF_o_pos(ix,iy,iz) = SIF_o_pos(ix,iy,iz) + 1
      end if

      return
      end