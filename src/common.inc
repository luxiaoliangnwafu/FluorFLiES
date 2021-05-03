! space definition
      integer size,ndivs(720),divs(720,700),mdiv
      integer ixmax,iymax,izmax,nzmin
      parameter(size=300)
      real xmax,ymax,zmax,res,epsi,zmin
      common /space/xmax,ymax,zmax,ndivs,divs,res,epsi,mdiv,
     &      ixmax,iymax,izmax,zmin,nzmin

! atmospheric common
      integer knz, knkd, knchi, knxr, knyr, knrdc, nlut, nlut1, nang
      parameter (knz=200, knkd=3, knchi=10)
      parameter (knxr=size, knyr=size, knrdc=700)
      parameter (nlut=2000, nlut1=nlut+1)
      integer klayer(0:knz)
      real zgrd(0:knz), zgrd_back(0:knz), chigrd(100)
      real zgrdm(0:knz), trulut(6, knchi, knz)
      real absg1d(knz, knkd), extt1d(knz, knchi), abst1d(knz)
      real pflut(0:nlut1, knz), salut(0:nlut1, knz)
      real uxrtab(knrdc), uyrtab(knrdc), uzrtab(knrdc)
      double precision prdcF(knxr, knyr, knrdc),prdcQ(knxr, knyr, knrdc)
      real wrr, fsang
      integer nz, nchi, nrdc
      common /catmosppi/prdcF,prdcQ, pflut, salut, fsang, trulut,
     & absg1d, extt1d, abst1d, zgrd, chigrd, uxrtab, uyrtab, uzrtab,
     & wrr, nz, nchi, nrdc, nang, klayer, zgrd_back, zgrdm
      save /catmosppi/

! canopy common

      integer nangc, nth,nph,nobj,iobj(6000),sobj(6000)
      integer diftype,dirflag,bound
      integer irefl(2,size,size)
      integer mc,mb,mf,nts
      real obj(6000,5), sbar(5)
      real angt(100),angp(100)
      real brf(2,700)
      real brfc(2,700),brfs(2,700),brff(2,700)
      real refl(2,size,size), feye(0:90,0:360), rfeye(0:90,0:360)
      real ap(size,size,100),apf(size,size),aps(size,size)
      real apd(size,size,100),apb(size,size,100),apfd(size,size)
      real apnp(100)
      real tfpr, cfpr, bfpr, ffpr, sfpr
      real ffdir(size,size),ffdif(size,size)
      real sfdir(size,size),sfdif(size,size)
      real gLAI,u(6),BAD(5),ir,fe,bp1,bp2
      real gtblc(0:180),gtblb(0:180),gtblf(0:180)
      real uxrc(knrdc), uyrc(knrdc), uzrc(knrdc)
      real rb
      parameter (rb = 0.5)
      common /canopy/diftype,irefl,brf,sbar,
     &  brfc,brfs,brff,refl,dirflag,ap,apf,aps,apd,apb,apfd,apnp,
     &  tfpr, cfpr, bfpr, ffpr, sfpr, feye, rfeye,
     &  ffdir,ffdif,sfdir,sfdif,bound,nth,nph,angt,angp,
     &  gLAI,u,ir,BAD,fe,nobj,iobj,sobj,obj,gtblc,gtblb,gtblf,
     &  bp1,bp2,
     &  uxrc, uyrc, uzrc, nangc, nts, mc, mb, mf

c scattering angle parameters

      real gmrc(19, 19, 37), gmtc(19, 19, 37)
      real gmrb(19, 19, 37), gmtb(19, 19, 37)
      real gmrf(19, 19, 37), gmtf(19, 19, 37)
      common /scatpahse/gmrc, gmtc, gmrb, gmtb, gmrf, gmtf

c for back trace
      integer nb,bg

      common /bug/nb, bg

! for running mode
      integer running_mode
      common /mode/running_mode

! collect Qin, QFin [400, 750]
! MBI, MBII, MFI, MFII, Qfs [640, 850]
      integer size_height
      parameter(size_height=50)
      real Qin(351), QFin(351), t_Q_QF(351)
      real*8 Qfs(211, 351), fbrf(2,700)
      real*8 Qfs_all(211)
      real WL_SPAN(100), areaTOC
      real MBI(211, 351), MBII(211, 351)
      real MFI(211, 351), MFII(211, 351)
      real WL_SIF_START, WL_SIF_END, CURRENT_WL
      real WL_PAR_START, WL_PAR_END
      real WL_START, WL_END
      real FQEI, FQEII, Fyeild(3000)
      real*8 F_scat(211), F_excited(211)
      real RF_COSQ0, RQ_COSQ0
      integer par0, par1, par2, WP_INDEX
      real WAIT_PHOTON(100000, 10)
      real PAR_LEAF(size,size,size_height)
      real HC, NAVO, LEAFAREA
      character*25 FILE_LABEL, FILE_LABEL_G
      character*25 FILE_O_LABEL
      real SIFarea(300,300), SIFplace(211, 3)
      real SIFplace_e(211,3)
      integer LEAFLIT(size,size,size_height)
      integer init_nscat
      real*8 F_leaflit(211), F_leafshade(211) 
      integer GRASS, C_G_SPECTRUM
      real gMBI(211, 351), gMBII(211, 351)
      real gMFI(211, 351), gMFII(211, 351)
      real PPFD_T, PPFD_G, PPFD_S, PPFD_F
      real*8 SIF_760_3D(size, size, size_height)
      real*8 SIF_680_3D(size, size, size_height)
      real*8 SIF_760_3D_e(size, size, size_height)
      real*8 SIF_680_3D_e(size, size, size_height)
      real*8 SIF_760_3D_s(size, size, size_height)
      real*8 SIF_680_3D_s(size, size, size_height)
      integer SIF_e_1st(size, size, size_height)
      integer SIF_o_pos(size, size, size_height)
      real*8 Qfs_e(211, 351), Qfs_s(211, 351)
      integer SIF_place_sca_num(211,3)
      real W_T, W_G, W_S, W_F
      real WF_T, WF_G, WF_S, WF_F

      common /sif/Qin, QFin, t_Q_QF, FQEI, FQEII,
     & WL_SIF_START, WL_SIF_END, WL_SPAN, CURRENT_WL,
     & WL_START, WL_END,
     & MBI, MBII, MFI, MFII, fbrf, Qfs, areaTOC,
     & Qfs_all,
     & Fyeild, F_scat, F_excited,
     & par0, par1, par2, WL_PAR_START, WL_PAR_END,
     & WAIT_PHOTON, WP_INDEX, PAR_LEAF,
     & HC, NAVO, LEAFAREA, FILE_LABEL, FILE_O_LABEL,
     & RF_COSQ0, RQ_COSQ0, SIFarea, SIFplace,
     & SIFplace_e,
     & LEAFLIT, F_leaflit, F_leafshade, init_nscat,
     & GRASS, gMBI, gMBII, gMFI, gMFII, FILE_LABEL_G,
     & PPFD_T, PPFD_G, PPFD_S, PPFD_F,
     & C_G_SPECTRUM, Qfs_e, Qfs_s, 
     & SIF_680_3D, SIF_760_3D,
     & SIF_680_3D_e, SIF_680_3D_s, 
     & SIF_760_3D_e, SIF_760_3D_s,
     & SIF_place_sca_num,
     & W_T, W_G, W_S, W_F,
     & WF_T, WF_G, WF_S, WF_F
