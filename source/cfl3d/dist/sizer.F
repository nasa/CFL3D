      subroutine sizer(mwork,mworki,nplots,minnode,iitot,intmax,
     .                 maxxe,mxbli,nsub1,lbcprd,lbcemb,lbcrad,
     .                 maxbl,maxgr,maxseg,maxcs,ncycmax,intmx,mxxe,
     .                 mptch,msub1,nmds,maxaes,maxsegdg,nnodes,nslave,
     .                 nmaster,myhost,myid,mycomm,nplots0,maxnode0,
     .                 mxbli0,lbcprd0,lbcemb0,lbcrad0,maxbl0,maxgr0,
     .                 maxseg0,maxcs0,ncycmax0,intmax0,nsub10,intmx0,
     .                 mxxe0,mptch0,msub10,ibufdim0,nbuf0,mxbcfil0,
     .                 nmds0,maxaes0,maxsegdg0,imode,ntr,bcfiles,
     .                 bou,nou)
c
c     $Id$
c
c***********************************************************************
c     Purpose: reads in a cfl3d input file (plus patch/ovelap files, if
c     appropriate) and determines the array sizes required for the
c     particular problem in the cfl3d input file.
c
c     for parallel mode, memory is allocated identically on all
c     processors, based on the maximum requirement on any of the noodes.
c
c     imode governs whether the routine is being used a a stand-alone
c     preprocessor, or as part of the dynamic memory allocation in 
c     cfl3d: 
c    
c     imode = 0  stand-alone preprocessor
c             1  part of the dynamic memory allocation in cfl3d
c***********************************************************************
c
c     parameter definitions
c
c     maxnode0 - maximum number of cpus for parallel processing
c     maxgr0   - maximum number of grids
c     maxbl0   - maximum number of blocks = maxg*( 1 + ncg )
c     maxseg0  - maximum number of segments specified per face
c     nplots0  - maximum number of data sets to output via plot3d or 
c                print options
c     mxbli0   - max no of block interps for matching bdys, 
c                including coarser levels
c     ncycmax0 - maximum number of time-steps (cycles)
c     intmx0   - maximum number of patch interfaces
c     msub10   - maximum number of "from" blocks on a patch interface
c     ibufdim0 - maximum internal buffer (array) size for storing 
c                output data
c     nbuf0    - maximum number of internal buffers
c     mxbcfil0 - maximum number of auxiliary files for setting 2000
c                series bc's
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*80 grid,plt3dg,plt3dq,output,residual,turbres,blomx,
     .             output2,printout,pplunge,ovrlap,patch,restrt,
     .             subres,subtur,grdmov,alphahist,errfile,preout,
     .             aeinp,aeout,sdhist
      character*80 bcfiles(mxbcfil0)
      character*120 bou(ibufdim0,nbuf0)
c
      integer xif1,etf1,xif2,etf2
      integer bcfilei,bcfilej,bcfilek
c
      dimension bcfilei(maxbl0,maxseg0,2),bcfilej(maxbl0,maxseg0,2),
     .          bcfilek(maxbl0,maxseg0,2)
      dimension mglevg(maxbl0),nemgl(maxbl0),ipl3dtmp(11*nplots0)
      dimension ifrom(msub10),xif1(msub10),etf1(msub10),xif2(msub10),
     .          etf2(msub10)
      dimension nou(nbuf0)
      dimension memblock(maxbl0),iv(maxbl0),icsinfo(maxcs0,9)
      dimension nblfine(maxbl0)
      dimension lw(65,maxbl0),lw2(43,maxbl0)
      dimension ltot(maxbl0)
      dimension iwrk(maxbl0,3)
      dimension jjmax1(nsub10),kkmax1(nsub10)
      dimension iiint1(nsub10),iiint2(nsub10)
      dimension n14(maxbl0)
      dimension mblk2nd(maxbl0)
      dimension mem_req_node(maxnode0),no_of_points(maxnode0)
      dimension dx(intmx0,msub10),dy(intmx0,msub10),
     .          dz(intmx0,msub10),dthetx(intmx0,msub10),
     .          dthety(intmx0,msub10),dthetz(intmx0,msub10)
      dimension inpr(nplots0,11)
      dimension inpl3d(nplots0,11)
      dimension lig(maxbl0),lbg(maxbl0)
      dimension ibpntsg(maxbl0,4),iipntsg(maxbl0)
      dimension lwdat(maxbl0,maxseg0,6)
      dimension nblk(2,mxbli0),limblk(2,6,mxbli0),
     .          isva(2,2,mxbli0),nblon(mxbli0)
      dimension iovrlp(maxbl0)
      dimension rkap0g(maxbl0,3),levelg(maxbl0),igridg(maxbl0),
     .          iflimg(maxbl0,3),ifdsg(maxbl0,3),iviscg(maxbl0,3),
     .          jdimg(maxbl0),kdimg(maxbl0),idimg(maxbl0),
     .          idiagg(maxbl0,3),nblcg(maxbl0),idegg(maxbl0,3),
     .          jsg(maxbl0),ksg(maxbl0),isg(maxbl0),jeg(maxbl0),
     .          keg(maxbl0),ieg(maxbl0),mit(5,maxbl0),
     .          ilamlog(maxbl0),ilamhig(maxbl0),jlamlog(maxbl0),
     .          jlamhig(maxbl0),klamlog(maxbl0),klamhig(maxbl0),
     .          iwfg(maxbl0,3)
      dimension utrans(maxbl0),vtrans(maxbl0),wtrans(maxbl0),
     .          omegax(maxbl0),omegay(maxbl0),omegaz(maxbl0),
     .          xorig(maxbl0),yorig(maxbl0),zorig(maxbl0),
     .          dxmx(maxbl0),dymx(maxbl0),dzmx(maxbl0),
     .          dthxmx(maxbl0),dthymx(maxbl0),dthzmx(maxbl0),
     .          thetax(maxbl0),thetay(maxbl0),thetaz(maxbl0),
     .          rfreqt(maxbl0),rfreqr(maxbl0),xorig0(maxbl0),
     .          yorig0(maxbl0),zorig0(maxbl0),time2(maxbl0),
     .          thetaxl(maxbl0),thetayl(maxbl0),thetazl(maxbl0),
     .          itrans(maxbl0),irotat(maxbl0),idefrm(maxbl0)
      dimension bcvali(maxbl0,maxseg0,7,2),bcvalj(maxbl0,maxseg0,7,2),
     .          bcvalk(maxbl0,maxseg0,7,2),nbci0(maxbl0),
     .          nbcidim(maxbl0),nbcj0(maxbl0),nbcjdim(maxbl0),
     .          nbck0(maxbl0),nbckdim(maxbl0),
     .          ibcinfo(maxbl0,maxseg0,7,2),
     .          jbcinfo(maxbl0,maxseg0,7,2),
     .          kbcinfo(maxbl0,maxseg0,7,2)
      dimension ncgg(maxgr0),nblg(maxgr0),iemg(maxgr0),
     .          inewgg(maxgr0)
      dimension iadvance(maxbl0),iforce(maxbl0)
      dimension iindx(intmx0,6*msub10+9),iindex(intmax0,2*nsub10+9),
     .          llimit(intmx0),iitmax(intmx0),mmcxie(intmx0),
     .          mmceta(intmx0),ncheck(maxbl0),iifit(intmx0),
     .          mblkpt(mxxe0),iic0(intmx0),iiorph(intmx0),
     .          iitoss(intmx0),ifiner(intmx0)
      dimension isav_pat(intmax0,17),isav_pat_b(intmax0,nsub10,6)
      dimension isav_dpat(intmx0,17),
     .          isav_dpat_b(intmx0,msub10,6)
      dimension isav_blk(2*mxbli0,17)
      dimension isav_prd(lbcprd0,12)
      dimension isav_emb(lbcemb0,12)
      dimension utrnsae(maxbl0,maxsegdg0),vtrnsae(maxbl0,maxsegdg0),
     .          wtrnsae(maxbl0,maxsegdg0),omgxae(maxbl0,maxsegdg0),
     .          omgyae(maxbl0,maxsegdg0),omgzae(maxbl0,maxsegdg0),
     .          xorgae(maxbl0,maxsegdg0),yorgae(maxbl0,maxsegdg0),
     .          zorgae(maxbl0,maxsegdg0),thtxae(maxbl0,maxsegdg0),
     .          thtyae(maxbl0,maxsegdg0),thtzae(maxbl0,maxsegdg0),
     .          rfrqtae(maxbl0,maxsegdg0),rfrqrae(maxbl0,maxsegdg0)
      dimension xorgae0(maxbl0,maxsegdg0),yorgae0(maxbl0,maxsegdg0),
     .          zorgae0(maxbl0,maxsegdg0),icouple(maxbl0,maxsegdg0)
      dimension icsi(maxbl0,maxsegdg0),icsf(maxbl0,maxsegdg0),
     .          jcsi(maxbl0,maxsegdg0),jcsf(maxbl0,maxsegdg0),
     .          kcsi(maxbl0,maxsegdg0),kcsf(maxbl0,maxsegdg0)
      dimension nsegdfrm(maxbl0),idfrmseg(maxbl0,maxsegdg0),
     .          iaesurf(maxbl0,maxsegdg0)
      dimension freq(nmds0,maxaes0),gmass(nmds0,maxaes0),
     .          gf0(2*nmds0,maxaes0),damp(nmds0,maxaes0),
     .          x0(2*nmds0,maxaes0),perturb(nmds0,maxaes0,4)
      dimension aesrfdat(5,maxaes0)
      dimension iskip(maxbl0),jskip(maxbl0),kskip(maxbl0)
c
      pointer (ip_bcfilei,bcfilei),(ip_bcfilej,bcfilej),
     .        (ip_bcfilek,bcfilek)
      pointer (ip_ifrom,ifrom),(ip_xif1,xif1),(ip_etf1,etf1),
     .        (ip_xif2,xif2),(ip_etf2,etf2)
      pointer (ip_icsinfo,icsinfo),(ip_iwrk,iwrk)
      pointer (ip_jjmax1,jjmax1),(ip_kkmax1,kkmax1),(ip_iiint1,iiint1),
     .        (ip_iiint2,iiint2),(ip_mglevg,mglevg),(ip_nemgl,nemgl),
     .        (ip_ipl3dtmp,ipl3dtmp)
      pointer (ip_mem_req_node,mem_req_node),(ip_nblfine,nblfine)
      pointer (ip_no_of_points,no_of_points),(ip_memblock,memblock)
      pointer (ip_lwdat,lwdat),(ip_iv,iv)
      pointer (ip_lw,lw),(ip_lw2,lw2)
      pointer (ip_nblk,nblk),(ip_limblk,limblk),
     .        (ip_isva,isva),(ip_nblon,nblon)
      pointer (ip_lig,lig),(ip_lbg,lbg),(ip_iovrlp,iovrlp)
      pointer (ip_ibpntsg,ibpntsg),(ip_iipntsg,iipntsg)
      pointer (ip_rkap0g,rkap0g),(ip_levelg,levelg),(ip_igridg,igridg),
     .        (ip_iflimg,iflimg),(ip_ifdsg,ifdsg),(ip_iviscg,iviscg),
     .        (ip_jdimg,jdimg),(ip_kdimg,kdimg),(ip_idimg,idimg),
     .        (ip_idiagg,idiagg),(ip_nblcg,nblcg),(ip_idegg,idegg),
     .        (ip_jsg,jsg),(ip_ksg,ksg),(ip_isg,isg),(ip_jeg,jeg),
     .        (ip_keg,keg),(ip_ieg,ieg),(ip_mit,mit),
     .        (ip_ilamlog,ilamlog),(ip_ilamhig,ilamhig),
     .        (ip_jlamlog,jlamlog),(ip_jlamhig,jlamhig),
     .        (ip_klamlog,klamlog),(ip_klamhig,klamhig),(ip_iwfg,iwfg)
      pointer (ip_utrans,utrans),(ip_vtrans,vtrans),(ip_wtrans,wtrans),
     .        (ip_omegax,omegax),(ip_omegay,omegay),(ip_omegaz,omegaz),
     .        (ip_xorig,xorig),(ip_yorig,yorig),(ip_zorig,zorig),
     .        (ip_dxmx,dxmx),(ip_dymx,dymx),
     .        (ip_dzmx,dzmx),(ip_dthxmx,dthxmx),(ip_dthymx,dthymx),
     .        (ip_dthzmx,dthzmx),(ip_thetax,thetax),(ip_thetay,thetay),
     .        (ip_thetaz,thetaz),(ip_rfreqt,rfreqt),
     .        (ip_rfreqr,rfreqr),(ip_xorig0,xorig0),(ip_yorig0,yorig0),
     .        (ip_zorig0,zorig0),(ip_time2,time2),(ip_thetaxl,thetaxl),
     .        (ip_thetayl,thetayl),(ip_thetazl,thetazl),
     .        (ip_itrans,itrans),(ip_irotat,irotat),(ip_idefrm,idefrm)
      pointer (ip_bcvali,bcvali),(ip_bcvalj,bcvalj),
     .        (ip_bcvalk,bcvalk),(ip_nbci0,nbci0),
     .        (ip_nbcidim,nbcidim),(ip_nbcj0,nbcj0),
     .        (ip_nbcjdim,nbcjdim),(ip_nbck0,nbck0),
     .        (ip_nbckdim,nbckdim),(ip_ibcinfo,ibcinfo),
     .        (ip_jbcinfo,jbcinfo),(ip_kbcinfo,kbcinfo)
      pointer (ip_ncgg,ncgg),(ip_nblg,nblg),(ip_iemg,iemg),
     .        (ip_inewgg,inewgg)
      pointer (ip_inpl3d,inpl3d),(ip_inpr,inpr)
      pointer (ip_iadvance,iadvance),(ip_iforce,iforce)
      pointer (ip_iindx,iindx),(ip_iindex,iindex),
     .        (ip_llimit,llimit),(ip_iitmax,iitmax),(ip_mmcxie,mmcxie),
     .        (ip_mmceta,mmceta),(ip_ncheck,ncheck),(ip_iifit,iifit),
     .        (ip_mblkpt,mblkpt),(ip_iic0,iic0),(ip_iiorph,iiorph),
     .        (ip_iitoss,iitoss),(ip_ifiner,ifiner)
      pointer (ip_dx,dx),(ip_dy,dy),(ip_dz,dz),
     .        (ip_dthetx,dthetx),(ip_dthety,dthety),
     .        (ip_dthetz,dthetz)
      pointer (ip_isav_blk,isav_blk)
      pointer (ip_isav_prd,isav_prd)
      pointer (ip_isav_pat,isav_pat),(ip_isav_pat_b,isav_pat_b)
      pointer (ip_isav_dpat,isav_dpat),(ip_isav_dpat_b,isav_dpat_b)
      pointer (ip_isav_emb,isav_emb)
      pointer (ip_mblk2nd,mblk2nd),(ip_n14,n14),(ip_ltot,ltot)
      pointer (ip_aesrfdat,aesrfdat),(ip_utrnsae,utrnsae),
     .        (ip_vtrnsae,vtrnsae),(ip_wtrnsae,wtrnsae),
     .        (ip_omgxae,omgxae),(ip_omgyae,omgyae),(ip_omgzae,omgzae),
     .        (ip_xorgae,xorgae),(ip_yorgae,yorgae),(ip_zorgae,zorgae),
     .        (ip_xorgae0,xorgae0),(ip_yorgae0,yorgae0),
     .        (ip_zorgae0,zorgae0),(ip_icouple,icouple),
     .        (ip_thtxae,thtxae),(ip_thtyae,thtyae),(ip_thtzae,thtzae),
     .        (ip_rfrqtae,rfrqtae),(ip_rfrqrae,rfrqrae),(ip_icsi,icsi),
     .        (ip_icsf,icsf),(ip_jcsi,jcsi),(ip_jcsf,jcsf),
     .        (ip_kcsi,kcsi),(ip_kcsf,kcsf),(ip_freq,freq),
     .        (ip_gmass,gmass),(ip_damp,damp),(ip_x0,x0),(ip_gf0,gf0),
     .        (ip_perturb,perturb),(ip_iskip,iskip),(ip_jskip,jskip),
     .        (ip_kskip,kskip),(ip_nsegdfrm,nsegdfrm),
     .        (ip_idfrmseg,idfrmseg),(ip_iaesurf,iaesurf)
c
      common /params/ lmaxgr,lmaxbl,lmxseg,lmaxcs,lnplts,lmxbli,lmaxxe,
     .                lnsub1,lintmx,lmxxe,liitot,isum,lncycm,
     .                isum_n,lminnode,isumi,isumi_n,lmptch,
     .                lmsub1,lintmax,libufdim,lnbuf,llbcprd,
     .                llbcemb,llbcrad,lnmds,lmaxaes,lnslave,lmxsegdg,
     .                lnmaster
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /ginfo/ jdim,kdim,idim,jj2,kk2,ii2,nblc,js,ks,is,je,ke,ie,
     .        lq,lqj0,lqk0,lqi0,lsj,lsk,lsi,lvol,ldtj,lx,ly,lz,lvis,
     .        lsnk0,lsni0,lq1,lqr,lblk,lxib,lsig,lsqtq,lg,
     .        ltj0,ltk0,lti0,lxkb,lnbl,lvj0,lvk0,lvi0,lbcj,lbck,lbci,
     .        lqc0,ldqc0,lxtbi,lxtbj,lxtbk,latbi,latbj,latbk,
c    .        lbcdj,lbcdk,lbcdi,lxib2,lux,lcmuv,lvolj0,lvolk0,lvoli0
     .        lbcdj,lbcdk,lbcdi,lxib2,lux,lcmuv,lvolj0,lvolk0,lvoli0,
     .        lxmdj,lxmdk,lxmdi,lvelg,ldeltj,ldeltk,ldelti,
     .        lxnm2,lynm2,lznm2,lxnm1,lynm1,lznm
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /maxiv/ ivmx
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /reyue/ reue,tinf,ivisc(3)
      common /unst/ time,cfltau,ntstep,ita,iunst
      common /twod/ i2d
      common /is_blockbc/ is_blk(5),ie_blk(5),ivolint
      common /is_perbc/ is_prd(5),ie_prd(5),nbcprd
      common /is_patch/ is_pat(5),ie_pat(5)
      common /is_embedbc/ is_emb(5),ie_emb(5),nbcemb
      common /filenam/ grid,plt3dg,plt3dq,output,residual,turbres,blomx,
     .                 output2,printout,pplunge,ovrlap,patch,restrt,
     .                 subres,subtur,grdmov,alphahist,errfile,preout,
     .                 aeinp,aeout,sdhist
      common /conversion/ radtodeg
      common /unit5/ iunit5
      common /moov/ movie,nframes,icall1,lhdr
      common /elastic_ss/ idef_ss
      common /memory/ memadd,memaddi
c
      icall  = 0
      memuse = 0
c
      call umalloc(ip_ltot,maxbl0,1,icall,memuse)
      call umalloc(ip_jjmax1,nsub10,1,icall,memuse)
      call umalloc(ip_kkmax1,nsub10,1,icall,memuse)
      call umalloc(ip_iiint1,nsub10,1,icall,memuse)
      call umalloc(ip_iiint2,nsub10,1,icall,memuse)
      call umalloc(ip_n14,maxbl0,1,icall,memuse)
      call umalloc(ip_iwrk,maxbl0*3,1,icall,memuse)
      call umalloc(ip_icsinfo,maxcs0*9,1,icall,memuse)
      call umalloc(ip_iv,maxbl0,1,icall,memuse)
      call umalloc(ip_memblock,maxbl0,1,icall,memuse)
      call umalloc(ip_lwdat,maxbl0*maxseg0*6,1,icall,memuse)
      call umalloc(ip_nblfine,maxbl0,1,icall,memuse)
      call umalloc(ip_mem_req_node,maxnode0,1,icall,memuse)
      call umalloc(ip_no_of_points,maxnode0,1,icall,memuse)
      call umalloc(ip_lw,65*maxbl0,1,icall,memuse)
      call umalloc(ip_lw2,43*maxbl0,1,icall,memuse)
      call umalloc(ip_nblk,2*mxbli0,1,icall,memuse)
      call umalloc(ip_limblk,12*mxbli0,1,icall,memuse)
      call umalloc(ip_isva,4*mxbli0,1,icall,memuse)
      call umalloc(ip_nblon,mxbli0,1,icall,memuse)
      call umalloc(ip_lig,maxbl0,1,icall,memuse)
      call umalloc(ip_lbg,maxbl0,1,icall,memuse)
      call umalloc(ip_iovrlp,maxbl0,1,icall,memuse)
      call umalloc(ip_ibpntsg,4*maxbl0,1,icall,memuse)
      call umalloc(ip_iipntsg,maxbl0,1,icall,memuse)
      call umalloc(ip_rkap0g,3*maxbl0,0,icall,memuse)
      call umalloc(ip_levelg,maxbl0,1,icall,memuse)
      call umalloc(ip_igridg,maxbl0,1,icall,memuse)
      call umalloc(ip_iflimg,3*maxbl0,1,icall,memuse)
      call umalloc(ip_ifdsg,3*maxbl0,1,icall,memuse)
      call umalloc(ip_iviscg,3*maxbl0,1,icall,memuse)
      call umalloc(ip_jdimg,maxbl0,1,icall,memuse)
      call umalloc(ip_kdimg,maxbl0,1,icall,memuse)
      call umalloc(ip_idimg,maxbl0,1,icall,memuse)
      call umalloc(ip_idiagg,3*maxbl0,1,icall,memuse)
      call umalloc(ip_nblcg,maxbl0,1,icall,memuse)
      call umalloc(ip_idegg,3*maxbl0,1,icall,memuse)
      call umalloc(ip_jsg,maxbl0,1,icall,memuse)
      call umalloc(ip_ksg,maxbl0,1,icall,memuse)
      call umalloc(ip_isg,maxbl0,1,icall,memuse)
      call umalloc(ip_jeg,maxbl0,1,icall,memuse)
      call umalloc(ip_keg,maxbl0,1,icall,memuse)
      call umalloc(ip_ieg,maxbl0,1,icall,memuse)
      call umalloc(ip_mit,5*maxbl0,1,icall,memuse)
      call umalloc(ip_jlamlog,maxbl0,1,icall,memuse)
      call umalloc(ip_klamlog,maxbl0,1,icall,memuse)
      call umalloc(ip_ilamlog,maxbl0,1,icall,memuse)
      call umalloc(ip_jlamhig,maxbl0,1,icall,memuse)
      call umalloc(ip_klamhig,maxbl0,1,icall,memuse)
      call umalloc(ip_ilamhig,maxbl0,1,icall,memuse)
      call umalloc(ip_iwfg,3*maxbl0,1,icall,memuse)
      call umalloc(ip_utrans,maxbl0,0,icall,memuse)
      call umalloc(ip_vtrans,maxbl0,0,icall,memuse)
      call umalloc(ip_wtrans,maxbl0,0,icall,memuse)
      call umalloc(ip_omegax,maxbl0,0,icall,memuse)
      call umalloc(ip_omegay,maxbl0,0,icall,memuse)
      call umalloc(ip_omegaz,maxbl0,0,icall,memuse)
      call umalloc(ip_xorig,maxbl0,0,icall,memuse)
      call umalloc(ip_yorig,maxbl0,0,icall,memuse)
      call umalloc(ip_zorig,maxbl0,0,icall,memuse)
      call umalloc(ip_dxmx,maxbl0,0,icall,memuse)
      call umalloc(ip_dymx,maxbl0,0,icall,memuse)
      call umalloc(ip_dzmx,maxbl0,0,icall,memuse)
      call umalloc(ip_dthxmx,maxbl0,0,icall,memuse)
      call umalloc(ip_dthymx,maxbl0,0,icall,memuse)
      call umalloc(ip_dthzmx,maxbl0,0,icall,memuse)
      call umalloc(ip_thetax,maxbl0,0,icall,memuse)
      call umalloc(ip_thetay,maxbl0,0,icall,memuse)
      call umalloc(ip_thetaz,maxbl0,0,icall,memuse)
      call umalloc(ip_rfreqt,maxbl0,0,icall,memuse)
      call umalloc(ip_rfreqr,maxbl0,0,icall,memuse)
      call umalloc(ip_xorig0,maxbl0,0,icall,memuse)
      call umalloc(ip_yorig0,maxbl0,0,icall,memuse)
      call umalloc(ip_zorig0,maxbl0,0,icall,memuse)
      call umalloc(ip_time2,maxbl0,0,icall,memuse)
      call umalloc(ip_thetaxl,maxbl0,0,icall,memuse)
      call umalloc(ip_thetayl,maxbl0,0,icall,memuse)
      call umalloc(ip_thetazl,maxbl0,0,icall,memuse)
      call umalloc(ip_itrans,maxbl0,1,icall,memuse)
      call umalloc(ip_irotat,maxbl0,1,icall,memuse)
      call umalloc(ip_idefrm,maxbl0,1,icall,memuse) 
      call umalloc(ip_bcvali,14*maxseg0*maxbl0,0,icall,memuse)
      call umalloc(ip_bcvalj,14*maxseg0*maxbl0,0,icall,memuse)
      call umalloc(ip_bcvalk,14*maxseg0*maxbl0,0,icall,memuse)
      call umalloc(ip_nbci0,maxbl0,1,icall,memuse)
      call umalloc(ip_nbcj0,maxbl0,1,icall,memuse)
      call umalloc(ip_nbck0,maxbl0,1,icall,memuse)
      call umalloc(ip_nbcidim,maxbl0,1,icall,memuse)
      call umalloc(ip_nbcjdim,maxbl0,1,icall,memuse)
      call umalloc(ip_nbckdim,maxbl0,1,icall,memuse)
      call umalloc(ip_ibcinfo,14*maxseg0*maxbl0,1,icall,memuse)
      call umalloc(ip_jbcinfo,14*maxseg0*maxbl0,1,icall,memuse)
      call umalloc(ip_kbcinfo,14*maxseg0*maxbl0,1,icall,memuse)
      call umalloc(ip_bcfilei,2*maxseg0*maxbl0,1,icall,memuse)
      call umalloc(ip_bcfilej,2*maxseg0*maxbl0,1,icall,memuse)
      call umalloc(ip_bcfilek,2*maxseg0*maxbl0,1,icall,memuse)
      call umalloc(ip_ncgg,maxgr0,1,icall,memuse)
      call umalloc(ip_nblg,maxgr0,1,icall,memuse)
      call umalloc(ip_iemg,maxgr0,1,icall,memuse)
      call umalloc(ip_inewgg,maxgr0,1,icall,memuse)
      call umalloc(ip_inpl3d,11*nplots0,1,icall,memuse)
      call umalloc(ip_inpr,11*nplots0,1,icall,memuse)
      call umalloc(ip_iadvance,maxbl0,1,icall,memuse)
      call umalloc(ip_iforce,maxbl0,1,icall,memuse)
      call umalloc(ip_iindex,intmax0*(6*nsub10+9),1,icall,memuse)
      call umalloc(ip_iindx,intmx0*(6*msub10+9),1,icall,memuse)
      call umalloc(ip_llimit,intmx0,1,icall,memuse)
      call umalloc(ip_iitmax,intmx0,1,icall,memuse)
      call umalloc(ip_mmcxie,intmx0,1,icall,memuse)
      call umalloc(ip_mmceta,intmx0,1,icall,memuse)
      call umalloc(ip_ncheck,maxbl0,1,icall,memuse)
      call umalloc(ip_iifit,intmx0,1,icall,memuse)
      call umalloc(ip_mblkpt,mxxe0,1,icall,memuse)
      call umalloc(ip_iic0,intmx0,1,icall,memuse)
      call umalloc(ip_iiorph,intmx0,1,icall,memuse)
      call umalloc(ip_iitoss,intmx0,1,icall,memuse)
      call umalloc(ip_ifiner,intmx0,1,icall,memuse)
      call umalloc(ip_dx,intmx0*msub10,0,icall,memuse)
      call umalloc(ip_dy,intmx0*msub10,0,icall,memuse)
      call umalloc(ip_dz,intmx0*msub10,0,icall,memuse)
      call umalloc(ip_dthetx,intmx0*msub10,0,icall,memuse)
      call umalloc(ip_dthety,intmx0*msub10,0,icall,memuse)
      call umalloc(ip_dthetz,intmx0*msub10,0,icall,memuse)
      call umalloc(ip_isav_blk,34*mxbli0,1,icall,memuse)
      call umalloc(ip_isav_prd,12*lbcprd0,1,icall,memuse)
      call umalloc(ip_isav_pat,17*intmax0,1,icall,memuse)
      call umalloc(ip_isav_pat_b,6*intmax0*nsub10,1,icall,memuse)
      call umalloc(ip_isav_dpat,17*intmx0,1,icall,memuse)
      call umalloc(ip_isav_dpat_b,6*intmx0*msub10,1,icall,memuse)
      call umalloc(ip_isav_emb,12*lbcemb0,1,icall,memuse)
      call umalloc(ip_mblk2nd,maxbl0,1,icall,memuse)
      call umalloc(ip_mglevg,maxbl0,1,icall,memuse)
      call umalloc(ip_nemgl,maxbl0,1,icall,memuse)
      call umalloc(ip_ipl3dtmp,11*nplots0,1,icall,memuse)
      call umalloc(ip_ifrom,msub10,1,icall,memuse)
      call umalloc(ip_xif1,msub10,1,icall,memuse)
      call umalloc(ip_etf1,msub10,1,icall,memuse)
      call umalloc(ip_xif2,msub10,1,icall,memuse)
      call umalloc(ip_etf2,msub10,1,icall,memuse)
      call umalloc(ip_utrnsae,maxbl0*maxsegdg0,0,icall,memuse)
      call umalloc(ip_vtrnsae,maxbl0*maxsegdg0,0,icall,memuse)
      call umalloc(ip_wtrnsae,maxbl0*maxsegdg0,0,icall,memuse)
      call umalloc(ip_omgxae,maxbl0*maxsegdg0,0,icall,memuse)
      call umalloc(ip_omgyae,maxbl0*maxsegdg0,0,icall,memuse)
      call umalloc(ip_omgzae,maxbl0*maxsegdg0,0,icall,memuse)
      call umalloc(ip_xorgae,maxbl0*maxsegdg0,0,icall,memuse)
      call umalloc(ip_yorgae,maxbl0*maxsegdg0,0,icall,memuse)
      call umalloc(ip_zorgae,maxbl0*maxsegdg0,0,icall,memuse)
      call umalloc(ip_xorgae0,maxbl0*maxsegdg0,0,icall,memuse)
      call umalloc(ip_yorgae0,maxbl0*maxsegdg0,0,icall,memuse)
      call umalloc(ip_zorgae0,maxbl0*maxsegdg0,0,icall,memuse)
      call umalloc(ip_icouple,maxbl0*maxsegdg0,1,icall,memuse)
      call umalloc(ip_thtxae,maxbl0*maxsegdg0,0,icall,memuse)
      call umalloc(ip_thtyae,maxbl0*maxsegdg0,0,icall,memuse)
      call umalloc(ip_thtzae,maxbl0*maxsegdg0,0,icall,memuse)
      call umalloc(ip_rfrqtae,maxbl0*maxsegdg0,0,icall,memuse)
      call umalloc(ip_rfrqrae,maxbl0*maxsegdg0,0,icall,memuse)
      call umalloc(ip_icsi,maxbl0*maxsegdg0,1,icall,memuse)
      call umalloc(ip_icsf,maxbl0*maxsegdg0,1,icall,memuse)
      call umalloc(ip_jcsi,maxbl0*maxsegdg0,1,icall,memuse)
      call umalloc(ip_jcsf,maxbl0*maxsegdg0,1,icall,memuse)
      call umalloc(ip_kcsi,maxbl0*maxsegdg0,1,icall,memuse)
      call umalloc(ip_kcsf,maxbl0*maxsegdg0,1,icall,memuse)
      call umalloc(ip_idfrmseg,maxbl0*maxsegdg0,1,icall,memuse)
      call umalloc(ip_iaesurf,maxbl0*maxsegdg0,1,icall,memuse)
      call umalloc(ip_nsegdfrm,maxbl0,1,icall,memuse)
      call umalloc(ip_freq,nmds0*maxaes0,0,icall,memuse)
      call umalloc(ip_gmass,nmds0*maxaes0,0,icall,memuse)
      call umalloc(ip_x0,2*nmds0*maxaes0,0,icall,memuse)
      call umalloc(ip_gf0,2*nmds0*maxaes0,0,icall,memuse)
      call umalloc(ip_damp,nmds0*maxaes0,0,icall,memuse)
      call umalloc(ip_perturb,4*nmds0*maxaes0,0,icall,memuse)
      call umalloc(ip_aesrfdat,5*maxaes0,0,icall,memuse)
      call umalloc(ip_iskip,maxbl0,1,icall,memuse)
      call umalloc(ip_jskip,maxbl0,1,icall,memuse)
      call umalloc(ip_kskip,maxbl0,1,icall,memuse)
c
c     ierrflg = -99 to indicate errors occur during memory
c     allocation phase
c
      ierrflg = -99
c
c     default to dimensions of 1
c
      lmaxgr   = 1
      lmaxbl   = 1
      lmxseg   = 1
      lmaxcs   = 1
      lnplts   = 1
      lmxbli   = 1
      lmaxxe   = 1
      lnsub1   = 1
      lintmx   = 1
      lmxxe    = 1
      liitot   = 1
      isum     = 1
      lncycm   = 1
      isum_h   = 1
      isum_n   = 1
      lminnode = 1
      isumi    = 1
      isumi_n  = 1
      lmptch   = 1
      lmsub1   = 1
      lintmax  = 1
      libufdim = 1
      lnbuf    = 1
      llbcprd  = 1
      llbcemb  = 1
      llbcrad  = 1
      lmaxaes  = 1
      lnmds    = 1
      lnslave  = 1
      lnmaster = 1
      lmxsegdg = 1
c
      read(iunit5,*)
      read(iunit5,'(a60)')grid
      read(iunit5,'(a60)')plt3dg
      read(iunit5,'(a60)')plt3dq
      read(iunit5,'(a60)')output
      read(iunit5,'(a60)')residual
      read(iunit5,'(a60)')turbres
      read(iunit5,'(a60)')blomx
      read(iunit5,'(a60)')output2
      read(iunit5,'(a60)')printout
      read(iunit5,'(a60)')pplunge
c     the ovrlap file is opened (if needed) as unit 21 in global
      read(iunit5,'(a60)')ovrlap
c     the patch file is opened (if needed) as unit 22 in global
      read(iunit5,'(a60)')patch
      read(iunit5,'(a60)')restrt
c
c***********************************************************************
c 
      gamma    = 1.4e0
      gm1      = gamma-1.0e0
      gp1      = gamma+1.0e0
      gm1g     = gm1/gamma
      gp1g     = gp1/gamma
      ggm1     = gamma*gm1
      radtodeg = 180.e0/3.1415927e0
c
c     output banner
c
      write(66,83)
      write(66,83)
      write(66,87)
      write(66,9900)
 9900 format(2(2h *),47h                  PRECFL3D - CFL3D PREPROCESSOR,
     .11h           ,4x,2(2h *))
      write(66,87)
      write(66,9990)
 9990 format(2(2h *),43h   VERSION 6.0 :  Computational Fluids Lab,,
     .15h Mail Stop 128,,4x,2(2h *),
     ./2(2h *),18x,41hNASA Langley Research Center, Hampton, VA,
     .3x,2(2h *),/2(2h *),18x,33hRelease Date:    August 14, 2001.,
     .11x,2(2h *))
      write(66,87)
      write(66,83)
      write(66,83)
   83 format(35(2h *))
   87 format(2(2h *),62x,2(2h *))
c
#ifdef CRAY_TIME
c     cray_time implies cray (always double precision)
      write(66,12) real(float(memuse))/1.e6
#else
#   ifdef DBLE_PRECSN
      write(66,12) real(float(memuse))/1.e6
#   else
      write(66,13) real(float(memuse))/1.e6
#   endif
#endif
   12 format(/,' memory allocation: ',f12.6,' Mbytes, double precision')
   13 format(/,' memory allocation: ',f12.6,' Mbytes, single precision')
c
      write(66,88)
   88 format(/19hinput/output files:)
c
      write(66,'(''  '',a60)')grid
      write(66,'(''  '',a60)')plt3dg
      write(66,'(''  '',a60)')plt3dq
      write(66,'(''  '',a60)')output
      write(66,'(''  '',a60)')residual
      write(66,'(''  '',a60)')turbres
      write(66,'(''  '',a60)')blomx
      write(66,'(''  '',a60)')output2
      write(66,'(''  '',a60)')printout
      write(66,'(''  '',a60)')pplunge
      write(66,'(''  '',a60)')ovrlap
      write(66,'(''  '',a60)')patch
      write(66,'(''  '',a60)')restrt
c
c***********************************************************************
c     read in global information
c***********************************************************************
c
c     output precfl3d info to unit 66
c
      iunit = 66
c
      do ll=1,nbuf0
         nou(ll) = 0
      end do
      icallgl = 0
      call global(myid,maxbl0,maxgr0,maxseg0,maxcs0,nplots0,mxbli0,
     .            bcvali,bcvalj,bcvalk,nbci0,nbcj0,nbck0,
     .            nbcidim,nbcjdim,nbckdim,ibcinfo,jbcinfo,
     .            kbcinfo,bcfilei,bcfilej,bcfilek,nblk,nbli,
     .            limblk,isva,nblon,rkap0g,nblock,levelg,
     .            igridg,iflimg,ifdsg,iviscg,jdimg,kdimg,
     .            idimg,idiagg,nblcg,idegg,jsg,ksg,isg,jeg,
     .            keg,ieg,mit,ilamlog,ilamhig,jlamlog,
     .            jlamhig,klamlog,klamhig,iwfg,utrans,vtrans,
     .            wtrans,omegax,omegay,omegaz,xorig,yorig,
     .            zorig,dxmx,dymx,dzmx,dthxmx,dthymx,
     .            dthzmx,thetax,thetay,thetaz,rfreqt,rfreqr,
     .            xorig0,yorig0,zorig0,itrans,irotat,idefrm,ngrid,
     .            ncgg,nblg,iemg,inewgg,iovrlp,ninter,
     .            nplot3d,inpl3d,ip3dsurf,nprint,inpr,
     .            iadvance,iforce,lfgm,ncs,icsinfo,ihstry,
     .            ncycmax0,iv,time2,thetaxl,thetayl,thetazl,
     .            intmax0,nsub10,iindex,lig,lbg,ibpntsg,
     .            iipntsg,icallgl,iunit,nou,bou,ibufdim0,nbuf0,
     .            mglevg,nemgl,ipl3dtmp,ntr,bcfiles,mxbcfil0,
     .            utrnsae,vtrnsae,wtrnsae,omgxae,
     .            omgyae,omgzae,xorgae,yorgae,zorgae,thtxae,
     .            thtyae,thtzae,rfrqtae,rfrqrae,icsi,icsf,jcsi,jcsf,
     .            kcsi,kcsf,freq,gmass,damp,x0,gf0,nmds0,maxaes0,
     .            aesrfdat,perturb,iskip,jskip,kskip,nsegdfrm,
     .            idfrmseg,iaesurf,maxsegdg0,xorgae0,yorgae0,zorgae0,
     .            icouple,iprnsurf)
c
c     read dynamic patch input data
c
      if (iunst.gt.0) then
            ioflag = 2
            imode2 = 1
            call global2(maxbl0,maxgr0,msub10,nintr,intmx0,ngrid,idimg,
     .                   jdimg,kdimg,levelg,ncgg,nblg,iindx,llimit,
     .                   iitmax,mmcxie,mmceta,ncheck,iifit,iic0,
     .                   iiorph,iitoss,ifiner,dx,dy,dz,dthetx,
     .                   dthety,dthetz,myid,mptch0,mxxe0,icallgl,iunit,
     .                   nou,bou,ibufdim0,nbuf0,ifrom,xif1,etf1,xif2,
     .                   etf2,igridg,iemg,nblock,ioflag,imode2)
      end if
c
c***********************************************************************
c
c     allocate blocks to nodes
c
c***********************************************************************
c
      call compg2n(nblock,ngrid,ncgg,nblg,idimg,jdimg,kdimg,
     .             nblcg,nnodes,iwrk,myid,myhost,mblk2nd,
     .             mycomm,maxgr0,maxbl0,ierrflg,ibufdim0,
     .             nbuf0,bou,nou)
c
      write(iunit,200)
      write(iunit,201) nblock,nnodes
      do i=1,nblock
         write(iunit,202) i,mblk2nd(i)
      end do
  200 format(/,1x,'BLOCK TO NODE MAPPING')
  201 format(1x,'no. of blocks = ',i4,/,1x,'no. of  nodes = ',i4,/,
     .1x,'block    node')
  202 format(i5,3x,i5)
c
c     check total number of nodes that are actually used
c
      nodel = 1
      do nbl=1,nblock
         if (mblk2nd(nbl) .gt. nodel) nodel = mblk2nd(nbl)
      end do
c
      lminnode = nodel
c
c***********************************************************************
c
c     evaluate requirements for permanent array size (w array)
c
c***********************************************************************
c
c     set up lw, lw2 pointer arrays
c
      icallpt = 0
c
      do nbl=1,nblock
         memblock(nbl) = 0
      end do
c
c     call pointers with myiduse=myhost and mpihost=0 to determine nstart
c     for sequential host
c
      myiduse = myhost
      mpihost = 0
      call pointers(lw,lw2,maxl,lembed,nstart_hseq,nwork,mwork,maxbl0,
     .              maxgr0,maxseg0,lwdat,levelg,igridg,iviscg,
     .              idimg,jdimg,kdimg,nblcg,itrans,irotat,idefrm,
     .              nbci0,nbcj0,nbck0,nbcidim,nbcjdim,nbckdim,
     .              ibcinfo,jbcinfo,kbcinfo,ngrid,ncgg,nblg,
     .              iemg,nblock,myhost,myiduse,mblk2nd,nou,bou,nbuf0,
     .              ibufdim0,nblfine,ilamlog,jlamlog,
     .              klamlog,ilamhig,jlamhig,klamhig,idegg,iwfg,
     .              idiagg,iflimg,ifdsg,rkap0g,jsg,ksg,isg,jeg,
     .              keg,ieg,memblock,icallpt,nmds0,maxaes0,mpihost)
c
c     output contents of lw() for sequential host
c
      write(66,*)'  '
      write(66,*)'***********************************************'
      write(66,*)
      write(66,*)'    PERMANENT STORAGE REQUIREMENTS - W ARRAY'
      write(66,*)'                SEQUENTIAL BUILD'
      write(66,*)
      write(66,*)'***********************************************'
      write(66,910)
c
      do nbl=1,nblock
         igrid = igridg(nbl)
         j = jdimg(nbl)
         k = kdimg(nbl)
         i = idimg(nbl)
         write(66,223)
         write(66,4)igrid,nbl,j,k,i
         do ii=1,63
            write(66,37)ii,nbl,lw(ii,nbl)
         end do
         do jj=1,6
            if(jj.eq.1) nseg = nbci0(nbl)
            if(jj.eq.2) nseg = nbcidim(nbl)
            if(jj.eq.3) nseg = nbcj0(nbl)
            if(jj.eq.4) nseg = nbcjdim(nbl)
            if(jj.eq.5) nseg = nbck0(nbl)
            if(jj.eq.6) nseg = nbckdim(nbl)
            do ii=1,nseg
               write(66,38)nbl,ii,jj,lwdat(nbl,ii,jj)
            end do
         end do
      end do
c
c     call pointers with myiduse=myhost and mpihost=1 to determine nstart
c     for parallel host
c
      myiduse = myhost
      mpihost = 1
      call pointers(lw,lw2,maxl,lembed,nstart_hmpi,nwork,mwork,maxbl0,
     .              maxgr0,maxseg0,lwdat,levelg,igridg,iviscg,
     .              idimg,jdimg,kdimg,nblcg,itrans,irotat,idefrm,
     .              nbci0,nbcj0,nbck0,nbcidim,nbcjdim,nbckdim,
     .              ibcinfo,jbcinfo,kbcinfo,ngrid,ncgg,nblg,
     .              iemg,nblock,myhost,myiduse,mblk2nd,nou,bou,nbuf0,
     .              ibufdim0,nblfine,ilamlog,jlamlog,
     .              klamlog,ilamhig,jlamhig,klamhig,idegg,iwfg,
     .              idiagg,iflimg,ifdsg,rkap0g,jsg,ksg,isg,jeg,
     .              keg,ieg,memblock,icallpt,nmds0,maxaes0,mpihost)
c
c     loop over parallel nodes
c
      do myiduse = 1,nnodes
         call pointers(lw,lw2,maxl,lembed,nstart,nwork,mwork,maxbl0,
     .                 maxgr0,maxseg0,lwdat,levelg,igridg,iviscg,
     .                 idimg,jdimg,kdimg,nblcg,itrans,irotat,idefrm,
     .                 nbci0,nbcj0,nbck0,nbcidim,nbcjdim,nbckdim,
     .                 ibcinfo,jbcinfo,kbcinfo,ngrid,ncgg,nblg,
     .                 iemg,nblock,myhost,myiduse,mblk2nd,nou,bou,nbuf0,
     .                 ibufdim0,nblfine,ilamlog,jlamlog,
     .                 klamlog,ilamhig,jlamhig,klamhig,idegg,iwfg,
     .                 idiagg,iflimg,ifdsg,rkap0g,jsg,ksg,isg,jeg,
     .                 keg,ieg,memblock,icallpt,nmds0,maxaes0,mpihost)
      end do
c
c     determine permanent w array requirement for nodes
c
      do ii = 1,nodel
c
         mem_req_node(ii) = 0
         no_of_points(ii) = 0
c
         do nbl = 1,nblock
c
c           determine if there are any finer embedded blocks for
c           global block nbl
c
            nfiner = 0
            do mbl = 1,nblock
               if (levelg(mbl) .gt. lglobal) then
                  if (nblcg(mbl).eq.nbl) then
                     nfiner = nfiner + 1
                     nblfine(nfiner) = mbl
                  end if
               end if
            end do
c
            nblc = nblcg(nbl)
            iaug = 0
            if (ii.eq.mblk2nd(nbl)) iaug = 1
            if (ii.eq.mblk2nd(nblc)) iaug = 1
            if (nfiner.gt.0) then
               do nf = 1,nfiner
                  if (ii.eq.mblk2nd(nblfine(nf))) iaug = 1
               end do
            end if
c
            if (iaug.gt.0) then
               npoints = (idimg(nbl)-1)*(jdimg(nbl)-1)*(kdimg(nbl)-1)
               mem_req_node(ii) = mem_req_node(ii) + memblock(nbl)
               no_of_points(ii) = no_of_points(ii) + npoints
            end if
c
         end do
      end do
c
      npts_max  = 0
      mem_w_max = nstart_hmpi
      do i = 1,nodel
         if (mem_req_node(i).gt.mem_w_max) mem_w_max = mem_req_node(i)
         if (no_of_points(i).gt.npts_max ) npts_max  = no_of_points(i)
      end do
c
      write(66,*)
      write(66,*)'***********************************************'
      write(66,*)
      write(66,*)'   PERMANENT STORAGE REQUIREMENTS - W ARRAY    '
      write(66,*)'            SUMMARY FOR ALL BUILDS'
      write(66,*)
      write(66,*)'***********************************************'
c
      ftot_seq = float(nstart_hseq)/1.e+06
      write(66,2001) real(ftot_seq)
c
      write(66,8425)
      write(66,8426)
      fmemsum = nstart_hmpi/1.e+06
      nptsum  = 0
      write(66,8427) myhost, npts_max, real(fmemsum)
      do i = 1,nodel
         write(66,8427) i, no_of_points(i), mem_req_node(i)/1.e+06
         fmemsum = fmemsum + mem_req_node(i)/1.e+06
         nptsum  = nptsum + no_of_points(i)
         mem_req_node(i) = 0
      end do
      write(66,8428) nptsum,real(fmemsum)
      write(66,8429) npts_max,real(float(mem_w_max))/1.e+06
c
      write(66,*)
c
      ftot_seq = float(nstart_hseq)
c
    4 format(5i6)
   37 format(2x,4h lw(,i2,1h,,i3,3h)= ,6x,i10)
   38 format(2x,7h lwdat(,i3,1h,,i2,1h,,i2,3h)= ,i10)
  223 format(/1x,5higrid,1x,5hblock,2x,4hjdim,2x,4hkdim,2x,4hidim)
  910 format(/52h summary of starting locations for block information,
     .15h on single node)
 2001 format(/,1x,41hmemory (mw) for w storage (sequential) = ,f9.4)
 3011 format(/1x,36havailable memory (mw) for wk storage,
     .23h in this preprocessor =,f9.5,/)
 8425 format(/,1x,33hmemory (mw) for w storage (nodes))
 8426 format(/,1x,33h node   total points  memory (mw))
 8427 format(i6,7x,i8,4x,f9.4)
 8428 format(/,5x,8htotal = ,i8,4x,f9.4)
 8429 format(/,3x,10hmaximum = ,i8,4x,f9.4)
c
c***********************************************************************
c
c     evaluate requirements for temporary real array size 
c     (work/wk array)
c
c***********************************************************************
c
      write(66,*)
      write(66,*)'***********************************************'
      write(66,*)
      write(66,*)'   TEMPORARY STORAGE REQUIREMENTS - WK ARRAY'
      write(66,*)'        SEQUENTIAL AND PARALLEL BUILDS'
      write(66,*)
      write(66,*)'***********************************************'
c
      do n=1,maxbl0
         ltot(n) = 0
      end do
c
      need  = lmaxbl
      needi = 1
      needi_node = 1
c
c     requirements of subroutine addx
c
      do 8000 iseq=1,mseq
      if (ncyc1(iseq).le.0) go to 8500
      if (iseq.gt.1 .or. iseq.eq.mseq) then
         if (iseq.le.mseq .and. iseq.ne.1) then
            do 7040 igrid=1,ngrid
            iem   = iemg(igrid)
            if (iem.gt.0) go to 7040
            nblf  = nblg(igrid)+(mseq-iseq)
            nblz  = nblf+1
            call lead(nblf,lw,lw2,maxbl0)
            iwk1  = 1
            iwk2  = iwk1+jdim*kdim*idim*5
            iwk3  = iwk2+jj2*kk2*ii2*5
            itemp = iwk3 + jdim*kk2*ii2*5
            need  = max(need,itemp)
            ii = mblk2nd(nblf)
            mem_req_node(ii) = max(mem_req_node(ii),itemp)
            write(66,*)' 1(addx)         itemp, need = ',itemp,need
            write(66,*)'       parallel: itemp, need = ',itemp,
     .      mem_req_node(ii)
 7040       continue
         end if
c
c     requirements of subroutine add2x
c
         if (lembed.eq.0 .or. iseq.ne.mseq) go to 333
         do 8100 level=mseq+1,maxl
         do 8040 igrid=1,ngrid
         nbl   = nblg(igrid)
         iem   = iemg(igrid)
         inewg = inewgg(igrid)
         if (iem.eq.0) go to 8040
         if (iseq.eq.1 .and. inewg.eq.0) go to 8040
         if (levelg(nbl).ne.level) go to 8040
         call lead(nbl,lw,lw2,maxbl0)
         jj2   = jdimg(nblc)
         kk2   = kdimg(nblc)
         ii2   = idimg(nblc)
         iwk1  = 1
         iwk2  = iwk1 + jdim*kdim*idim*5
         iwk3  = iwk2 + jj2*kk2*ii2*5
         iwk4  = iwk3 + jdim*kk2*ii2
         itemp = iwk4 + jdim*kdim*ii2
         need  = max(need,itemp)
         ii = mblk2nd(nbl)
         mem_req_node(ii) = max(mem_req_node(ii),itemp)
         write(66,*)' 2(add2x)        itemp, need = ',itemp,need
         write(66,*)'       parallel: itemp, need = ',itemp,
     .   mem_req_node(ii)
 8040    continue
 8100    continue
      end if
c
  333 levt  = levelt(iseq)
      levb  = levelb(iseq)
      level = levelt(iseq)
      lglobal = lfgm-(mseq-iseq)
 7000 continue
      do 6 i=1,level-levb+1
      nsm(i) = ngam
    6 continue
      nit    = mit(level-levb+1,iseq)
      kxpand = 1
 9000 continue
      ilc    = 0
c
      do 2000 nbl=1,nblock
      if (level.ne.levelg(nbl)) go to 2000
      call lead(nbl,lw,lw2,maxbl0)
c
      nit1 = nit + min(1,mgflag)
      if (level.eq.levb) nit1 = nit
      if (level.gt.lglobal) nit1 = nit+1
      if (level.eq.lglobal) nit1 = nit
      if (nit1.gt.0) then
         do 6501 ntime=1,nit1
         if (ilc.eq.0) then
            nsm(level-levb+1) = nsm(level-levb+1)-1
            if (nsm(level-levb+1).lt.0) nsm(level-levb+1) = ngam-1
            ilc = 1
         end if
         igrid = igridg(nbl)
         iem   = iemg(igrid)
         ifluxa = 0
         if (level.ge.lglobal .and.
     .       level.ne.levt) ifluxa = min(1,iconsf)
         if (ifluxa.gt.0) then
            if (ntime.eq.1) then
            call cntfa(nbl,nwfa,maxbl0,maxseg0,nblcg,ieg,isg,jdimg,
     .                 kdimg,idimg,nbcj0,nbck0,nbci0,nbcjdim,nbckdim,
     .                 nbcidim,jbcinfo,kbcinfo,ibcinfo,nblock)
            end if
         else
            nwfa = 0
         end if
c
         lres      = nwfa  + 1
         lwj0      = lres  + jdim*kdim*(idim-1)*5
         lwk0      = lwj0  + kdim*idim*22
         lwi0      = lwk0  + jdim*idim*22
         lvmuk     = lwi0  + kdim*jdim*22
         lvmuj     = lvmuk + 2*(jdim-1)*(idim-1)
         lvmui     = lvmuj + 2*(kdim-1)*(idim-1)
         ltot(nbl) = lvmui + 2*(kdim-1)*(jdim-1) + maxbl0
         itemp=0
c
         if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2)then
            if (ivisc(3).eq.2 .or. ivisc(2).eq.2 .or. ivisc(1).eq.2)then
c              Baldwin-Lomax
               inmx = max(jdim,kdim,idim)
               iwk1 = 1
               iwk2 = iwk1+inmx
               iwk3 = iwk2+inmx
               iwk4 = iwk3+inmx
               iwk5 = iwk4+inmx
               iwk6 = iwk5+inmx
               iwk7 = iwk6+inmx
               iwk8 = iwk7+inmx
               iwk9 = iwk8+inmx
               itemp= iwk9+inmx
            else if (ivisc(3).eq.4 .or. ivisc(2).eq.4 .or. 
     .               ivisc(1).eq.4) then
c              Baldwin-Barth
               iwk1=1
               iwk5=iwk1+(jdim+1)*(kdim+1)*(idim+1)
               inmx=(kdim-1)*(jdim-1)
               iwk6=iwk5+inmx
               iwk7=iwk6+inmx
               iwk8=iwk7+inmx
               iwk9=iwk8+inmx
               iwk10=iwk9+inmx
               inmx=(jdim-1)*(kdim-1)
               iwk11=iwk10+inmx
               iwk12=iwk11+inmx
               iwk13=iwk12+inmx
               iwk14=iwk13+inmx
               iwk15=iwk14+inmx
               inmx=(kdim-1)*(idim-1)
               iwk16=iwk15+inmx
               iwk17=iwk16+inmx
               iwk18=iwk17+inmx
               iwk19=iwk18+inmx
               itemp=iwk19+inmx
            else if (ivisc(3).eq.5 .or. ivisc(2).eq.5 .or. 
     .               ivisc(1).eq.5) then
c              Spalart
               iwk1=1
               iwk3=iwk1+(jdim+1)*(kdim+1)*(idim+1)
               iwk4=iwk3+(jdim+1)*(kdim+1)*(idim+1)
               inmx=(kdim-1)*(jdim-1)
               iwk5=iwk4+inmx
               iwk6=iwk5+inmx
               iwk7=iwk6+inmx
               iwk8=iwk7+inmx
               iwk9=iwk8+inmx
               inmx=(jdim-1)*(kdim-1)
               iwk10=iwk9+inmx
               iwk11=iwk10+inmx
               iwk12=iwk11+inmx
               iwk13=iwk12+inmx
               iwk14=iwk13+inmx
               inmx=(kdim-1)*(idim-1)
               iwk15=iwk14+inmx
               iwk16=iwk15+inmx
               iwk17=iwk16+inmx
               iwk18=iwk17+inmx
               itemp=iwk18+inmx
            else if (ivisc(3).ge.6 .or. ivisc(2).ge.6 .or.
     .               ivisc(1).ge.6) then
c              Two-equation turb models
               iwk1=1
               iwk4=iwk1+(jdim+1)*(kdim+1)*(idim+1)*2
               iwk5=iwk4+(jdim+1)*(kdim+1)*(idim+1)
               inmx=(kdim-1)*(jdim-1)
               iwk6=iwk5+inmx
               iwk7=iwk6+inmx
               iwk8=iwk7+inmx
               iwk9=iwk8+inmx
               iwk10=iwk9+inmx
               iwk11=iwk10+inmx
               iwk12=iwk11+inmx
               iwk13=iwk12+inmx
               iwk14=iwk13+inmx
               inmx=(jdim-1)*(kdim-1)
               iwk15=iwk14+inmx
               iwk16=iwk15+inmx
               iwk17=iwk16+inmx
               iwk18=iwk17+inmx
               iwk19=iwk18+inmx
               iwk20=iwk19+inmx
               iwk21=iwk20+inmx
               iwk22=iwk21+inmx
               iwk23=iwk22+inmx
               inmx=(kdim-1)*(idim-1)
               iwk24=iwk23+inmx
               iwk25=iwk24+inmx
               iwk26=iwk25+inmx
               iwk27=iwk26+inmx
               iwk28=iwk27+inmx
               iwk29=iwk28+inmx
               iwk30=iwk29+inmx
               iwk31=iwk30+inmx
               iwk32=iwk31+inmx
               inmx=(jdim-1)*(kdim-1)*(idim-1)*2
               iwk33=iwk32+inmx
               inmx=(jdim+1)*(kdim+1)*(idim+1)
               itemp=iwk33+inmx
            end if
         end if
c
         need      = max(need,itemp)
         ii = mblk2nd(nbl)
         mem_req_node(ii) = max(mem_req_node(ii),itemp)
         write(66,*)' 3(turbs)        itemp, need = ',itemp,need
         write(66,*)'       parallel: itemp, need = ',itemp,
     .   mem_req_node(ii)
         itempt = itemp
c
         nroomf = 0
         nroomi = 0
c
         if (ivisc(1).ge.4 .or. ivisc(2).ge.4 .or. ivisc(3).ge.4) then
c
c        space for recursive-box minimum distance routine,
c        used only for advanced turb models
c        nroomf is work space needed for floating pt values,
c        nroomi is for integer values
c
         call cntsurf(nsurf,maxbl0,maxgr0,maxseg0,ngrid,nblg,
     .                nbci0,nbcj0,nbck0,nbcidim,nbcjdim,
     .                nbckdim,ibcinfo,jbcinfo,kbcinfo,2004)
         minbox=sqrt(float(nsurf))
         minbox=max(minbox,50)
         nbb = 3*nsurf/minbox
         do 2010 igrid = 1,ngrid
         nbll = nblg(igrid)
         ntempf =  9*nsurf + 7*nbb
     .          +  4*jdimg(nbll)*kdimg(nbll)*idimg(nbll)
         if (ivmx.eq.4) then
            ntempi = 15*nsurf + 2*nbb
     .             +    jdimg(nbll)*kdimg(nbll)*idimg(nbll)
         else
            ntempi = 11*nsurf + 2*nbb
     .             +    jdimg(nbll)*kdimg(nbll)*idimg(nbll)
         end if
         nroomf = max(nroomf,ntempf)
         nroomi = max(nroomi,ntempi)
 2010    continue
c        additional memory allocation for temporary storage of smin at
c        grid points on all levels for blocks stored on same node
c        as current block nbl (mpi version)
         nadd_node = 0
         do 2225 nbll=1,nblock
         if (mblk2nd(nbll).eq.mblk2nd(nbl)) then
            nadd_node = nadd_node + jdimg(nbll)*kdimg(nbll)*idimg(nbll)
         end if
 2225    continue
c        additional memory allocation for temporary storage of smin at
c        grid points on all levels for all blocks (sequential version)
         nadd   = 0
         do 2226 nbll= 1,nblock
         nadd = nadd + jdimg(nbll)*kdimg(nbll)*idimg(nbll)
 2226    continue
         nroomf0 = nroomf
         nroomf = nroomf0 + nadd
         nroomf_node = nroomf0 + nadd_node
c
c        nworkxs, nworkixs needed to store surface data from
c        subroutines getpts and getptsbb; nworkxs and nworkixs
c        must be compatable with dimension statements for arrays
c        xs and ixs in subroutines getpts, getptsbb, collect_surf,
c        and collect_surfbb
c
         nworkxs  = 3*nsurf
         if (ivmx.eq.4) then
            nworkixs = 4*nsurf
         else
            nworkixs = 0
         end if
c        additional integer work space
         nwklsmin = maxbl0
         nwkireq1 = maxbl0*maxseg0*6
         nwkireq2 = maxbl0*maxseg0*6
         nroom       = nroomf      + nworkxs
         nroom_node  = nroomf_node + nworkxs
         nroomi      = nroomi      + nworkixs + nwklsmin
     .               + nwkireq1    + nwkireq2
         nroomi_node = nroomi      + nworkixs + nwklsmin
     .               + nwkireq1    + nwkireq2
         itemp      = nroom
         itemp_node = nroom_node
         need      = max(need,itemp)
         needi     = max(needi,nroomi)
         needi_node= max(needi_node,nroomi_node)
         ii = mblk2nd(nbl)
         mem_req_node(ii) = max(mem_req_node(ii),itemp_node)
            write(66,*)' 4(findmin_new)  itemp, need = ',itemp,need
         write(66,*)'       parallel: itemp, need = ',itemp_node,
     .   mem_req_node(ii)
         end if
c
         nv = 35
         if (idefrm(nbl) .gt. 0) nv=41
         maxmem = max(jdim*kdim*nv,jdim*idim*nv)
         itemp  = maxmem + ltot(nbl)
         itempt = itempt + ltot(nbl)
         itemp  = max(itemp,itempt)
         need   = max(need,itemp)
         ii = mblk2nd(nbl)
         mem_req_node(ii) = max(mem_req_node(ii),itemp)
         write(66,*)' 6(resid)        itemp, need = ',itemp,need
         write(66,*)'       parallel: itemp, need = ',itemp,
     .   mem_req_node(ii)
c
         if (ntime.eq.nit1 .and. level.ne.levb) then
            if (level.le.lglobal) then
               itemp = lwj0+jdim*kdim*idim*5
               need  = max(need,itemp)
               ii = mblk2nd(nbl)
               mem_req_node(ii) = max(mem_req_node(ii),itemp)
               write(66,*)' 7(collq)        itemp, need = ',itemp,need
               write(66,*)'       parallel: itemp, need = ',itemp,
     .         mem_req_node(ii)
            else
               itemp = lwj0+jdim*kdim*idim*5
               need  = max(need,itemp)
               ii = mblk2nd(nbl)
               mem_req_node(ii) = max(mem_req_node(ii),itemp)
               write(66,*)' 8(coll2q)       itemp, need = ',itemp,need
               write(66,*)'       parallel: itemp, need = ',itemp,
     .         mem_req_node(ii)
            end if
         end if
c
 6501    continue
      end if
 2000 continue
c
      if (level.eq.lglobal) go to 671
 2222 if (kxpand.eq.-1) go to 2012
      if (level.eq.levb) go to 7000
c
      level  = level-1
      ntime  = 0
      nit    = mit(level-levb+1,iseq)
      if (level.gt.levb) go to 9000 
      nit    = mit(1,iseq) + mtt
      kxpand = -1
      go to 9000
 2012 continue
c
      if (mgflag.ne.0) then
c
      do 7020 nbl=1,nblock
      if (level.ne.levelg(nbl)) go to 7020
      if (level.lt.lglobal) then
         iwk1  = 1
         iwk2  = iwk1 + jdim*kdim*idim*5
         iwk3  = iwk2 + jj2*kk2*ii2*5
         itemp = iwk3 + jdim*kk2*ii2*5
         need  = max(need,itemp)
         ii = mblk2nd(nbl)
         mem_req_node(ii) = max(mem_req_node(ii),itemp)
         write(66,*)' 9(addx)         itemp, need = ',itemp,need
         write(66,*)'       parallel: itemp, need = ',itemp,
     .   mem_req_node(ii)
      else
         if (mgflag.gt.1) then
            jc = jdimg(nbl)
            kc = kdimg(nbl)
            ic = idimg(nbl)
            do 9638 nblz=1,nblock
            if (nblz.eq.nbl) go to 9638
            nblcc = nblcg(nblz)
            if (nblcc.eq.nbl) then
               igridc = igridg(nblz)
                  if (iemg(igridc).gt.0) then
                     call lead(nblz,lw,lw2,maxbl0)
                     iwk1  = 1
                     iwk2  = iwk1 + jdim*kdim*idim*5
                     iwk3  = iwk2 + jc*kc*ic*5
                     iwk4  = iwk3 + jdim*kc*ic
                     itemp = iwk4 + jdim*kdim*ic
                     need  = max(need,itemp)
                     ii = mblk2nd(nblz)
                     mem_req_node(ii) = max(mem_req_node(ii),itemp)
                     write(66,*)'10(add2x)        itemp, need = ',
     .               itemp,need
                     write(66,*)'       parallel: itemp, need = ',
     .               itemp,mem_req_node(ii)
                  end if
            end if
 9638    continue
         end if 
      end if
 7020 continue
      end if
c
      level = level + 1
      if (level.eq.levt) go to 7000
      nit   = mtt
      if (nsm(level-levb+1).gt.0) then
         nit    = mit(level-levb+1,iseq)+mtt
         kxpand = 1
      end if
      ntime = 0
      if (nit.eq.0) go to 2012
      go to 9000
  671 continue
 8000 continue
 8500 continue
c
      do 1007 igrid=1,ngrid
      nbl  = nblg(igrid)
      call lead(nbl,lw,lw2,maxbl0)
      itemp = jdim*kdim*15
      need  = max(need,itemp)
      ii = mblk2nd(nbl)
      mem_req_node(ii) = max(mem_req_node(ii),itemp)
      write(66,*)' 11(cellvol)     itemp, need = ',itemp,need
      write(66,*)'       parallel: itemp, need = ',itemp,
     .mem_req_node(ii)
 1007 continue
c
      do 2007 igrid=1,ngrid
      nbl  = nblg(igrid)
      call lead(nbl,lw,lw2,maxbl0)
      itemp = (jdim-1)*(kdim-1)*(idim-1)*4
      need  = max(need,itemp)
      ii = mblk2nd(nbl)
      mem_req_node(ii) = max(mem_req_node(ii),itemp)
      write(66,*)' 12(dird)        itemp, need = ',itemp,need
      write(66,*)'       parallel: itemp, need = ',itemp,
     .mem_req_node(ii)
 2007 continue
c
      do 3000 igrid=1,ngrid
      nbl  = nblg(igrid)
      call lead(nbl,lw,lw2,maxbl0)
      ndim = max(kdim,idim)
      if (ifdsg(nbl,2).eq.0.or.idiag(2).eq.0) then
         itemp = ltot(nbl)   + jdim*ndim*145
         need  = max(need,itemp)
         ii = mblk2nd(nbl)
         mem_req_node(ii) = max(mem_req_node(ii),itemp)
         write(66,*)' 13(5x5 J)       itemp, need = ',itemp,need
         write(66,*)'       parallel: itemp, need = ',itemp,
     .   mem_req_node(ii)
      else
         itemp = ltot(nbl)   + jdim*ndim*35
         need  = max(need,itemp)
         ii = mblk2nd(nbl)
         mem_req_node(ii) = max(mem_req_node(ii),itemp)
         write(66,*)' 14(diagonal J)  itemp, need = ',itemp,need
         write(66,*)'       parallel: itemp, need = ',itemp,
     .   mem_req_node(ii)
      end if
c
      if (ifdsg(nbl,3).eq.0.or.idiag(3).eq.0) then
         itemp = ltot(nbl)   + jdim*ndim*145
         need  = max(need,itemp)
         ii = mblk2nd(nbl)
         mem_req_node(ii) = max(mem_req_node(ii),itemp)
         write(66,*)' 15(5x5 K)       itemp, need = ',itemp,need
         write(66,*)'       parallel: itemp, need = ',itemp,
     .   mem_req_node(ii)
      else
         itemp = ltot(nbl)   + jdim*ndim*35
         need  = max(need,itemp)
         ii = mblk2nd(nbl)
         mem_req_node(ii) = max(mem_req_node(ii),itemp)
         write(66,*)' 16(diagonal K)  itemp, need = ',itemp,need
         write(66,*)'       parallel: itemp, need = ',itemp,
     .   mem_req_node(ii)
      end if
c
      if (abs(i2d).ne.1) then
      if (ifdsg(nbl,1).eq.0.or.idiag(1).eq.0) then
         itemp = ltot(nbl)   + jdim*ndim*145
         need  = max(need,itemp)
         ii = mblk2nd(nbl)
         mem_req_node(ii) = max(mem_req_node(ii),itemp)
         write(66,*)' 17(5x5 I)       itemp, need = ',itemp,need
         write(66,*)'       parallel: itemp, need = ',itemp,
     .   mem_req_node(ii)
      else
         itemp = ltot(nbl)   + jdim*ndim*35
         need  = max(need,itemp)
         ii = mblk2nd(nbl)
         mem_req_node(ii) = max(mem_req_node(ii),itemp)
         write(66,*)' 18(diagonal I)  itemp, need = ',itemp,need
         write(66,*)'       parallel: itemp, need = ',itemp,
     .   mem_req_node(ii)
c
      end if
      end if
c
      itemp = ltot(nbl) + jdim*kdim*6 + jdim*kdim*idim*5
      need  = max(need,itemp)
      ii = mblk2nd(nbl)
      mem_req_node(ii) = max(mem_req_node(ii),itemp)
      write(66,*)' 19(metric)      itemp, need = ',itemp,need
      write(66,*)'       parallel: itemp, need = ',itemp,
     .mem_req_node(ii)
c
      iuns = max(itrans(nbl),irotat(nbl),idefrm(nbl))
      if (iuns.gt.0) then
c
         itemp = ltot(nbl) + jdim*kdim*idim*3 + jdim*kdim*3*2
     .         + kdim*idim*3*2 + jdim*idim*3*2
         need  = max(need,itemp)
         ii = mblk2nd(nbl)
         mem_req_node(ii) = max(mem_req_node(ii),itemp)
         write(66,*)' 20(trans)       itemp, need = ',itemp,need
         write(66,*)'       parallel: itemp, need = ',itemp,
     .   mem_req_node(ii)
c
         itemp = ltot(nbl) + 2*jdim*kdim*idim*3 + jdim*kdim*3*2 
     .         + kdim*idim*3*2 + jdim*idim*3*2
         need  = max(need,itemp)
         ii = mblk2nd(nbl)
         mem_req_node(ii) = max(mem_req_node(ii),itemp)
         write(66,*)' 21(rotate)      itemp, need = ',itemp,need
         write(66,*)'       parallel: itemp, need = ',itemp,
     .   mem_req_node(ii)
c
         itemp = ltot(nbl) + 2*jdim*kdim*idim*3 + jdim*kdim*6 
     .         + jdim*kdim*idim*5 + jdim*kdim*3*2
     .         + kdim*idim*3*2 + jdim*idim*3*2
         need  = max(need,itemp)
         ii = mblk2nd(nbl)
         mem_req_node(ii) = max(mem_req_node(ii),itemp)
         write(66,*)' 22(tmetric)     itemp, need = ',itemp,need
         write(66,*)'       parallel: itemp, need = ',itemp,
     .   mem_req_node(ii)
c
         itemp = ltot(nbl) + jdim*kdim*idim*3 + jdim*kdim*3*2
     .         + kdim*idim*3*2 + jdim*idim*3*2
         need  = max(need,itemp)
         ii = mblk2nd(nbl)
         mem_req_node(ii) = max(mem_req_node(ii),itemp)
         write(66,*)' 23(xtbatb)      itemp, need = ',itemp,need
         write(66,*)'       parallel: itemp, need = ',itemp,
     .   mem_req_node(ii)
      end if 
c
 3000 continue
c
c     requirements for plot3d output
c
      do n=1,nplot3d
         nbl = inpl3d(n,1)
         i1 = inpl3d(n,3)
         i2 = inpl3d(n,4)
         i3 = inpl3d(n,5)
         j1 = inpl3d(n,6)
         j2 = inpl3d(n,7)
         j3 = inpl3d(n,8)
         k1 = inpl3d(n,9)
         k2 = inpl3d(n,10)
         k3 = inpl3d(n,11)
         if (inpl3d(n,2).le.0) then
c           plot3d requirements
            call lead(nbl,lw,lw2,maxbl0)
            jdw = (j2-j1)/j3 + 1
            kdw = (k2-k1)/k3 + 1
            idw = (i2-i1)/i3 + 1
            itemp = jdw*kdw*idw*14 + jdim*kdim*idim*2
         else if (inpl3d(n,2).eq.1 .or. inpl3d(n,2).gt.2) then
c           plot3c requirements
            call lead(nbl,lw,lw2,maxbl0)
            i2 = min(idim-1,i2)
            j2 = min(jdim-1,j2)
            k2 = min(kdim-1,k2)
            i1 = min(idim-1,i1)
            j1 = min(jdim-1,j1)
            k1 = min(kdim-1,k1)
            jdw = (j2-j1)/j3 + 1
            kdw = (k2-k1)/k3 + 1
            idw = (i2-i1)/i3 + 1
            itemp = jdw*kdw*idw*9 + jdim*kdim*idim
         else
c           plot3t requirements
            if (ivmx .gt. 1) then
               call lead(nbl,lw,lw2,maxbl0)
               i2 = min(idim-1,i2)
               j2 = min(jdim-1,j2)
               k2 = min(kdim-1,k2)
               i1 = min(idim-1,i1)
               j1 = min(jdim-1,j1)
               k1 = min(kdim-1,k1)
               jdw = (j2-j1)/j3 + 1
               kdw = (k2-k1)/k3 + 1
               idw = (i2-i1)/i3 + 1
               if (ivmx.eq.8 .or. ivmx.eq.9 .or. ivmx.ge.11) then
                  itemp = jdw*kdw*idw*9 + jdim*kdim*idim + jdim*kdim*9
               else
                  itemp = jdw*kdw*idw*9 + jdim*kdim*idim 
     .                  + (jdim-1)*(kdim-1)*(idim-1)*9 + jdim*kdim*9
               end if
            end if
         end if
         need  = max(need,itemp)
         ii = mblk2nd(nbl)
         mem_req_node(ii) = max(mem_req_node(ii),itemp)
         write(66,*)' 23(plot3d)      itemp, need = ',itemp,need
         write(66,*)'       parallel: itemp, need = ',itemp,
     .   mem_req_node(ii)
      end do
c
c     requirements for print output
c
      do n=1,nprint
         nbl = inpr(n,1)
         i1 = inpr(n,3)
         i2 = inpr(n,4)
         i3 = inpr(n,5)
         j1 = inpr(n,6)
         j2 = inpr(n,7)
         j3 = inpr(n,8)
         k1 = inpr(n,9)
         k2 = inpr(n,10)
         k3 = inpr(n,11)
         if (inpr(n,2).eq.0) then
c           plot3d requirements
            call lead(nbl,lw,lw2,maxbl0)
            jdw = (j2-j1)/j3 + 1
            kdw = (k2-k1)/k3 + 1
            idw = (i2-i1)/i3 + 1
            itemp = jdw*kdw*idw*14 + jdim*kdim*idim*2
         else
c           plot3c requirements
            call lead(nbl,lw,lw2,maxbl0)
            i2 = min(idim-1,i2)
            j2 = min(jdim-1,j2)
            k2 = min(kdim-1,k2)
            i1 = min(idim-1,i1)
            j1 = min(jdim-1,j1)
            k1 = min(kdim-1,k1)
            jdw = (j2-j1)/j3 + 1
            kdw = (k2-k1)/k3 + 1
            idw = (i2-i1)/i3 + 1
            itemp = jdw*kdw*idw*9 + jdim*kdim*idim
         end if
         need  = max(need,itemp)
         ii = mblk2nd(nbl)
         mem_req_node(ii) = max(mem_req_node(ii),itemp)
         write(66,*)' 23(print)       itemp, need = ',itemp,need
         write(66,*)'       parallel: itemp, need = ',itemp,
     .   mem_req_node(ii)
      end do
c
c     determine i/o buffer size - defaults to the value that is
c     used for this sizing routine, which should be large enough 
c     for all but the baldwin-lomax output file
c
      libufdim = ibufdim0
c
      do 3005 igrid=1,ngrid
      nbl  = nblg(igrid)
      call lead(nbl,lw,lw2,maxbl0)
      n14(nbl) = 0
      ivisc(1) = iviscg(nbl,1)
      ivisc(2) = iviscg(nbl,2)
      ivisc(3) = iviscg(nbl,3)
      if (ivisc(3).eq.2 .or. ivisc(3).eq.3 .or. ivisc(2).eq.2 .or.
     .    ivisc(2).eq.3 .or. ivisc(1).eq.2 .or. ivisc(1).eq.3)then
c        determine i/o buffer size needed for Baldwin-Lomax
         ipw = 2
         call cntblmx(nbl,jdim,kdim,idim,iovrlp(nbl),ipw,n14(nbl),
     .                maxbl0,maxseg0,nbcj0,nbck0,nbci0,nbcjdim,
     .                nbckdim,nbcidim,jbcinfo,kbcinfo,
     .                ibcinfo)
c        need to add 1 to avoid spurious warning for last output entry
         n14_plus = n14(nbl) + 1
         libufdim = max(n14_plus,libufdim)
      end if
 3005 continue
c
c     the following memory allocation is needed for asynchronous
c     message passing of the 1-1, patch, periodic, and embeded 
c     interface data. note that this memory has no correspondance
c     in the sequential build of the code.
c
      icount_pat  = 0
      icount_blk  = 0
      icount_prd  = 0
      icount_emb  = 0
      icount_dpat = 0
c     icount_rad...placeholder for future bc2006 implementation
      icount_rad  = 0
      do levl = 1,levt
         is_blk(levl) = icount_blk + 1
         is_pat(levl) = icount_pat + 1
         is_prd(levl) = icount_prd + 1
         is_emb(levl) = icount_emb + 1
         ie_blk(levl) = is_blk(levl) - 1
         ie_pat(levl) = is_pat(levl) - 1
         ie_prd(levl) = is_prd(levl) - 1
         ie_emb(levl) = is_emb(levl) - 1
         do 6909 nbl=1,nblock
            if (levl.ne.levelg(nbl)) go to 6909
            icount_pat1 = icount_pat
            call pre_patch(nbl,lw,icount_pat,ninter,
     .                     iindex,intmax0,nsub10,isav_pat,
     .                     isav_pat_b,jjmax1,kkmax1,
     .                     iiint1,iiint2,maxbl0,jdimg,kdimg,idimg,
     .                     ierrflg)
            icount_blk1 = icount_blk
            call pre_blockbc(nbl,lw,icount_blk,idimg,
     .                       jdimg,kdimg,isav_blk,nblk,nbli,
     .                       limblk,isva,nblon,mxbli0,nou,bou,nbuf0,
     .                       ibufdim0,myid,maxbl0,ierrflg)
            icount_prd1 = icount_prd
            call pre_period(nbl,lw,lw2,icount_prd,
     .                      maxbl0,maxseg0,lbcprd0,
     .                      nbcj0,nbck0,nbci0,nbcjdim,nbckdim,
     .                      nbcidim,jbcinfo,kbcinfo,ibcinfo,
     .                      igridg,jdimg,kdimg,idimg,isav_prd,
     .                      is_prd,ie_prd,nbcprd,nou,bou,nbuf0,ibufdim0,
     .                      bcvali,bcvalj,bcvalk,myid,nblg,maxgr0,
     .                      ierrflg)
            icount_emb1 = icount_emb
            call pre_embed(nbl,lw,lw2,icount_emb,
     .                     maxbl0,maxseg0,lbcemb0,
     .                     nbcj0,nbck0,nbci0,nbcjdim,nbckdim,
     .                     nbcidim,jbcinfo,kbcinfo,ibcinfo,iemg,
     .                     igridg,jdimg,kdimg,idimg,isav_emb,
     .                     is_emb,ie_emb,nbcemb,nou,bou,nbuf0,ibufdim0,
     .                     myid,maxgr0,ierrflg)
            if (iunst.gt.0) then
c              pre_patch for dynamic interfaces
               call pre_patch(nbl,lw,icount_dpat,nintr,
     .                        iindx,intmx0,msub10,isav_dpat,
     .                        isav_dpat_b,jjmax1,kkmax1,
     .                        iiint1,iiint2,maxbl0,jdimg,kdimg,idimg,
     .                        ierrflg)
            end if
            if (icount_blk .gt. icount_blk1) then
               ie_blk(levl) = icount_blk
            end if
            if (icount_pat .gt. icount_pat1) then
               ie_pat(levl) = icount_pat
            end if
            if (icount_prd .gt. icount_prd1) then
               ie_prd(levl) = icount_prd
            end if
            if (icount_emb .gt. icount_emb1) then
               ie_emb(levl) = icount_emb
            end if
 6909    continue
      end do
c
      llbcprd = max(icount_prd,1)
      llbcemb = max(icount_emb,1)
      llbcrad = max(icount_rad,1)
c
      if (icount_blk.gt.0  .or. icount_pat.gt.0 .or.
     .    icount_prd.gt.0  .or. icount_emb.gt.0 .or.
     .    icount_dpat.gt.0 .or. icount_rad.gt.0) then 
         write(66,*)' -async are requirements for',
     .   ' asynchronous message passing, by node'
      end if

c 
c     for 1-1 interfaces
c
      if (icount_blk .gt. 0) then
         do myid0 = 1,nodel
            itemp = 1
c
c           for wk(kqintl) allocation
            do lcnt = 1,icount_blk
c              ic_blk is current (to) block
c              in_blk is neighbor (from) block
               ic_blk  = isav_blk(lcnt,4)
               in_blk  = isav_blk(lcnt,5)
               nd_dest = mblk2nd(ic_blk)
               nd_srce = mblk2nd(in_blk)
               if (nd_srce.ne.myid0) then
                  if (nd_dest.eq.myid0) then
                  n = isav_blk(lcnt,1)
                  jface = isav_blk(lcnt,6)
                  idimn = idimg(in_blk)
                  jdimn = jdimg(in_blk)
                  kdimn = kdimg(in_blk)
                  if (jface.eq.1) maxdims = jdimn*kdimn
                  if (jface.eq.2) maxdims = kdimn*idimn
                  if (jface.eq.3) maxdims = jdimn*idimn
                     itemp = itemp + maxdims*13
                     if (ivmx .ge. 2)
     .                  itemp = itemp + maxdims*2
                     if (ivmx .ge. 4)
     .                  itemp = itemp + maxdims*4
                     if (ivmx .ge. 1)
c                       for cell volumes
     .                  itemp = itemp + maxdims*2
                  end if
               end if
            end do
c
c           for wk(ktl) allocation
c
            iskipz = 0
            if (iskipz .eq. 0) then
c
            do lcnt = 1,icount_blk
c              ic_blk is current (to) block
c              in_blk is neighbor (from) block
               ic_blk  = isav_blk(lcnt,4)
               in_blk  = isav_blk(lcnt,5)
               nd_dest = mblk2nd(ic_blk)
               nd_srce = mblk2nd(in_blk)
               if (nd_dest.ne.myid0) then
                  if (nd_srce.eq.myid0) then
                     jface = isav_blk(lcnt,6)
                     idimn  = idimg(in_blk)
                     jdimn  = jdimg(in_blk)
                     kdimn  = kdimg(in_blk)
                     if (jface.eq.1) maxdims = jdimn*kdimn
                     if (jface.eq.2) maxdims = kdimn*idimn
                     if (jface.eq.3) maxdims = jdimn*idimn
                     itemp = itemp + maxdims*13
                     if (ivmx .ge. 2)
     .                  itemp = itemp + maxdims*2
                     if (ivmx .ge. 4)
     .                  itemp = itemp + maxdims*4
                     if (ivmx .ge. 1)
c                       for cell volumes
     .                  itemp = itemp + maxdims*2
                  end if
               end if
            end do
c
            end if
c
            mem_req_node(myid0) = max(mem_req_node(myid0),itemp)
            write(66,*)
     .      ' 22(1:1-async)   itemp, need = ',itemp,mem_req_node(myid0)
         end do
      end if
c
c     for patched interfaces
c
      if (icount_pat .gt. 0) then
         do myid0 = 1,nodel
            itemp = 1
c
c           for wk(kqintl) allocation
c           itb is "to" block
c           ifb is "from" block
            do intr = 1,icount_pat
               lmax_val = isav_pat(intr,2)
               do nf = 1,lmax_val
                  itb  = isav_pat(intr,1)
                  ifb  = isav_pat_b(intr,nf,1)
                  ityp = isav_pat_b(intr,nf,2)
                  if (mblk2nd(ifb).ne.myid0) then
                     if (mblk2nd(itb).eq.myid0) then
                        if (ityp .eq.1) then
                           jmax1  = jdimg(ifb)
                           kmax1  = kdimg(ifb)
                        else if (ityp .eq.2) then
                           jmax1  = kdimg(ifb)
                           kmax1  = idimg(ifb)
                        else if (ityp .eq.3) then
                           jmax1  = jdimg(ifb)
                           kmax1  = idimg(ifb)
                        end if
                        itemp = itemp + jmax1*kmax1*5*2
                        if (ivmx .ge. 2)
     .                     itemp = itemp + jmax1*kmax1*1*2
                        if (ivmx .ge. 4)
     .                     itemp = itemp + jmax1*kmax1*2*2
                     end if
                  end if
               end do
            end do
c
c           for wk(ktl) allocation
c           itb is "to" block
c           ifb is "from" block
            do intr = 1,icount_pat
               lmax_val = isav_pat(intr,2)
               do nf = 1,lmax_val
                  itb  = isav_pat(intr,1)
                  ifb  = isav_pat_b(intr,nf,1)
                  ityp = isav_pat_b(intr,nf,2)
                  if (mblk2nd(ifb).eq.myid0) then
                     if (mblk2nd(itb).ne.myid0) then
                        if (ityp .eq.1) then
c                          i=const. patch
                           jmax1  = jdimg(ifb)
                           kmax1  = kdimg(ifb)
                        else if (ityp .eq.2) then
c                          j=const. patch
                           jmax1  = kdimg(ifb)
                           kmax1  = idimg(ifb)
                        else if (ityp .eq.3) then
c                          k=const. patch
                           jmax1  = jdimg(ifb)
                           kmax1  = idimg(ifb)
                        end if
                        itemp = itemp + jmax1*kmax1*5*2
                        if (ivmx .ge. 2)
     .                     itemp = itemp + jmax1*kmax1*1*2
                        if (ivmx .ge. 4)
     .                     itemp = itemp + jmax1*kmax1*2*2
                     end if
                  end if
               end do
            end do
c
c           for subroutine int2
c           itb is "to" block
c           ifb is "from" block
            itemp1 = 0
            do intr = 1,icount_pat
               lmax_val = isav_pat(intr,2) 
               do nf = 1,lmax_val
                  itb  = isav_pat(intr,1)
                  ifb  = isav_pat_b(intr,nf,1)
                  ityp = isav_pat_b(intr,nf,2)
                  if (ityp .eq.1) then
                     jmax1  = jdimg(ifb)
                     kmax1  = kdimg(ifb)
                  else if (ityp .eq.2) then
                     jmax1  = kdimg(ifb)
                     kmax1  = idimg(ifb)
                  else if (ityp .eq.3) then
                     jmax1  = jdimg(ifb)
                     kmax1  = idimg(ifb)
                  end if
                  itemp1 = max(itemp1,4*jmax1*kmax1)
               end do
            end do
c
            itemp = itemp + itemp1
            mem_req_node(myid0) = max(mem_req_node(myid0),itemp)
            write(66,*)
     .      ' 23(patch-async) itemp, need = ',itemp,mem_req_node(myid0)
         end do
      end if
c
c     for periodic interfaces
c
      if (icount_prd .gt. 0) then
         do myid0 = 1,nodel
            itemp = 1
c
c           for wk(kqintl) allocation
c           nbll is current block
c           nblp is periodic block
            do lcnt = 1,icount_prd
               nbll    = isav_prd(lcnt,1)
               nseg    = isav_prd(lcnt,11)
               nface   = isav_prd(lcnt,2)
               ldata   = lwdat(nbll,nseg,nface)
               nblp    = isav_prd(lcnt,12)
               nd_recv = mblk2nd(nbll)
               nd_srce = mblk2nd(nblp)
               if (nd_recv.eq.myid0) then
                  if (nd_srce.ne.myid0) then
                     idimp = idimg(nblp)
                     jdimp = jdimg(nblp)
                     kdimp = kdimg(nblp)
                     if (nface.eq.1.or.nface.eq.2) maxdims = jdimp*kdimp
                     if (nface.eq.3.or.nface.eq.4) maxdims = kdimp*idimp
                     if (nface.eq.5.or.nface.eq.6) maxdims = jdimp*idimp
                     itemp = itemp + maxdims*13
                     if (ivmx .ge. 2)
     .                  itemp = itemp + maxdims*2
                     if (ivmx .ge. 4)
     .                  itemp = itemp + maxdims*4
                  end if
               end if
            end do
c
c           for wk(ktl) allocation
c           nbll is current block
c           nblp is periodic block
            do lcnt = 1,icount_prd
               nbll    = isav_prd(lcnt,1)
               nblp    = isav_prd(lcnt,12)
               nseg    = isav_prd(lcnt,11)
               nface   = isav_prd(lcnt,2)
               ldata   = lwdat(nbll,nseg,nface)
               nd_recv = mblk2nd(nbll)
               nd_srce = mblk2nd(nblp)
               if (nd_srce.eq.myid0) then
                  if (nd_recv.ne.myid0) then
                     nface  = isav_prd(lcnt,2)
                     ista   = isav_prd(lcnt,3)
                     iend   = isav_prd(lcnt,4)
                     jsta   = isav_prd(lcnt,5)
                     jend   = isav_prd(lcnt,6)
                     ksta   = isav_prd(lcnt,7)
                     kend   = isav_prd(lcnt,8)
                     idimp  = idimg(nblp)
                     jdimp  = jdimg(nblp)
                     kdimp  = kdimg(nblp)
                     if (nface.eq.1.or.nface.eq.2) maxdims = jdimp*kdimp
                     if (nface.eq.3.or.nface.eq.4) maxdims = kdimp*idimp
                     if (nface.eq.5.or.nface.eq.6) maxdims = jdimp*idimp
                     itemp = itemp + maxdims*13
                     if (ivmx .ge. 2)
     .                  itemp = itemp + maxdims*2
                     if (ivmx .ge. 4)
     .                  itemp = itemp + maxdims*4
                  end if
               end if
            end do
c
            mem_req_node(myid0) = max(mem_req_node(myid0),itemp)
            write(66,*)
     .      ' 24(period-async)itemp, need = ',itemp,mem_req_node(myid0)
         end do
      end if
c
c     for embedded interfaces
c
      if (icount_emb .gt. 0) then
         do myid0 = 1,nodel
            itemp = 1
c
c           for wk(kqintl) allocation
c           nblf is finer (embedded) block
c           nblc is coarser block
            do lcnt = 1,icount_emb
               nblf    = isav_emb(lcnt,1)
               nface   = isav_emb(lcnt,2)
               nblc    = isav_emb(lcnt,9)
               nd_recv = mblk2nd(nblf)
               nd_srce = mblk2nd(nblc)
               if (nd_recv.eq.myid0) then
                  if (nd_srce.ne.myid0) then
                     idimc = idimg(nblc)
                     jdimc = jdimg(nblc)
                     kdimc = kdimg(nblc)
                     if (nface.eq.1.or.nface.eq.2) maxdims = jdimc*kdimc
                     if (nface.eq.3.or.nface.eq.4) maxdims = kdimc*idimc
                     if (nface.eq.5.or.nface.eq.6) maxdims = jdimc*idimc
                     itemp = itemp + maxdims*15
                     if (ivmx .ge. 2)
     .                  itemp = itemp + maxdims*3
                     if (ivmx .ge. 4)
     .                  itemp = itemp + maxdims*6
                  end if
               end if
            end do
c
c           for wk(ktl) allocation
c           nblf is finer (embedded) block
c           nblc is coarser block
            do lcnt = 1,icount_emb
               nblf    = isav_emb(lcnt,1)
               nblc    = isav_emb(lcnt,9)
               nface   = isav_emb(lcnt,2)
               nd_recv = mblk2nd(nblf)
               nd_srce = mblk2nd(nblc)
               if (nd_srce.eq.myid0) then
                  if (nd_recv.ne.myid0) then
                     is     = isav_emb(lcnt,3)
                     ie     = isav_emb(lcnt,4)
                     js     = isav_emb(lcnt,5)
                     je     = isav_emb(lcnt,6)
                     ks     = isav_emb(lcnt,7)
                     ke     = isav_emb(lcnt,8)
                     idimc  = idimg(nblc)
                     jdimc  = jdimg(nblc)
                     kdimc  = kdimg(nblc)
                     if (nface.eq.1.or.nface.eq.2) maxdims = jdimc*kdimc
                     if (nface.eq.3.or.nface.eq.4) maxdims = kdimc*idimc
                     if (nface.eq.5.or.nface.eq.6) maxdims = jdimc*idimc
                     itemp = itemp + maxdims*15
                     if (ivmx .ge. 2)
     .                  itemp = itemp + maxdims*3
                     if (ivmx .ge. 4)
     .                  itemp = itemp + maxdims*6
                  end if
               end if
            end do
c
            mem_req_node(myid0) = max(mem_req_node(myid0),itemp)
            write(66,*)
     .      ' 25(embed-async) itemp, need = ',itemp,mem_req_node(myid0)
         end do
      end if
c
c     for dynamic grids
c
c     itb is "to" block
c     ifb is "from" block
      if (iunst.gt.0) then
         do myid0=1,nodel
            itemp = 1
            do intr=1,icount_dpat
               lmax_val=isav_dpat(intr,2)
               do nf=1,lmax_val
                  itb = isav_dpat(intr,1)
                  ifb = isav_dpat_b(intr,nf,1)
                  mdim1  = isav_dpat_b(intr,nf,3)
                  ndim1  = isav_dpat_b(intr,nf,4)
                  if (mblk2nd(ifb).ne.myid0) then
                     if (mblk2nd(itb).eq.myid0) then
                        itemp = itemp + 3*mdim1*ndim1
                     end if
                  end if
               end do
            end do
         mem_req_node(myid0) = max(mem_req_node(myid0),itemp)
         write(66,*)
     .   ' 26(dyn-async)   itemp, need = ',itemp,mem_req_node(myid0)
         end do
      end if
c
c     also need space for int2 in sequential case
c
      if (ninter.gt.0) then
         need  = max(need,itemp1)
         write(66,*)' 27(int2)        itemp, need = ',itemp1,need
      end if
c
c     add space for timesave array in mgblk 
c
      need = need + lmaxbl
      do ii = 1,nodel
         mem_req_node(ii) = mem_req_node(ii) + lmaxbl
      end do
c
      mem_wk_max = mem_req_node(1)
      do myid0 = 2, nodel
         if (mem_req_node(myid0).gt.mem_wk_max)
     &       mem_wk_max = mem_req_node(myid0)
      end do
      tot_nod = mem_w_max + mem_wk_max
c
      write(66,*)
      write(66,*)'***********************************************'
      write(66,*)
      write(66,*)'      TEMPORARY STORAGE REQUIREMENTS - WK'
      write(66,*)'            SUMMARY FOR ALL BUILDS'
      write(66,*)
      write(66,*)'***********************************************'
c
      write(66,8436) real(float(need))/1.e+06
c
      write(66,8435)
      write(66,8426)
      fmemsum = 0
      nptsum  = 0
      do i = 1,nodel
         write(66,8427) i, no_of_points(i), mem_req_node(i)/1.e+06
         fmemsum = fmemsum + mem_req_node(i)/1.e+06
         nptsum  = nptsum + no_of_points(i)
         mem_req_node(i) = 0
      end do
      write(66,8429) npts_max,real(float(mem_wk_max))/1.e+06
c
 8435 format(/,1x,34hmemory (mw) for wk storage (nodes))
 8436 format(/,1x,42hmemory (mw) for wk storage (sequential) = ,f9.4)
c
c***********************************************************************
c
c     evaluate requirements for temporary integer array size 
c     (iwork/iwk array)
c
c***********************************************************************
c
c     in mgbl, to call global
      needi_glo = lmaxbl*3 + 11*lnplts
      needi  = max(needi,needi_glo)
c
c     in mgbl, to call trnsfr_vals
      needi_trn = 53*lmaxbl        + 117      + 9*lmaxbl*lmxsegdg
     .          + 48*lmaxbl*lmxseg + 19*lmxbli + 22*lnplts
     .          + 4*lmaxgr         + 9*lmaxcs
      needi  = max(needi,needi_trn)
c
c     in mgbl, to call pointers
      needi_pnt = lmaxbl*2
      needi  = max(needi,needi_pnt)
c
c     in mgbl, to call setup
      needi_set = 3*lmaxgr + lmaxbl*8 
      needi  = max(needi,needi_set)
c
c     in mgbl, to call mgblk
      needi_mgblk = lmaxbl*7*3 + lmaxbl*8
      needi  = max(needi,needi_mgblk)
c
c     in mgbl, to call qout
      needi_qout = 3*lnplts + lmaxbl
      needi  = max(needi,needi_qout)
c
c     in mgbl, to call yplusout
      needi_yout = 12*lmaxbl
      needi  = max(needi,needi_yout)
c
c     in setup, to call findmin_new (nroomi already calculated above)
      itempi = nroomi + needi_set
      needi  = max(needi,itempi)
c
c     in setslave, to call setcorner
      if (iunst.gt.1 .or. idef_ss.gt.0) then
         call cntsurf(ns2004,maxbl0,maxgr0,maxseg0,ngrid,nblg,
     .                nbci0,nbcj0,nbck0,nbcidim,nbcjdim,
     .                nbckdim,ibcinfo,jbcinfo,kbcinfo,2004)
         call cntsurf(ns1005,maxbl0,maxgr0,maxseg0,ngrid,nblg,
     .                nbci0,nbcj0,nbck0,nbcidim,nbcjdim,
     .                nbckdim,ibcinfo,jbcinfo,kbcinfo,1005)
         call cntsurf(ns1006,maxbl0,maxgr0,maxseg0,ngrid,nblg,
     .                nbci0,nbcj0,nbck0,nbcidim,nbcjdim,
     .                nbckdim,ibcinfo,jbcinfo,kbcinfo,1006)
         nsurf = ns2004 + ns1005 + ns1006
         itempi = nsurf + needi_set
         needi  = max(needi,itempi)
      end if
c
c     in mgblk, to output movie data
c 
      imovie = 0
      if (abs(movie).gt. 0) then     
         imovie = lmaxbl + 3*lnplts 
      end if
c
c     in mgblk, to call dynptch
      if (iunst .gt. 0) then
         itempi = 65*lmaxbl + 3*lintmx + needi_mgblk + imovie
         needi  = max(needi,itempi)
      end if
c
c     in mgblk, to call pre_bc
      itempi = 65*lmaxbl + 4*lnsub1 + needi_mgblk  + imovie
      needi  = max(needi,itempi)
c
c     in mgblk, to call qout
      itempi = lmaxbl + 3*lnplts + needi_mgblk + imovie
      needi  = max(needi,itempi)
c
c     in mgblk, to call resetg
      itempi = lmaxbl + needi_mgblk + imovie
      needi  = max(needi,itempi)
c
c     in mgblk, to call bc_period
      itempi = 30*llbcprd + needi_mgblk + imovie
      needi  = max(needi,itempi)
c
c     in mgblk, to call bc_embed
      itempi = 18*llbcemb + needi_mgblk + imovie
      needi  = max(needi,itempi)
c
c     in mgblk, to call bc_blkint
      itempi = 30*lmxbli + needi_mgblk + imovie
      needi  = max(needi,itempi)
c
c     in mgblk, to call bc_patch
      itempi = 21*lintmax*lnsub1 + needi_mgblk + imovie
      needi  = max(needi,itempi)
c
c     in qout, to call bc_period
      itempi = 30*llbcprd + needi_qout
      needi  = max(needi,itempi)
c
c     in qout, to call bc_embed
      itempi = 18*llbcemb + needi_qout
      needi  = max(needi,itempi)
c
c     in qout, to call bc_blkint
      itempi = 30*lmxbli + needi_qout
      needi  = max(needi,itempi)
c
c     in qout, to call bc_patch
      itempi = 21*lintmax*lnsub1 + needi_qout
      needi  = max(needi,itempi)
c
c     need to add 1 for safety
c
      needi = needi + 1
c
c     currently, mpi code requires same integer work space
c     as sequential code
c
      needi_node = needi
c
c***********************************************************************
c
c     set various work array sizes
c
c***********************************************************************
c
c     sequential build real and integer work array sizes
      isum = ftot_seq+need
      isumi = max(needi,1)
c     parallel build real and integer work array sizes (per compute node)
      isum_n = tot_nod
      isumi_n = max(needi_node,1)
c
c***********************************************************************
c
c     determine number of slave points for deforming mesh case
c
c***********************************************************************
c
      if (iunst.gt.1 .or. idef_ss.gt.0) then
      do inode = 1,nodel
         ivert = 0
         do nbl = 1,nblock
            iskp = iskip(nbl)
            jskp = jskip(nbl)
            kskp = kskip(nbl)
            if (mblk2nd(nbl).eq.inode) then
c              i=1 subfaces
               do nseg=1,nbci0(nbl)
                  ii     = 1
                  ibctyp = ibcinfo(nbl,nseg,1,ii)
                  if (abs(ibctyp).ne.2004 .or.
     .                abs(ibctyp).ne.1005 .or.
     .                abs(ibctyp).ne.1006) then
                      js = ibcinfo(nbl,nseg,2,ii)
                      je = ibcinfo(nbl,nseg,3,ii)
                      ks = ibcinfo(nbl,nseg,4,ii)
                      ke = ibcinfo(nbl,nseg,5,ii)
                      do j=js,je,jskp
                         do k=ks,ke,kskp
                            ivert = ivert + 1
                         end do
                      end do
                  end if
               end do
c              i=idim subfaces
               do nseg=1,nbcidim(nbl)
                  ii     = 2
                  ibctyp = ibcinfo(nbl,nseg,1,ii)
                  if (abs(ibctyp).ne.2004 .or.
     .                abs(ibctyp).ne.1005 .or.
     .                abs(ibctyp).ne.1006) then
                      js = ibcinfo(nbl,nseg,2,ii)
                      je = ibcinfo(nbl,nseg,3,ii)
                      ks = ibcinfo(nbl,nseg,4,ii)
                      ke = ibcinfo(nbl,nseg,5,ii)
                      do j=js,je,jskp
                         do k=ks,ke,kskp
                            ivert = ivert + 1
                         end do
                      end do
                  end if
               end do
               if (i2d .eq. 0) then
c                 j=1 subfaces
                  do nseg=1,nbcj0(nbl)
                     jj     = 1
                     jbctyp = jbcinfo(nbl,nseg,1,jj)
                     if (abs(jbctyp).ne.2004 .or.
     .                   abs(jbctyp).ne.1005 .or.
     .                   abs(jbctyp).ne.1006) then
                         is = jbcinfo(nbl,nseg,2,jj)
                         ie = jbcinfo(nbl,nseg,3,jj)
                         ks = jbcinfo(nbl,nseg,4,jj)
                         ke = jbcinfo(nbl,nseg,5,jj)
                         do i=is,ie,iskp
                            do k=ks,ke,kskp
                               ivert = ivert + 1
                            end do
                         end do
                     end if
                  end do
c                 j=jdim subfaces
                  do nseg=1,nbcjdim(nbl)
                     jj     = 2
                     jbctyp = jbcinfo(nbl,nseg,1,jj)
                     if (abs(jbctyp).ne.2004 .or.
     .                   abs(jbctyp).ne.1005 .or.
     .                   abs(jbctyp).ne.1006) then
                         is = jbcinfo(nbl,nseg,2,jj)
                         ie = jbcinfo(nbl,nseg,3,jj)
                         ks = jbcinfo(nbl,nseg,4,jj)
                         ke = jbcinfo(nbl,nseg,5,jj)
                         do i=is,ie,iskp
                            do k=ks,ke,kskp
                               ivert = ivert + 1
                            end do
                         end do
                     end if
                  end do
c                 k=1 subfaces
                  do nseg=1,nbck0(nbl)
                     kk     = 1
                     kbctyp = kbcinfo(nbl,nseg,1,kk)
                     if (abs(kbctyp).ne.2004 .or.
     .                   abs(kbctyp).ne.1005 .or.
     .                   abs(kbctyp).ne.1006) then
                         is = kbcinfo(nbl,nseg,2,kk)
                         ie = kbcinfo(nbl,nseg,3,kk)
                         js = kbcinfo(nbl,nseg,4,kk)
                         je = kbcinfo(nbl,nseg,5,kk)
                         do i=is,ie,iskp
                            do j=js,je,jskp
                               ivert = ivert + 1
                            end do
                         end do
                     end if
                  end do
c                 k=kdim subfaces
                  do nseg=1,nbckdim(nbl)
                     kk     = 2
                     kbctyp = kbcinfo(nbl,nseg,1,kk)
                     if (abs(kbctyp).ne.2004 .or.
     .                   abs(kbctyp).ne.1005 .or.
     .                   abs(kbctyp).ne.1006) then
                         is = kbcinfo(nbl,nseg,2,kk)
                         ie = kbcinfo(nbl,nseg,3,kk)
                         js = kbcinfo(nbl,nseg,4,kk)
                         je = kbcinfo(nbl,nseg,5,kk)
                         do i=is,ie,iskp
                            do j=js,je,jskp
                               ivert = ivert + 1
                            end do
                         end do
                     end if
                  end do
               end if
            end if
         end do
c
         lnslave = max(ivert,lnslave)
c
      end do
      end if
c
c     add any extra memory input by the user
c     (memadd, memaddi default to zero)
c
      if (nnodes.eq.1) then
         isum    = isum  + memadd
         isumi   = isumi + memaddi
      else
         isum_n  = isum_n  + memadd
         isumi_n = isumi_n + memaddi
      end if
c
c***********************************************************************
c
c     set the auxiliary array sizes and print them out
c
c***********************************************************************
c
      if (nnodes.eq.1) then
         mwork   = isum
         mworki  = isumi
      else
         mwork   = isum_n
         mworki  = isumi_n
      end if
      nplots   = lnplts
      minnode  = lminnode
      iitot    = liitot
      intmax   = lintmax
      maxxe    = lmaxxe
      mxbli    = lmxbli
      nsub1    = lnsub1
      lbcprd   = llbcprd
      lbcemb   = llbcemb
      lbcrad   = llbcrad
      maxbl    = lmaxbl
      maxgr    = lmaxgr
      maxseg   = lmxseg
      maxcs    = lmaxcs
      ncycmax  = lncycm
      intmx    = lintmx
      mxxe     = lmxxe
      mptch    = lmptch
      msub1    = lmsub1
      ibufdim  = libufdim
      nbuf     = nbuf0
      nmds     = lnmds
      maxaes   = lmaxaes
      nslave   = lnslave
      nmaster  = lnmaster
      maxsegdg = lmxsegdg
c
      write(66,*)
      write(66,*)'***********************************************'
      write(66,*)
      write(66,*)'    PARAMETER SIZES REQUIRED FOR THIS CASE:'
      write(66,*)
      write(66,*)'***********************************************'
      write(66,*)
      if (imode.eq.0) then
         write(66,'('' mwork    = '',i10,'' (sequential)'')') isum
         write(66,'('' mworki   = '',i10,'' (sequential)'')') isumi
         write(66,'('' mwork    = '',i10,'' (parallel)'')') isum_n
         write(66,'('' mworki   = '',i10,'' (parallel)'')') isumi_n
      else
         write(66,'('' mwork    = '',i10)') mwork
         write(66,'('' mworki   = '',i10)') mworki
      end if
      write(66,'('' nplots   = '',i10)') nplots
      write(66,'('' minnode  = '',i10)') minnode
      write(66,'('' iitot    = '',i10)') iitot
      write(66,'('' intmax   = '',i10)') intmax
      write(66,'('' maxxe    = '',i10)') maxxe
      write(66,'('' mxbli    = '',i10)') mxbli
      write(66,'('' nsub1    = '',i10)') nsub1
      write(66,'('' lbcprd   = '',i10)') lbcprd
      write(66,'('' lbcemb   = '',i10)') lbcemb
      write(66,'('' lbcrad   = '',i10)') lbcrad
      write(66,'('' maxbl    = '',i10)') maxbl
      write(66,'('' maxgr    = '',i10)') maxgr
      write(66,'('' maxseg   = '',i10)') maxseg
      write(66,'('' maxcs    = '',i10)') maxcs
      write(66,'('' ncycmax  = '',i10)') ncycmax
      write(66,'('' intmx    = '',i10)') intmx
      write(66,'('' mxxe     = '',i10)') mxxe
      write(66,'('' mptch    = '',i10)') mptch
      write(66,'('' msub1    = '',i10)') msub1
      write(66,'('' nmds     = '',i10)') nmds
      write(66,'('' maxaes   = '',i10)') maxaes
      write(66,'('' ibufdim  = '',i10)') ibufdim
      write(66,'('' nbuf     = '',i10)') nbuf
      write(66,'('' nslave   = '',i10)') nslave
      write(66,'('' nmaster  = '',i10)') nmaster
      write(66,'('' maxsegdg = '',i10)') maxsegdg
      write(66,*)
c
c     stop if any of the parameters are non-zero
c
      ichek = 1
c
      if (mwork   .le. 0 .or.
     .   mworki   .le. 0 .or.
     .   nplots   .le. 0 .or.
     .   minnode  .le. 0 .or.
     .   iitot    .le. 0 .or.
     .   intmax   .le. 0 .or.
     .   maxxe    .le. 0 .or.
     .   mxbli    .le. 0 .or.
     .   nsub1    .le. 0 .or.
     .   lbcprd   .le. 0 .or.
     .   lbcemb   .le. 0 .or.
     .   lbcrad   .le. 0 .or.
     .   maxbl    .le. 0 .or.
     .   maxgr    .le. 0 .or.
     .   maxseg   .le. 0 .or.
     .   maxcs    .le. 0 .or.
     .   ncycmax  .le. 0 .or.
     .   intmx    .le. 0 .or.
     .   mxxe     .le. 0 .or.
     .   mptch    .le. 0 .or.
     .   msub1    .le. 0 .or.
     .   ibufdim  .le. 0 .or.
     .   nbuf     .le. 0 .or.
     .   nmds     .le. 0 .or.
     .   maxaes   .le. 0 .or.
     .   nslave   .le. 0 .or.
     .   nmaster  .le. 0 .or.
     .   maxsegdg .le. 0) ichek = 0
      if (ichek.eq.0) then
         nou(1) = min(nou(1)+1,ibufdim0)
         write(bou(nou(1),1),'(''error in routine sizer - one of'',
     .                     ''the parameters listed above is zero'')')
         ierrflg = -99
         call termn8(myid,ierrflg,ibufdim0,nbuf0,bou,nou)
      end if
c
c***********************************************************************
c
c     output storage requirement summary
c
c***********************************************************************
c
      write(66,*)
      write(66,*)
     .'**************************************************************'
      write(66,*)
      write(66,*)'   SUMMARY OF STORAGE REQUIREMENTS - W + WK ARRAYS'
      write(66,*)
      write(66,*)' sequential version:'
      write(66,*)
      write(66,770) int(real(ftot_seq))
      write(66,775) need
      write(66,776) needi
c
  770 format('          permanent array w   requires ',i10,' (words)')
  775 format('          temporary array wk  requires ',i10,' (words)')
  776 format('          temporary array iwk requires ',i10,' (words)')
c
      write(66,*)
      write(66,*)' parallel version, per node:'
      write(66,*)
      write(66,770) mem_w_max
      write(66,775) mem_wk_max
      write(66,776) needi_node
c
      write(66,*)'  '
      write(66,*)'  '
      write(66,800) isum
 800  format(' >>> Estimate for mwork      (sequential)     = ',
     .i10,' <<<')
c
      write(66,1801) isumi
 1801 format(/,' >>> Estimate for mworki     (sequential)     = ',
     .i10,' <<<')
c
      write(66,871) isum_n
 871  format(/,' >>> Estimate for mwork  (per node, parallel) = ',
     .i10,' <<<')
c
      write(66,1803) isumi_n
 1803 format(/,' >>> Estimate for mworki (per node, parallel) = ',
     .i10,' <<<')
c
      write(66,872) nodel
 872  format(/,' >>> Parallel code sized for',i4,' nodes,',
     .' min. (+host)       <<<')
c
c     calculate best possible speedup compared to single-node case
c     (note: this assumes all nodes are of equal performance...i.e.
c     heterogeneous clusters excluded). Under this assumption, the
c     best possible speedup is simply a function of the number of
c     points in each block
c
      maxpt = 0
      nptt  = 0
      do ind=1,nodel
         npt  = no_of_points(ind)
         nptt = nptt + npt
         if (npt .gt. maxpt) maxpt = npt
      end do
c
      write(66,*)
      write(66,*)
     .'**************************************************************'
      write(66,*)
c
c     output ideal speedup vs no. of nodes
c
      if (imode .eq. 0) then
         open(67,file='ideal_speedup.dat',form='formatted',
     .   status='unknown')
         write(67,101)
 101     format(1x,'compute nodes  ideal speedup')
         do nnn=1,ngrid
            call compg2n(nblock,ngrid,ncgg,nblg,idimg,jdimg,kdimg,
     .                   nblcg,nnn,iwrk,myid,myhost,mblk2nd,
     .                   mycomm,maxgr0,maxbl0,ierrflg,ibufdim,
     .                   nbuf,bou,nou)
            do ii = 1,nnn
               mem_req_node(ii) = 0
               no_of_points(ii) = 0
            end do
            do i = 1,nblock
               npoints = (idimg(i)-1)*(jdimg(i)-1)*(kdimg(i)-1)
               nod = mblk2nd(i)
               mem_req_node(nod) = mem_req_node(nod) + memblock(i)
               no_of_points(nod) = no_of_points(nod) + npoints
            end do
            npts_max  = no_of_points(1)
            mem_w_max = mem_req_node(1)
            do i = 2,nnn
               if (mem_req_node(i).gt.mem_w_max)
     .            mem_w_max = mem_req_node(i)
               if (no_of_points(i).gt.npts_max )
     .            npts_max  = no_of_points(i)
            end do
            maxpt = 0
            nptt  = 0
            do ind=1,nnn
               npt  = no_of_points(ind)
               nptt = nptt + npt
               if (npt .gt. maxpt) maxpt = npt
            end do
            speedopt = float(nptt)/float(maxpt)
            write(67,102) nnn,real(speedopt)
 102        format(10x,i4,8x,f7.2)
         end do
c
         if (imode .eq. 0) then
            write(66,876)
 876        format(/,' an estimate of ideal parallel speedup',
     .             ' has been put in file ideal_speedup.dat',/)
         end if
c
      end if
c
      if (nodel .lt. nnodes) then
         write(66,873)
 873     format(/,' NOTE: the value of nodes on the command line was',
     .   ' overspecified',/,' for the number of global-level grids')
         write(11,873)
      end if
c
      rewind(iunit5)
c
c     free up memory used by sizing routine
c
      ifree = 1
      if (ifree.gt.0) then
         call ufree(ip_ltot)
         call ufree(ip_jjmax1)
         call ufree(ip_kkmax1)
         call ufree(ip_iiint1)
         call ufree(ip_iiint2)    
         call ufree(ip_n14)
         call ufree(ip_iwrk)
         call ufree(ip_icsinfo)
         call ufree(ip_iv)
         call ufree(ip_memblock)
         call ufree(ip_lwdat)
         call ufree(ip_nblfine)
         call ufree(ip_mem_req_node)
         call ufree(ip_no_of_points)
         call ufree(ip_lw)
         call ufree(ip_lw2)
         call ufree(ip_nblk)
         call ufree(ip_limblk)
         call ufree(ip_isva)
         call ufree(ip_nblon)
         call ufree(ip_lig)
         call ufree(ip_lbg)
         call ufree(ip_iovrlp)
         call ufree(ip_ibpntsg)
         call ufree(ip_iipntsg)
         call ufree(ip_rkap0g)
         call ufree(ip_levelg)
         call ufree(ip_igridg)
         call ufree(ip_iflimg)
         call ufree(ip_ifdsg)
         call ufree(ip_iviscg)
         call ufree(ip_jdimg)
         call ufree(ip_kdimg)
         call ufree(ip_idimg)
         call ufree(ip_idiagg)
         call ufree(ip_nblcg)
         call ufree(ip_idegg)
         call ufree(ip_jsg)
         call ufree(ip_ksg)
         call ufree(ip_isg)
         call ufree(ip_jeg)
         call ufree(ip_keg)
         call ufree(ip_ieg)
         call ufree(ip_mit)
         call ufree(ip_jlamlog)
         call ufree(ip_klamlog)
         call ufree(ip_ilamlog)
         call ufree(ip_jlamhig)
         call ufree(ip_klamhig)
         call ufree(ip_ilamhig)
         call ufree(ip_iwfg)
         call ufree(ip_utrans)
         call ufree(ip_vtrans)
         call ufree(ip_wtrans)
         call ufree(ip_omegax)
         call ufree(ip_omegay)
         call ufree(ip_omegaz)
         call ufree(ip_xorig)
         call ufree(ip_yorig)
         call ufree(ip_zorig)
         call ufree(ip_dxmx)
         call ufree(ip_dymx)
         call ufree(ip_dzmx)
         call ufree(ip_dthxmx)
         call ufree(ip_dthymx)
         call ufree(ip_dthzmx)
         call ufree(ip_thetax)
         call ufree(ip_thetay)
         call ufree(ip_thetaz)
         call ufree(ip_rfreqt)
         call ufree(ip_rfreqr)
         call ufree(ip_xorig0)
         call ufree(ip_yorig0)
         call ufree(ip_zorig0)
         call ufree(ip_time2)
         call ufree(ip_thetaxl)
         call ufree(ip_thetayl)
         call ufree(ip_thetazl)
         call ufree(ip_itrans)
         call ufree(ip_irotat)
         call ufree(ip_idefrm)
         call ufree(ip_bcvali)
         call ufree(ip_bcvalj)
         call ufree(ip_bcvalk)
         call ufree(ip_nbci0)
         call ufree(ip_nbcj0)
         call ufree(ip_nbck0)
         call ufree(ip_nbcidim)
         call ufree(ip_nbcjdim)
         call ufree(ip_nbckdim)
         call ufree(ip_ibcinfo)
         call ufree(ip_jbcinfo)
         call ufree(ip_kbcinfo)
         call ufree(ip_bcfilei)
         call ufree(ip_bcfilej)
         call ufree(ip_bcfilek)
         call ufree(ip_ncgg)
         call ufree(ip_nblg)
         call ufree(ip_iemg)
         call ufree(ip_inewgg)
         call ufree(ip_inpl3d)
         call ufree(ip_inpr)
         call ufree(ip_iadvance)
         call ufree(ip_iforce)
         call ufree(ip_iindex)
         call ufree(ip_iindx)
         call ufree(ip_llimit)
         call ufree(ip_iitmax)
         call ufree(ip_mmcxie)
         call ufree(ip_mmceta)
         call ufree(ip_ncheck)
         call ufree(ip_iifit)
         call ufree(ip_mblkpt)
         call ufree(ip_iic0)
         call ufree(ip_iiorph)
         call ufree(ip_iitoss)
         call ufree(ip_ifiner)
         call ufree(ip_dx)
         call ufree(ip_dy)
         call ufree(ip_dz)
         call ufree(ip_dthetx)
         call ufree(ip_dthety)
         call ufree(ip_dthetz)
         call ufree(ip_isav_blk)
         call ufree(ip_isav_prd)
         call ufree(ip_isav_pat)
         call ufree(ip_isav_pat_b)
         call ufree(ip_isav_dpat)
         call ufree(ip_isav_dpat_b)
         call ufree(ip_isav_emb)
         call ufree(ip_mblk2nd)
         call ufree(ip_mglevg)
         call ufree(ip_nemgl)
         call ufree(ip_ipl3dtmp)
         call ufree(ip_ifrom)
         call ufree(ip_xif1)
         call ufree(ip_etf1)
         call ufree(ip_xif2)
         call ufree(ip_etf2)
         call ufree(ip_utrnsae)
         call ufree(ip_vtrnsae)
         call ufree(ip_wtrnsae)
         call ufree(ip_omgxae)
         call ufree(ip_omgyae)
         call ufree(ip_omgzae)
         call ufree(ip_xorgae)
         call ufree(ip_yorgae)
         call ufree(ip_zorgae)
         call ufree(ip_xorgae0)
         call ufree(ip_yorgae0)
         call ufree(ip_zorgae0)
         call ufree(ip_icouple)
         call ufree(ip_thtxae)
         call ufree(ip_thtyae)
         call ufree(ip_thtzae)
         call ufree(ip_rfrqtae)
         call ufree(ip_rfrqrae)
         call ufree(ip_icsi)
         call ufree(ip_icsf)
         call ufree(ip_jcsi)
         call ufree(ip_jcsf)
         call ufree(ip_kcsi)
         call ufree(ip_kcsf)
         call ufree(ip_idfrmseg)
         call ufree(ip_iaesurf)
         call ufree(ip_nsegdfrm)
         call ufree(ip_freq)
         call ufree(ip_gmass)
         call ufree(ip_x0)
         call ufree(ip_gf0)
         call ufree(ip_damp)
         call ufree(ip_perturb)
         call ufree(ip_aesrfdat)
         call ufree(ip_iskip)
         call ufree(ip_jskip)
         call ufree(ip_kskip)
      end if
c
      if (imode .eq. 1) then
         write(66,'(/,'' memory for precfl3d has been deallocated'')')
      else
         write(6,'(/,''precfl3d has completed successfully'')')
         write(6,'(''precfl3d information has been put in'',
     .  '' file precfl3d.out'',/)')
      end if
c
      close(66)
c
      return
      end
      subroutine cntfa(nbl,nwfa,maxbl,maxseg,nblcg,ieg,isg,jdimg,kdimg,
     .                 idimg,nbcj0,nbck0,nbci0,nbcjdim,nbckdim,nbcidim,
     .                 jbcinfo,kbcinfo,ibcinfo,nblock)
c***********************************************************************
c     Purpose: count number of flux accumulations needed
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension jdimg(maxbl),kdimg(maxbl),idimg(maxbl),nblcg(maxbl),
     .          isg(maxbl),ieg(maxbl),nbci0(maxbl),nbcidim(maxbl),
     .          nbcj0(maxbl),nbcjdim(maxbl),nbck0(maxbl),nbckdim(maxbl),
     .          ibcinfo(maxbl,maxseg,7,2),jbcinfo(maxbl,maxseg,7,2),
     .          kbcinfo(maxbl,maxseg,7,2)
c
      nwfa = 0
c
      do 6500 nblc=1,nblock
      if (nbl.eq.nblc) go to 6500
c
      nblcc  = nblcg(nblc)
      if (nblcc.eq.nbl) then
         jfm1 = jdimg(nblc)-1
         kfm1 = kdimg(nblc)-1
         ifm1 = idimg(nblc)-1
         nsi  = ifm1/(ieg(nblc)-isg(nblc))
c
         do 802 nseg=1,nbcj0(nblc)
         if (jbcinfo(nblc,nseg,1,1).eq.21) then
            do 100 l=1,5
            do 100 i=1,ifm1,nsi
            do 100 k=1,kfm1,2
            nwfa = nwfa+1
  100       continue
         end if
  802    continue
c
         do 803 nseg=1,nbcjdim(nblc)
         if (jbcinfo(nblc,nseg,1,2).eq.21) then
            do 200 l=1,5
            do 200 i=1,ifm1,nsi
            do 200 k=1,kfm1,2
            nwfa = nwfa+1
  200       continue
         end if
  803    continue
c
         do 804 nseg=1,nbck0(nblc)
         if (kbcinfo(nblc,nseg,1,1).eq.21) then
            do 300 l=1,5
            do 300 i=1,ifm1,nsi
            do 300 j=1,jfm1,2
            nwfa = nwfa+1
  300       continue
         end if
  804    continue
c
         do 805 nseg=1,nbckdim(nblc)
         if (kbcinfo(nblc,nseg,1,2).eq.21) then
            do 400 l=1,5
            do 400 i=1,ifm1,nsi
            do 400 j=1,jfm1,2
            nwfa = nwfa+1
  400       continue
         end if
  805    continue
c
         do 806 nseg=1,nbci0(nblc)
         if (ibcinfo(nblc,nseg,1,1).eq.21) then
            do 500 l=1,5
            do 500 k=1,kfm1,nsi
            do 500 j=1,jfm1,2
            nwfa = nwfa+1
  500       continue
         end if
  806    continue
c
         do 807 nseg=1,nbcidim(nblc)
         if (ibcinfo(nblc,nseg,1,2).eq.21) then
            do 600 l=1,5
            do 600 k=1,kfm1,nsi
            do 600 j=1,jfm1,2
            nwfa = nwfa+1
  600       continue
         end if
  807    continue
      end if
 6500 continue
      return
      end
      subroutine cntblmx(nbl,jdim,kdim,idim,iover,iprint,n14,maxbl,
     .                   maxseg,nbcj0,nbck0,nbci0,nbcjdim,nbckdim,
     .                   nbcidim,jbcinfo,kbcinfo,ibcinfo)
c***********************************************************************
c     Purpose: count the number of lines (n14) output to the baldwin-
c     lomax output file, for use in determining the required size of
c     parameter ibufdim
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension nbci0(maxbl),nbcidim(maxbl),nbcj0(maxbl),
     .          nbcjdim(maxbl),nbck0(maxbl),nbckdim(maxbl),
     .          ibcinfo(maxbl,maxseg,7,2),jbcinfo(maxbl,maxseg,7,2),
     .          kbcinfo(maxbl,maxseg,7,2)
c
      common /reyue/ reue,tinf,ivisc(3)
c
c     WARNING: ilfreq1, ilfreq2, ipfreq1, and ipfreq2 must
c     be defined exactly as in the cfl3d subroutine blomax
c
      ilfreq1 = 1
      ilfreq2 = 5
      ipfreq1 = 10
      ipfreq2 = 10
      if (idim.eq.2 .or. jdim.eq.2 .or. kdim.eq.2) then
         ilfreq2 = 1
         ipfreq2 = 1
      end if
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c 
c     determine which walls to use, based on ibcinfo, jbcinfo, kbcinfo
c
      ibcjmin=0
      ibcjmax=0
      do mseg=1,nbcj0(nbl)
        if (abs(jbcinfo(nbl,mseg,1,1)) .eq. 2004) ibcjmin=1 
      end do
      do mseg=1,nbcjdim(nbl)
        if (abs(jbcinfo(nbl,mseg,1,2)) .eq. 2004) ibcjmax=1 
      end do
c
      ibckmin=0
      ibckmax=0
      do mseg=1,nbck0(nbl)
        if (abs(kbcinfo(nbl,mseg,1,1)) .eq. 2004) ibckmin=1
      end do
      do mseg=1,nbckdim(nbl)
        if (abs(kbcinfo(nbl,mseg,1,2)) .eq. 2004) ibckmax=1
      end do
c
      ibcimin=0
      ibcimax=0
      do mseg=1,nbci0(nbl)
        if (abs(ibcinfo(nbl,mseg,1,1)) .eq. 2004) ibcimin=1
      end do
      do mseg=1,nbcidim(nbl)
        if (abs(ibcinfo(nbl,mseg,1,2)) .eq. 2004) ibcimax=1
      end do
c
      if (ivisc(3) .gt. 1) then
c
         if (ibckmin.eq.1 .or. ibckmax.eq.0) then
c
         ihead = 0
         do 1000 i=1,idim1
         do 1000 j=1,jdim1
c 
         kloop  = .80*kdim
         inmax  = kloop
c
         if (iover.eq.1) then
           if (iprint.ge.1) then
             if (j.eq.j/ilfreq1*ilfreq1.and.i.eq.i/ilfreq2*ilfreq2) then
               if (ihead.eq.0) then
                  n14 = n14 + 1
               end if   
               n14 = n14 + 1
               ihead = ihead+1
             end if
           end if
         end if 
c
         if (iprint.ge.1) then
         if (j.eq.(j/ilfreq1)*ilfreq1.and.i.eq.(i/ilfreq2)*ilfreq2) then
         if (ihead.eq.0) then
            n14 = n14 + 1
         end if
         ihead = ihead+1
         n14 = n14 + 1
         end if
         end if
c
         if (iprint.gt.1) then
         if (j.eq.(j/ipfreq1)*ipfreq1.and.i.eq.(i/ipfreq2)*ipfreq2) then
         ihead = 0
         n14 = n14 + 1
         do 1415 in=1,inmax
         n14 = n14 + 1
 1415    continue
         end if
         end if
 1000    continue
c 
         end if
c
         if (ibckmax.eq.1) then
c
         ihead = 0
         do 10 i=1,idim1
         do 10 j=1,jdim1
c
         kloop  = .80*kdim
         kloop  = kdim1 - kloop + 1
         kloop  = max(1,kloop)
         inmax  = kloop
         inmax1 = inmax + 1
c
         if (iover.eq.1) then
           if (iprint.ge.1) then
             if (j.eq.j/ilfreq1*ilfreq1.and.i.eq.i/ilfreq2*ilfreq2) then
                if (ihead.eq.0) then
                   n14 = n14 + 1
                end if
                n14 = n14 + 1
                ihead = ihead+1
             end if
           end if
         end if
c
         if (iprint.ge.1) then
         if (j.eq.(j/ilfreq1)*ilfreq1.and.i.eq.(i/ilfreq2)*ilfreq2) then
         if (ihead.eq.0) then
            n14 = n14 + 1
         end if
         ihead = ihead+1
         n14 = n14 + 1
         end if
         end if
         if (iprint.gt.1) then
         if (j.eq.(j/ipfreq1)*ipfreq1.and.i.eq.(i/ipfreq2)*ipfreq2) then
         ihead = 0
         n14 = n14 + 1
         do 95 in=kdim1,inmax1,-1
         n14 = n14 + 1
   95    continue
         end if
         end if
c
   10    continue
c
         end if
c
      end if
c
      if (ivisc(2) .gt. 1) then
c
         if (ibcjmin.eq.1 .or. ibcjmax.eq.0) then
c
         ihead = 0
         do 7000 i=1,idim1
         do 7000 k=1,kdim1
         jloop  = .80*jdim
         jloop1 = jloop-1
         inmax  = jloop
c
         if (iover.eq.1) then
           if (iprint.ge.1) then
             if (k.eq.k/ilfreq1*ilfreq1.and.i.eq.i/ilfreq2*ilfreq2) then
                if (ihead.eq.0) then
                   n14 = n14 + 1
                end if
                n14 = n14 + 1
                ihead = ihead+1
              end if
           end if
         end if 
c
         if (iprint.ge.1) then
         if (k.eq.(k/ilfreq1)*ilfreq1.and.i.eq.(i/ilfreq2)*ilfreq2) then
         if (ihead.eq.0) then
            n14 = n14 + 1
         end if
         ihead = ihead+1
         n14 = n14 + 1
         end if
         end if
c
         if (iprint.gt.1) then
         if (k.eq.(k/ipfreq1)*ipfreq1.and.i.eq.(i/ipfreq2)*ipfreq2) then
         ihead = 0
         n14 = n14 + 1
         do 1615 in=1,inmax
         n14 = n14 + 1
 1615    continue
         end if
         end if
c
 7000    continue
c
         end if
c
         if (ibcjmax.eq.1) then
c
         ihead = 0
         do 7610 i=1,idim1
         do 7610 k=1,kdim1
         jloop  = .80*jdim
         jloop  = jdim1 - jloop + 1
         jloop  = max(1,jloop)
         inmax  = jloop
         inmax1 = inmax + 1
c
         if (iover.eq.1) then
           if (iprint.ge.1) then
             if (k.eq.k/ilfreq1*ilfreq1.and.i.eq.i/ilfreq2*ilfreq2) then
                if (ihead.eq.0) then
                   n14 = n14 + 1
                end if
                n14 = n14 + 1
                ihead = ihead+1
             end if
           end if
         end if
c
         if (iprint.ge.1) then
         if (k.eq.(k/ilfreq1)*ilfreq1.and.i.eq.(i/ilfreq2)*ilfreq2) then
         if (ihead.eq.0) then
            n14 = n14 + 1
         end if
         ihead = ihead+1
         n14 = n14 + 1
         end if
         end if
c
         if (iprint.gt.1) then
         if (k.eq.(k/ipfreq1)*ipfreq1.and.i.eq.(i/ipfreq2)*ipfreq2) then
         ihead = 0
         n14 = n14 + 1
         do 495 in=jdim1,inmax1,-1
         n14 = n14 + 1
  495    continue
         end if
         end if
 7610    continue
c
         end if
c
      end if
c
      if (ivisc(1) .gt. 1) then
c
         if (ibcimin.eq.1 .or. ibcimax.eq.0) then
c
         ihead = 0
         do 2000 k=1,kdim1
         do 2000 j=1,jdim1
c 
         iloop  = .80*idim
         inmax  = iloop
c
         if (iover.eq.1) then
           if (iprint.ge.1) then
             if (j.eq.j/ilfreq1*ilfreq1.and.k.eq.k/ilfreq2*ilfreq2) then
               if (ihead.eq.0) then
                  n14 = n14 + 1
               end if
               n14 = n14 + 1
               ihead = ihead+1
             end if
           end if
         end if 
c
         if (iprint.ge.1) then
         if (j.eq.(j/ilfreq1)*ilfreq1.and.k.eq.(k/ilfreq2)*ilfreq2) then
         if (ihead.eq.0) then
            n14 = n14 + 1
         end if
         ihead = ihead+1
         n14 = n14 + 1
         end if
         end if
c
         if (iprint.gt.1) then
         if (j.eq.(j/ipfreq1)*ipfreq1.and.k.eq.(k/ipfreq2)*ipfreq2) then
            ihead = 0
         n14 = n14 + 1
         do 1515 in=1,inmax
         n14 = n14 + 1
 1515    continue
         end if
         end if
c
 2000    continue
c
         end if
c
         if (ibcimax .eq. 1) then
c
         ihead = 0
         do 2007 k=1,kdim1
         do 2007 j=1,jdim1
c 
         iloop  = .80*idim
         iloop  = idim1-iloop+1
         iloop  = max(1,iloop)
         inmax  = iloop
         inmax1 = inmax+1
c
         if (iover.eq.1) then
           if (iprint.ge.1) then
             if (j.eq.j/ilfreq1*ilfreq1.and.k.eq.k/ilfreq2*ilfreq2) then
               if (ihead.eq.0) then
                  n14 = n14 + 1
               end if
               n14 = n14 + 1
               ihead = ihead+1
	     end if
           end if
         end if 
c
         if (iprint.ge.1) then
         if (j.eq.(j/ilfreq1)*ilfreq1.and.k.eq.(k/ilfreq2)*ilfreq2) then
         if (ihead.eq.0) then
            n14 = n14 + 1
         end if
         ihead = ihead+1
         n14 = n14 + 1
         end if
         end if
c
         if (iprint.gt.1) then
         if (j.eq.(j/ipfreq1)*ipfreq1.and.k.eq.(k/ipfreq2)*ipfreq2) then
         ihead = 0
         n14 = n14 + 1
         do 6515 in=idim1,inmax1,-1
         n14 = n14 + 1
 6515    continue
         end if
         end if
 2007    continue
c
         end if
c
      end if
c
      return
      end 
