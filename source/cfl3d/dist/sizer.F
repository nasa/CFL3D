c  ---------------------------------------------------------------------------
c  CFL3D is a structured-grid, cell-centered, upwind-biased, Reynolds-averaged
c  Navier-Stokes (RANS) code. It can be run in parallel on multiple grid zones
c  with point-matched, patched, overset, or embedded connectivities. Both
c  multigrid and mesh sequencing are available in time-accurate or
c  steady-state modes.
c
c  Copyright 2001 United States Government as represented by the Administrator
c  of the National Aeronautics and Space Administration. All Rights Reserved.
c 
c  The CFL3D platform is licensed under the Apache License, Version 2.0 
c  (the "License"); you may not use this file except in compliance with the 
c  License. You may obtain a copy of the License at 
c  http://www.apache.org/licenses/LICENSE-2.0. 
c 
c  Unless required by applicable law or agreed to in writing, software 
c  distributed under the License is distributed on an "AS IS" BASIS, WITHOUT 
c  WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the 
c  License for the specific language governing permissions and limitations 
c  under the License.
c  ---------------------------------------------------------------------------
c
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
c     processors, based on the maximum requirement on any of the nodes.
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
     .             aeinp,aeout,sdhist,avgg,avgq
      character*80 bcfiles(mxbcfil0)
      character*120 bou(ibufdim0,nbuf0)
c
      integer xif1,etf1,xif2,etf2
      integer bcfilei,bcfilej,bcfilek
      integer stats
c
      dimension nou(nbuf0)
c
      allocatable :: aesrfdat(:,:)
      allocatable :: bcfilei(:,:,:)
      allocatable :: bcfilej(:,:,:)
      allocatable :: bcfilek(:,:,:)
      allocatable :: bcvali(:,:,:,:)
      allocatable :: bcvalj(:,:,:,:)
      allocatable :: bcvalk(:,:,:,:)
      allocatable :: damp(:,:)
      allocatable :: dthetx(:,:)
      allocatable :: dthety(:,:)
      allocatable :: dthetz(:,:)
      allocatable :: dthxmx(:)
      allocatable :: dthymx(:)
      allocatable :: dthzmx(:)
      allocatable :: dx(:,:)
      allocatable :: dxmx(:)
      allocatable :: dy(:,:)
      allocatable :: dymx(:)
      allocatable :: dz(:,:)
      allocatable :: dzmx(:)
      allocatable :: etf1(:)
      allocatable :: etf2(:)
      allocatable :: freq(:,:)
      allocatable :: gf0(:,:)
      allocatable :: gmass(:,:)
      allocatable :: iadvance(:)
      allocatable :: iaesurf(:,:)
      allocatable :: ibcinfo(:,:,:,:)
      allocatable :: ibpntsg(:,:)
      allocatable :: icouple(:,:)
      allocatable :: icsf(:,:)
      allocatable :: icsi(:,:)
      allocatable :: icsinfo(:,:)
      allocatable :: idefrm(:)
      allocatable :: idegg(:,:)
      allocatable :: idfrmseg(:,:)
      allocatable :: idiagg(:,:)
      allocatable :: idimg(:)
      allocatable :: ieg(:)
      allocatable :: iemg(:)
      allocatable :: ifdsg(:,:)
      allocatable :: ifiner(:)
      allocatable :: iflimg(:,:)
      allocatable :: iforce(:)
      allocatable :: ifrom(:)
      allocatable :: igridg(:)
      allocatable :: iic0(:)
      allocatable :: iifit(:)
      allocatable :: iiint1(:)
      allocatable :: iiint2(:)
      allocatable :: iindex(:,:)
      allocatable :: iindx(:,:)
      allocatable :: iiorph(:)
      allocatable :: iipntsg(:)
      allocatable :: iitmax(:)
      allocatable :: iitoss(:)
      allocatable :: ilamhig(:)
      allocatable :: ilamlog(:)
      allocatable :: inewgg(:)
      allocatable :: inpl3d(:,:)
      allocatable :: inpr(:,:)
      allocatable :: iovrlp(:)
      allocatable :: ipl3dtmp(:)
      allocatable :: irotat(:)
      allocatable :: isav_blk(:,:)
      allocatable :: isav_dpat(:,:)
      allocatable :: isav_dpat_b(:,:,:)
      allocatable :: isav_emb(:,:)
      allocatable :: isav_pat(:,:)
      allocatable :: isav_pat_b(:,:,:)
      allocatable :: isav_prd(:,:)
      allocatable :: isg(:)
      allocatable :: iskip(:,:)
      allocatable :: isva(:,:,:)
      allocatable :: itrans(:)
      allocatable :: iv(:)
      allocatable :: iviscg(:,:)
      allocatable :: iwfg(:,:)
      allocatable :: iwrk(:,:)
      allocatable :: jbcinfo(:,:,:,:)
      allocatable :: jcsf(:,:)
      allocatable :: jcsi(:,:)
      allocatable :: jdimg(:)
      allocatable :: jeg(:)
      allocatable :: jjmax1(:)
      allocatable :: jlamhig(:)
      allocatable :: jlamlog(:)
      allocatable :: jsg(:)
      allocatable :: jskip(:,:)
      allocatable :: kbcinfo(:,:,:,:)
      allocatable :: kcsf(:,:)
      allocatable :: kcsi(:,:)
      allocatable :: kdimg(:)
      allocatable :: keg(:)
      allocatable :: kkmax1(:)
      allocatable :: klamhig(:)
      allocatable :: klamlog(:)
      allocatable :: ksg(:)
      allocatable :: kskip(:,:)
      allocatable :: lbg(:)
      allocatable :: levelg(:)
      allocatable :: lig(:)
      allocatable :: limblk(:,:,:)
      allocatable :: llimit(:)
      allocatable :: ltot(:)
      allocatable :: lw(:,:)
      allocatable :: lw2(:,:)
      allocatable :: lwdat(:,:,:)
      allocatable :: mblk2nd(:)
      allocatable :: mblkpt(:)
      allocatable :: mem_req_node(:)
      allocatable :: memblock(:)
      allocatable :: mglevg(:)
      allocatable :: mit(:,:)
      allocatable :: mmceta(:)
      allocatable :: mmcxie(:)
      allocatable :: n14(:)
      allocatable :: nbci0(:)
      allocatable :: nbcidim(:)
      allocatable :: nbcj0(:)
      allocatable :: nbcjdim(:)
      allocatable :: nbck0(:)
      allocatable :: nbckdim(:)
      allocatable :: nblcg(:)
      allocatable :: nblfine(:)
      allocatable :: nblg(:)
      allocatable :: nblk(:,:)
      allocatable :: nblon(:)
      allocatable :: ncgg(:)
      allocatable :: ncheck(:)
      allocatable :: nemgl(:)
      allocatable :: no_of_points(:)
      allocatable :: nsegdfrm(:)
      allocatable :: omegax(:)
      allocatable :: omegay(:)
      allocatable :: omegaz(:)
      allocatable :: omgxae(:,:)
      allocatable :: omgyae(:,:)
      allocatable :: omgzae(:,:)
      allocatable :: perturb(:,:,:)
      allocatable :: rfreqr(:)
      allocatable :: rfreqt(:)
      allocatable :: rfrqrae(:,:)
      allocatable :: rfrqtae(:,:)
      allocatable :: rkap0g(:,:)
      allocatable :: thetax(:)
      allocatable :: thetaxl(:)
      allocatable :: thetay(:)
      allocatable :: thetayl(:)
      allocatable :: thetaz(:)
      allocatable :: thetazl(:)
      allocatable :: thtxae(:,:)
      allocatable :: thtyae(:,:)
      allocatable :: thtzae(:,:)
      allocatable :: time2(:)
      allocatable :: utrans(:)
      allocatable :: utrnsae(:,:)
      allocatable :: vtrans(:)
      allocatable :: vtrnsae(:,:)
      allocatable :: wtrans(:)
      allocatable :: wtrnsae(:,:)
      allocatable :: x0(:,:)
      allocatable :: xif1(:)
      allocatable :: xif2(:)
      allocatable :: xorgae(:,:)
      allocatable :: xorgae0(:,:)
      allocatable :: xorig(:)
      allocatable :: xorig0(:)
      allocatable :: yorgae(:,:)
      allocatable :: yorgae0(:,:)
      allocatable :: yorig(:)
      allocatable :: yorig0(:)
      allocatable :: zorgae(:,:)
      allocatable :: zorgae0(:,:)
      allocatable :: zorig(:)
      allocatable :: zorig0(:)
c
      common /des/ cdes,ides,cddes
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
     .        lbcdj,lbcdk,lbcdi,lxib2,lux,lcmuv,lvolj0,lvolk0,lvoli0,
     .        lxmdj,lxmdk,lxmdi,lvelg,ldeltj,ldeltk,ldelti,
     .        lxnm2,lynm2,lznm2,lxnm1,lynm1,lznm1,lqavg
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /maxiv/ ivmx
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /reyue/ reue,tinf,ivisc(3)
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /twod/ i2d
      common /is_blockbc/ is_blk(5),ie_blk(5),ivolint
      common /is_perbc/ is_prd(5),ie_prd(5),nbcprd
      common /is_patch/ is_pat(5),ie_pat(5),ipatch1st
      common /is_embedbc/ is_emb(5),ie_emb(5),nbcemb
      common /filenam/ grid,plt3dg,plt3dq,output,residual,turbres,blomx,
     .                 output2,printout,pplunge,ovrlap,patch,restrt,
     .                 subres,subtur,grdmov,alphahist,errfile,preout,
     .                 aeinp,aeout,sdhist,avgg,avgq
      common /conversion/ radtodeg
      common /unit5/ iunit5
      common /moov/movie,nframes,icall1,lhdr,icoarsemovie,i2dmovie
      common /elastic_ss/ idef_ss
      common /memory/ memadd,memaddi
      common /turbconv/ cflturb(7),edvislim,iturbprod,nsubturb,nfreeze,
     .                  iwarneddy,itime2read,itaturb,tur1cut,tur2cut,
     .                  iturbord,tur1cutlev,tur2cutlev
      common /deformz/ beta1,beta2,alpha1,alpha2,isktyp,negvol,meshdef,
     .                 nsprgit,ndgrd,ndwrt 
      common /curvat/ isarc2d,sarccr3,ieasmcc2d,isstrc,sstrc_crc,
     .        isar,crot,isarc3d
      common /avgdata/ xnumavg,iteravg,xnumavg2,ipertavg,iclcd,isubit_r
c
      icall  = 0
      memuse = 0
c
      allocate( aesrfdat(5,maxaes0), stat=stats )
      call umalloc(5*maxaes0,0,'aesrfdat',memuse,stats)
      allocate( bcfilei(maxbl0,maxseg0,2), stat=stats )
      call umalloc(maxbl0*maxseg0*2,1,'bcfilei',memuse,stats)
      allocate( bcfilej(maxbl0,maxseg0,2), stat=stats )
      call umalloc(maxbl0*maxseg0*2,1,'bcfilej',memuse,stats)
      allocate( bcfilek(maxbl0,maxseg0,2), stat=stats )
      call umalloc(maxbl0*maxseg0*2,1,'bcfilek',memuse,stats)
      allocate( bcvali(maxbl0,maxseg0,12,2), stat=stats )
      call umalloc(maxbl0*maxseg0*12*2,0,'bcvali',memuse,stats)
      allocate( bcvalj(maxbl0,maxseg0,12,2), stat=stats )
      call umalloc(maxbl0*maxseg0*12*2,0,'bcvalj',memuse,stats)
      allocate( bcvalk(maxbl0,maxseg0,12,2), stat=stats )
      call umalloc(maxbl0*maxseg0*12*2,0,'bcvalk',memuse,stats)
      allocate( damp(nmds0,maxaes0), stat=stats )
      call umalloc(nmds0*maxaes0,0,'damp',memuse,stats)
      allocate( dthetx(intmx0,msub10), stat=stats )
      call umalloc(intmx0*msub10,0,'dthetx',memuse,stats)
      allocate( dthety(intmx0,msub10), stat=stats )
      call umalloc(intmx0*msub10,0,'dthety',memuse,stats)
      allocate( dthetz(intmx0,msub10), stat=stats )
      call umalloc(intmx0*msub10,0,'dthetz',memuse,stats)
      allocate( dthxmx(maxbl0), stat=stats )
      call umalloc(maxbl0,0,'dthxmx',memuse,stats)
      allocate( dthymx(maxbl0), stat=stats )
      call umalloc(maxbl0,0,'dthymx',memuse,stats)
      allocate( dthzmx(maxbl0), stat=stats )
      call umalloc(maxbl0,0,'dthzmx',memuse,stats)
      allocate( dx(intmx0,msub10), stat=stats )
      call umalloc(intmx0*msub10,0,'dx',memuse,stats)
      allocate( dxmx(maxbl0), stat=stats )
      call umalloc(maxbl0,0,'dxmx',memuse,stats)
      allocate( dy(intmx0,msub10), stat=stats )
      call umalloc(intmx0*msub10,0,'dy',memuse,stats)
      allocate( dymx(maxbl0), stat=stats )
      call umalloc(maxbl0,0,'dymx',memuse,stats)
      allocate( dz(intmx0,msub10), stat=stats )
      call umalloc(intmx0*msub10,0,'dz',memuse,stats)
      allocate( dzmx(maxbl0), stat=stats )
      call umalloc(maxbl0,0,'dzmx',memuse,stats)
      allocate( etf1(msub10), stat=stats )
      call umalloc(msub10,1,'etf1',memuse,stats)
      allocate( etf2(msub10), stat=stats )
      call umalloc(msub10,1,'etf2',memuse,stats)
      allocate( freq(nmds0,maxaes0), stat=stats )
      call umalloc(nmds0*maxaes0,0,'freq',memuse,stats)
      allocate( gf0(2*nmds0,maxaes0), stat=stats )
      call umalloc(2*nmds0*maxaes0,0,'gf0',memuse,stats)
      allocate( gmass(nmds0,maxaes0), stat=stats )
      call umalloc(nmds0*maxaes0,0,'gmass',memuse,stats)
      allocate( iadvance(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'iadvance',memuse,stats)
      allocate( iaesurf(maxbl0,maxsegdg0), stat=stats )
      call umalloc(maxbl0*maxsegdg0,1,'iaesurf',memuse,stats)
      allocate( ibcinfo(maxbl0,maxseg0,7,2), stat=stats )
      call umalloc(maxbl0*maxseg0*7*2,1,'ibcinfo',memuse,stats)
      allocate( ibpntsg(maxbl0,4), stat=stats )
      call umalloc(maxbl0*4,1,'ibpntsg',memuse,stats)
      allocate( icouple(maxbl0,maxsegdg0), stat=stats )
      call umalloc(maxbl0*maxsegdg0,1,'icouple',memuse,stats)
      allocate( icsf(maxbl0,maxsegdg0), stat=stats )
      call umalloc(maxbl0*maxsegdg0,1,'icsf',memuse,stats)
      allocate( icsi(maxbl0,maxsegdg0), stat=stats )
      call umalloc(maxbl0*maxsegdg0,1,'icsi',memuse,stats)
      allocate( icsinfo(maxcs0,9), stat=stats )
      call umalloc(maxcs0*9,1,'icsinfo',memuse,stats)
      allocate( idefrm(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'idefrm',memuse,stats)
      allocate( idegg(maxbl0,3), stat=stats )
      call umalloc(maxbl0*3,1,'idegg',memuse,stats)
      allocate( idfrmseg(maxbl0,maxsegdg0), stat=stats )
      call umalloc(maxbl0*maxsegdg0,1,'idfrmseg',memuse,stats)
      allocate( idiagg(maxbl0,3), stat=stats )
      call umalloc(maxbl0*3,1,'idiagg',memuse,stats)
      allocate( idimg(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'idimg',memuse,stats)
      allocate( ieg(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'ieg',memuse,stats)
      allocate( iemg(maxgr0), stat=stats )
      call umalloc(maxgr0,1,'iemg',memuse,stats)
      allocate( ifdsg(maxbl0,3), stat=stats )
      call umalloc(maxbl0*3,1,'ifdsg',memuse,stats)
      allocate( ifiner(intmx0), stat=stats )
      call umalloc(intmx0,1,'ifiner',memuse,stats)
      allocate( iflimg(maxbl0,3), stat=stats )
      call umalloc(maxbl0*3,1,'iflimg',memuse,stats)
      allocate( iforce(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'iforce',memuse,stats)
      allocate( ifrom(msub10), stat=stats )
      call umalloc(msub10,1,'ifrom',memuse,stats)
      allocate( igridg(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'igridg',memuse,stats)
      allocate( iic0(intmx0), stat=stats )
      call umalloc(intmx0,1,'iic0',memuse,stats)
      allocate( iifit(intmx0), stat=stats )
      call umalloc(intmx0,1,'iifit',memuse,stats)
      allocate( iiint1(nsub10), stat=stats )
      call umalloc(nsub10,1,'iiint1',memuse,stats)
      allocate( iiint2(nsub10), stat=stats )
      call umalloc(nsub10,1,'iiint2',memuse,stats)
      allocate( iindex(intmax0,6*nsub10+9), stat=stats )
      call umalloc(intmax0*(6*nsub10+9),1,'iindex',memuse,stats)
      allocate( iindx(intmx0,6*msub10+9), stat=stats )
      call umalloc(intmx0*(6*msub10+9),1,'iindx',memuse,stats)
      allocate( iiorph(intmx0), stat=stats )
      call umalloc(intmx0,1,'iiorph',memuse,stats)
      allocate( iipntsg(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'iipntsg',memuse,stats)
      allocate( iitmax(intmx0), stat=stats )
      call umalloc(intmx0,1,'iitmax',memuse,stats)
      allocate( iitoss(intmx0), stat=stats )
      call umalloc(intmx0,1,'iitoss',memuse,stats)
      allocate( ilamhig(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'ilamhig',memuse,stats)
      allocate( ilamlog(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'ilamlog',memuse,stats)
      allocate( inewgg(maxgr0), stat=stats )
      call umalloc(maxgr0,1,'inewgg',memuse,stats)
      allocate( inpl3d(nplots0,11), stat=stats )
      call umalloc(nplots0*11,1,'inpl3d',memuse,stats)
      allocate( inpr(nplots0,11), stat=stats )
      call umalloc(nplots0*11,1,'inpr',memuse,stats)
      allocate( iovrlp(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'iovrlp',memuse,stats)
      allocate( ipl3dtmp(11*nplots0), stat=stats )
      call umalloc(11*nplots0,1,'ipl3dtmp',memuse,stats)
      allocate( irotat(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'irotat',memuse,stats)
      allocate( isav_blk(2*mxbli0,17), stat=stats )
      call umalloc(2*mxbli0*17,1,'isav_blk',memuse,stats)
      allocate( isav_dpat(intmx0,17), stat=stats )
      call umalloc(intmx0*17,1,'isav_dpat',memuse,stats)
      allocate( isav_dpat_b(intmx0,msub10,6), stat=stats )
      call umalloc(intmx0*msub10*6,1,'isav_dpat_b',memuse,stats)
      allocate( isav_emb(lbcemb0,12), stat=stats )
      call umalloc(lbcemb0*12,1,'isav_emb',memuse,stats)
      allocate( isav_pat(intmax0,17), stat=stats )
      call umalloc(intmax0*17,1,'isav_pat',memuse,stats)
      allocate( isav_pat_b(intmax0,nsub10,6), stat=stats )
      call umalloc(intmax0*nsub10*6,1,'isav_pat_b',memuse,stats)
      allocate( isav_prd(lbcprd0,12), stat=stats )
      call umalloc(lbcprd0*12,1,'isav_prd',memuse,stats)
      allocate( isg(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'isg',memuse,stats)
      allocate( iskip(maxbl0,500), stat=stats )
      call umalloc(500*maxbl0,1,'iskip',memuse,stats)
      allocate( isva(2,2,mxbli0), stat=stats )
      call umalloc(2*2*mxbli0,1,'isva',memuse,stats)
      allocate( itrans(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'itrans',memuse,stats)
      allocate( iv(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'iv',memuse,stats)
      allocate( iviscg(maxbl0,3), stat=stats )
      call umalloc(maxbl0*3,1,'iviscg',memuse,stats)
      allocate( iwfg(maxbl0,3), stat=stats )
      call umalloc(maxbl0*3,1,'iwfg',memuse,stats)
      allocate( iwrk(maxbl0,3), stat=stats )
      call umalloc(maxbl0*3,1,'iwrk',memuse,stats)
      allocate( jbcinfo(maxbl0,maxseg0,7,2), stat=stats )
      call umalloc(maxbl0*maxseg0*7*2,1,'jbcinfo',memuse,stats)
      allocate( jcsf(maxbl0,maxsegdg0), stat=stats )
      call umalloc(maxbl0*maxsegdg0,1,'jcsf',memuse,stats)
      allocate( jcsi(maxbl0,maxsegdg0), stat=stats )
      call umalloc(maxbl0*maxsegdg0,1,'jcsi',memuse,stats)
      allocate( jdimg(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'jdimg',memuse,stats)
      allocate( jeg(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'jeg',memuse,stats)
      allocate( jjmax1(nsub10), stat=stats )
      call umalloc(nsub10,1,'jjmax1',memuse,stats)
      allocate( jlamhig(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'jlamhig',memuse,stats)
      allocate( jlamlog(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'jlamlog',memuse,stats)
      allocate( jsg(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'jsg',memuse,stats)
      allocate( jskip(maxbl0,500), stat=stats )
      call umalloc(500*maxbl0,1,'jskip',memuse,stats)
      allocate( kbcinfo(maxbl0,maxseg0,7,2), stat=stats )
      call umalloc(maxbl0*maxseg0*7*2,1,'kbcinfo',memuse,stats)
      allocate( kcsf(maxbl0,maxsegdg0), stat=stats )
      call umalloc(maxbl0*maxsegdg0,1,'kcsf',memuse,stats)
      allocate( kcsi(maxbl0,maxsegdg0), stat=stats )
      call umalloc(maxbl0*maxsegdg0,1,'kcsi',memuse,stats)
      allocate( kdimg(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'kdimg',memuse,stats)
      allocate( keg(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'keg',memuse,stats)
      allocate( kkmax1(nsub10), stat=stats )
      call umalloc(nsub10,1,'kkmax1',memuse,stats)
      allocate( klamhig(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'klamhig',memuse,stats)
      allocate( klamlog(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'klamlog',memuse,stats)
      allocate( ksg(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'ksg',memuse,stats)
      allocate( kskip(maxbl0,500), stat=stats )
      call umalloc(500*maxbl0,1,'kskip',memuse,stats)
      allocate( lbg(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'lbg',memuse,stats)
      allocate( levelg(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'levelg',memuse,stats)
      allocate( lig(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'lig',memuse,stats)
      allocate( limblk(2,6,mxbli0), stat=stats )
      call umalloc(2*6*mxbli0,1,'limblk',memuse,stats)
      allocate( llimit(intmx0), stat=stats )
      call umalloc(intmx0,1,'llimit',memuse,stats)
      allocate( ltot(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'ltot',memuse,stats)
      allocate( lw(65,maxbl0), stat=stats )
      call umalloc(65*maxbl0,1,'lw',memuse,stats)
      allocate( lw2(43,maxbl0), stat=stats )
      call umalloc(43*maxbl0,1,'lw2',memuse,stats)
      allocate( lwdat(maxbl0,maxseg0,6), stat=stats )
      call umalloc(maxbl0*maxseg0*6,1,'lwdat',memuse,stats)
      allocate( mblk2nd(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'mblk2nd',memuse,stats)
      allocate( mblkpt(mxxe0), stat=stats )
      call umalloc(mxxe0,1,'mblkpt',memuse,stats)
      allocate( mem_req_node(maxnode0), stat=stats )
      call umalloc(maxnode0,1,'mem_req_node',memuse,stats)
      allocate( memblock(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'memblock',memuse,stats)
      allocate( mglevg(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'mglevg',memuse,stats)
      allocate( mit(5,maxbl0), stat=stats )
      call umalloc(5*maxbl0,1,'mit',memuse,stats)
      allocate( mmceta(intmx0), stat=stats )
      call umalloc(intmx0,1,'mmceta',memuse,stats)
      allocate( mmcxie(intmx0), stat=stats )
      call umalloc(intmx0,1,'mmcxie',memuse,stats)
      allocate( n14(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'n14',memuse,stats)
      allocate( nbci0(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'nbci0',memuse,stats)
      allocate( nbcidim(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'nbcidim',memuse,stats)
      allocate( nbcj0(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'nbcj0',memuse,stats)
      allocate( nbcjdim(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'nbcjdim',memuse,stats)
      allocate( nbck0(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'nbck0',memuse,stats)
      allocate( nbckdim(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'nbckdim',memuse,stats)
      allocate( nblcg(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'nblcg',memuse,stats)
      allocate( nblfine(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'nblfine',memuse,stats)
      allocate( nblg(maxgr0), stat=stats )
      call umalloc(maxgr0,1,'nblg',memuse,stats)
      allocate( nblk(2,mxbli0), stat=stats )
      call umalloc(2*mxbli0,1,'nblk',memuse,stats)
      allocate( nblon(mxbli0), stat=stats )
      call umalloc(mxbli0,1,'nblon',memuse,stats)
      allocate( ncgg(maxgr0), stat=stats )
      call umalloc(maxgr0,1,'ncgg',memuse,stats)
      allocate( ncheck(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'ncheck',memuse,stats)
      allocate( nemgl(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'nemgl',memuse,stats)
      allocate( no_of_points(maxnode0), stat=stats )
      call umalloc(maxnode0,1,'no_of_points',memuse,stats)
      allocate( nsegdfrm(maxbl0), stat=stats )
      call umalloc(maxbl0,1,'nsegdfrm',memuse,stats)
      allocate( omegax(maxbl0), stat=stats )
      call umalloc(maxbl0,0,'omegax',memuse,stats)
      allocate( omegay(maxbl0), stat=stats )
      call umalloc(maxbl0,0,'omegay',memuse,stats)
      allocate( omegaz(maxbl0), stat=stats )
      call umalloc(maxbl0,0,'omegaz',memuse,stats)
      allocate( omgxae(maxbl0,maxsegdg0), stat=stats )
      call umalloc(maxbl0*maxsegdg0,0,'omgxae',memuse,stats)
      allocate( omgyae(maxbl0,maxsegdg0), stat=stats )
      call umalloc(maxbl0*maxsegdg0,0,'omgyae',memuse,stats)
      allocate( omgzae(maxbl0,maxsegdg0), stat=stats )
      call umalloc(maxbl0*maxsegdg0,0,'omgzae',memuse,stats)
      allocate( perturb(nmds0,maxaes0,4), stat=stats )
      call umalloc(nmds0*maxaes0*4,0,'perturb',memuse,stats)
      allocate( rfreqr(maxbl0), stat=stats )
      call umalloc(maxbl0,0,'rfreqr',memuse,stats)
      allocate( rfreqt(maxbl0), stat=stats )
      call umalloc(maxbl0,0,'rfreqt',memuse,stats)
      allocate( rfrqrae(maxbl0,maxsegdg0), stat=stats )
      call umalloc(maxbl0*maxsegdg0,0,'rfrqrae',memuse,stats)
      allocate( rfrqtae(maxbl0,maxsegdg0), stat=stats )
      call umalloc(maxbl0*maxsegdg0,0,'rfrqtae',memuse,stats)
      allocate( rkap0g(maxbl0,3), stat=stats )
      call umalloc(maxbl0*3,0,'rkap0g',memuse,stats)
      allocate( thetax(maxbl0), stat=stats )
      call umalloc(maxbl0,0,'thetax',memuse,stats)
      allocate( thetaxl(maxbl0), stat=stats )
      call umalloc(maxbl0,0,'thetaxl',memuse,stats)
      allocate( thetay(maxbl0), stat=stats )
      call umalloc(maxbl0,0,'thetay',memuse,stats)
      allocate( thetayl(maxbl0), stat=stats )
      call umalloc(maxbl0,0,'thetayl',memuse,stats)
      allocate( thetaz(maxbl0), stat=stats )
      call umalloc(maxbl0,0,'thetaz',memuse,stats)
      allocate( thetazl(maxbl0), stat=stats )
      call umalloc(maxbl0,0,'thetazl',memuse,stats)
      allocate( thtxae(maxbl0,maxsegdg0), stat=stats )
      call umalloc(maxbl0*maxsegdg0,0,'thtxae',memuse,stats)
      allocate( thtyae(maxbl0,maxsegdg0), stat=stats )
      call umalloc(maxbl0*maxsegdg0,0,'thtyae',memuse,stats)
      allocate( thtzae(maxbl0,maxsegdg0), stat=stats )
      call umalloc(maxbl0*maxsegdg0,0,'thtzae',memuse,stats)
      allocate( time2(maxbl0), stat=stats )
      call umalloc(maxbl0,0,'time2',memuse,stats)
      allocate( utrans(maxbl0), stat=stats )
      call umalloc(maxbl0,0,'utrans',memuse,stats)
      allocate( utrnsae(maxbl0,maxsegdg0), stat=stats )
      call umalloc(maxbl0*maxsegdg0,0,'utrnsae',memuse,stats)
      allocate( vtrans(maxbl0), stat=stats )
      call umalloc(maxbl0,0,'vtrans',memuse,stats)
      allocate( vtrnsae(maxbl0,maxsegdg0), stat=stats )
      call umalloc(maxbl0*maxsegdg0,0,'vtrnsae',memuse,stats)
      allocate( wtrans(maxbl0), stat=stats )
      call umalloc(maxbl0,0,'wtrans',memuse,stats)
      allocate( wtrnsae(maxbl0,maxsegdg0), stat=stats )
      call umalloc(maxbl0*maxsegdg0,0,'wtrnsae',memuse,stats)
      allocate( x0(2*nmds0,maxaes0), stat=stats )
      call umalloc(2*nmds0*maxaes0,0,'x0',memuse,stats)
      allocate( xif1(msub10), stat=stats )
      call umalloc(msub10,1,'xif1',memuse,stats)
      allocate( xif2(msub10), stat=stats )
      call umalloc(msub10,1,'xif2',memuse,stats)
      allocate( xorgae(maxbl0,maxsegdg0), stat=stats )
      call umalloc(maxbl0*maxsegdg0,0,'xorgae',memuse,stats)
      allocate( xorgae0(maxbl0,maxsegdg0), stat=stats )
      call umalloc(maxbl0*maxsegdg0,0,'xorgae0',memuse,stats)
      allocate( xorig(maxbl0), stat=stats )
      call umalloc(maxbl0,0,'xorig',memuse,stats)
      allocate( xorig0(maxbl0), stat=stats )
      call umalloc(maxbl0,0,'xorig0',memuse,stats)
      allocate( yorgae(maxbl0,maxsegdg0), stat=stats )
      call umalloc(maxbl0*maxsegdg0,0,'yorgae',memuse,stats)
      allocate( yorgae0(maxbl0,maxsegdg0), stat=stats )
      call umalloc(maxbl0*maxsegdg0,0,'yorgae0',memuse,stats)
      allocate( yorig(maxbl0), stat=stats )
      call umalloc(maxbl0,0,'yorig',memuse,stats)
      allocate( yorig0(maxbl0), stat=stats )
      call umalloc(maxbl0,0,'yorig0',memuse,stats)
      allocate( zorgae(maxbl0,maxsegdg0), stat=stats )
      call umalloc(maxbl0*maxsegdg0,0,'zorgae',memuse,stats)
      allocate( zorgae0(maxbl0,maxsegdg0), stat=stats )
      call umalloc(maxbl0*maxsegdg0,0,'zorgae0',memuse,stats)
      allocate( zorig(maxbl0), stat=stats )
      call umalloc(maxbl0,0,'zorig',memuse,stats)
      allocate( zorig0(maxbl0), stat=stats )
      call umalloc(maxbl0,0,'zorig0',memuse,stats)
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
      pi       = 4.*atan(1.0)
      radtodeg = 180.e0/pi
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
 9990 format(2(2h *),43h   VERSION 6.7 :  Computational Fluids Lab,,
     .15h Mail Stop 128,,4x,2(2h *),
     ./2(2h *),18x,41hNASA Langley Research Center, Hampton, VA,
     .3x,2(2h *),/2(2h *),18x,33hRelease Date:  February  1, 2017.,
     .11x,2(2h *))
      write(66,87)
      write(66,83)
      write(66,83)
   83 format(35(2h *))
   87 format(2(2h *),62x,2(2h *))
c
#ifdef CRAY
c     cray always double precision
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
         do ii=1,65
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
            iwk2  = iwk1+jdim*kdim*idim*7
            iwk3  = iwk2+jj2*kk2*ii2*7
            itemp = iwk3 + jdim*kk2*ii2*7
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
         iwk2  = iwk1 + jdim*kdim*idim*7
         iwk3  = iwk2 + jj2*kk2*ii2*7
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
         if (ivmx.ge.2)then
            if (ivmx .eq. 2)then
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
            else if (ivmx .eq. 4) then
c              Baldwin-Barth
               iwk1=1
               if (iturbord .eq. 1) then
                 iwk5=iwk1+(jdim+1)*(kdim+1)*(idim+1-(2*i2d))
               else
                 iwk5=iwk1+(jdim+3)*(kdim+3)*(idim+3-(4*i2d))
               end if
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
            else if (ivmx .eq. 5) then
c              Spalart
               iwk1=1
               if (iturbord .eq. 1) then
                 iwk3=iwk1+(jdim+1)*(kdim+1)*(idim+1-(2*i2d))
               else
                 iwk3=iwk1+(jdim+3)*(kdim+3)*(idim+3-(4*i2d))
               end if
               iwk4=iwk3+(jdim+1)*(kdim+1)*(idim+1-(2*i2d))
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
               iwk19=iwk18+inmx
               if (isarc2d .eq. 1) then
                 inmx=(jdim+1)*(kdim+1)*(idim-1)*4
               else if (isarc3d .eq. 1) then
                 inmx=(jdim+1)*(kdim+1)*(idim-1)*6
               else
                 inmx=0
               end if
               itemp=iwk19+inmx
            else if (ivmx .ge. 6 .and. ivmx .le. 16) then
c              Two-equation turb models
               iwk1=1
               if (iturbord .eq. 1) then
                 iwk4=iwk1+(jdim+1)*(kdim+1)*(idim+1-(2*i2d))*2
               else
                 iwk4=iwk1+(jdim+3)*(kdim+3)*(idim+3-(4*i2d))*2
               end if
               iwk5=iwk4+(jdim+1)*(kdim+1)*(idim+1-(2*i2d))
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
               inmx=(jdim+1)*(kdim+1)*(idim+1-(2*i2d))
c allocate for gradients of vorticity and turbulence - k-zeta only:
               iwk34=iwk33+inmx
               if (ivmx .eq. 15) then
               inmx=jdim
               else
               inmx=0
               end if
               iwk35=iwk34+inmx
               iwk36=iwk35+inmx
               iwk37=iwk36+inmx
               iwk38=iwk37+inmx
               iwk39=iwk38+inmx
               iwk40=iwk39+inmx
               iwk41=iwk40+inmx
               iwk42=iwk41+inmx
               iwk43=iwk42+inmx
               iwk44=iwk43+inmx
               iwk45=iwk44+inmx
               iwk46=iwk45+inmx
               iwk47=iwk46+inmx
               iwk48=iwk47+inmx
               iwk49=iwk48+inmx
               iwk50=iwk49+inmx
               iwk51=iwk50+inmx
               iwk52=iwk51+inmx
               iwk53=iwk52+inmx
               iwk54=iwk53+inmx
               iwk55=iwk54+inmx
               iwk56=iwk55+inmx
               iwk57=iwk56+inmx
               iwk58=iwk57+inmx
               iwk59=iwk58+inmx
               iwk60=iwk59+inmx
               iwk61=iwk60+inmx
               iwk62=iwk61+inmx
               iwk63=iwk62+inmx
               iwk64=iwk63+inmx
               iwk65=iwk64+inmx
               iwk66=iwk65+inmx
               iwk67=iwk66+inmx
               iwk68=iwk67+inmx
               iwk69=iwk68+inmx
               iwk70=iwk69+inmx
               iwk71=iwk70+inmx
               iwk72=iwk71+inmx
               iwk73=iwk72+inmx
               iwk74=iwk73+inmx
               iwk75=iwk74+inmx
               iwk76=iwk75+inmx
               if (ivmx .eq. 15) then
               inmx=(jdim+1)*(kdim+1)*(idim+1)
               else
               inmx=0
               end if
               iwk77=iwk76+inmx
               iwk78=iwk77+inmx
               iwk79=iwk78+inmx
               if (ivmx .eq. 15) then
               inmx=(jdim+1)*(kdim+1)*(idim+1)*2
               else
               inmx=0
               end if
               iwk80=iwk79+inmx
c allocate for easmcc2d
               if (ieasmcc2d .eq. 1) then
                 inmx=(jdim+1)*(kdim+1)*(idim-1)*4
               else if (isstrc .eq. 2) then
                 inmx=(jdim+1)*(kdim+1)*(idim-1)*6
               else if (ivmx .eq. 16) then
                 inmx=(jdim+1)*(kdim+1)*(idim-1)*4
               else
                 inmx=0
               end if
               iwk81=iwk80+inmx
c allocate for DES and DDES
               if (ides .ne. 0) then
                 inmx=(jdim-1)*(kdim-1)*(idim-1)
               else
                 inmx=0
               end if
               iwk82=iwk81+inmx
               if (ides .eq. 3) then
                 inmx=(jdim-1)*(kdim-1)*(idim-1)
               else
                 inmx=0
               end if
               itemp=iwk82+inmx
            else if (ivmx .eq. 30) then
c   three-equation turb model
               iwk1=1
               if (iturbord .eq. 1) then
                 iwk4=iwk1+(jdim+1)*(kdim+1)*(idim+1-(2*i2d))*3
                 iex=0
                 iex2=-i2d
               else
                 iwk4=iwk1+(jdim+3)*(kdim+3)*(idim+3-(4*i2d))*3
                 iex=1
                 iex2=1-(2*i2d)
               end if
               iwk5=iwk4+(jdim+1)*(kdim+1)*(idim+1-(2*i2d))
               iex3=-i2d
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
               inmx=(jdim-1)*(kdim-1)*(idim-1)*3
               iwk33=iwk32+inmx
               inmx=(jdim+1)*(kdim+1)*(idim+1-(2*i2d))
               iwk34=iwk33+inmx
               inmx=(kdim-1)*(jdim-1)
               iwk35=iwk34+inmx
               iwk36=iwk35+inmx
               iwk37=iwk36+inmx
               iwk38=iwk37+inmx
               inmx=(jdim-1)*(kdim-1)
               iwk39=iwk38+inmx
               iwk40=iwk39+inmx
               iwk41=iwk40+inmx
               iwk42=iwk41+inmx
               inmx=(kdim-1)*(idim-1)
               iwk43=iwk42+inmx
               iwk44=iwk43+inmx
               iwk45=iwk44+inmx
               itemp=iwk45+inmx
            else if (ivmx .eq. 40) then
c   four-equation turb model
               iwk1=1
               if (iturbord .eq. 1) then
                 iwk4=iwk1+(jdim+1)*(kdim+1)*(idim+1-(2*i2d))*4
                 iex=0
                 iex2=-i2d
               else
                 iwk4=iwk1+(jdim+3)*(kdim+3)*(idim+3-(4*i2d))*4
                 iex=1
                 iex2=1-(2*i2d)
               end if
               iwk5=iwk4+(jdim+1)*(kdim+1)*(idim+1-(2*i2d))
               iex3=-i2d
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
               inmx=(jdim-1)*(kdim-1)*(idim-1)*4
               iwk33=iwk32+inmx
               inmx=(jdim+1)*(kdim+1)*(idim+1-(2*i2d))
               iwk34=iwk33+inmx
               inmx=(kdim-1)*(jdim-1)
               iwk35=iwk34+inmx
               iwk36=iwk35+inmx
               iwk37=iwk36+inmx
               iwk38=iwk37+inmx
               inmx=(jdim-1)*(kdim-1)
               iwk39=iwk38+inmx
               iwk40=iwk39+inmx
               iwk41=iwk40+inmx
               iwk42=iwk41+inmx
               inmx=(kdim-1)*(idim-1)
               iwk43=iwk42+inmx
               iwk44=iwk43+inmx
               iwk45=iwk44+inmx
               iwk46=iwk45+inmx
               inmx=(kdim-1)*(jdim-1)
               iwk47=iwk46+inmx
               iwk48=iwk47+inmx
               iwk49=iwk48+inmx
               iwk50=iwk49+inmx
               inmx=(jdim-1)*(kdim-1)
               iwk51=iwk50+inmx
               iwk52=iwk51+inmx
               iwk53=iwk52+inmx
               iwk54=iwk53+inmx
               inmx=(kdim-1)*(idim-1)
               iwk55=iwk54+inmx
               iwk56=iwk55+inmx
               iwk57=iwk56+inmx
               iwk58=iwk57+inmx
c allocate for DES and DDES
               if (ides .ne. 0) then
                 inmx=(jdim-1)*(kdim-1)*(idim-1)
               else
                 inmx=0
               end if
               iwk59=iwk58+inmx
               if (ides .eq. 3) then
                 inmx=(jdim-1)*(kdim-1)*(idim-1)
               else
                 inmx=0
               end if
               itemp=iwk59+inmx
            else if(ivmx.eq.72) then 
               itemp = 0
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
         call cntsurf(ns2004,maxbl0,maxgr0,maxseg0,ngrid,nblg,
     .                nbci0,nbcj0,nbck0,nbcidim,nbcjdim,
     .                nbckdim,ibcinfo,jbcinfo,kbcinfo,2004)
         call cntsurf(ns2014,maxbl0,maxgr0,maxseg0,ngrid,nblg,
     .                nbci0,nbcj0,nbck0,nbcidim,nbcjdim,
     .                nbckdim,ibcinfo,jbcinfo,kbcinfo,2014)
         call cntsurf(ns2024,maxbl0,maxgr0,maxseg0,ngrid,nblg,
     .                nbci0,nbcj0,nbck0,nbcidim,nbcjdim,
     .                nbckdim,ibcinfo,jbcinfo,kbcinfo,2024)
         call cntsurf(ns2034,maxbl0,maxgr0,maxseg0,ngrid,nblg,
     .                nbci0,nbcj0,nbck0,nbcidim,nbcjdim,
     .                nbckdim,ibcinfo,jbcinfo,kbcinfo,2034)
         call cntsurf(ns2016,maxbl0,maxgr0,maxseg0,ngrid,nblg,
     .                nbci0,nbcj0,nbck0,nbcidim,nbcjdim,
     .                nbckdim,ibcinfo,jbcinfo,kbcinfo,2016)
         nsurf=ns2004+ns2014+ns2024+ns2034+ns2016
         minbox=sqrt(float(nsurf))
         minbox=max(minbox,50)
         nbb = 3*nsurf/minbox
         do 2010 igrid = 1,ngrid
         nbll = nblg(igrid)
         ntempf =  9*nsurf + 7*nbb
     .          +  4*jdimg(nbll)*kdimg(nbll)*idimg(nbll)
         if (ivmx.eq.4 .or. ivmx.eq.25) then
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
         nworkxs  = 4*nsurf
         if (ivmx.eq.4 .or. ivmx.eq.25) then
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
         maxmem = max(jdim*kdim*nv,jdim*idim*nv)
         if (idefrm(nbl) .gt. 0) then
           nv=41
           maxmem = max(jdim*kdim*nv,jdim*idim*nv)
           mem_updatedg =     16*idim*jdim*kdim +
     .                             30*jdim*kdim +
     .                             24*jdim*idim +
     .                             24*kdim*idim
           maxmem = maxmem+mem_updatedg 
         end if
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
c     requirements for qoutavg (plot3davg) output for iteravg>0
c
      if(iteravg.gt.0)then
      do n=1,nplot3d
         nbl = inpl3d(n,1)
         call lead(nbl,lw,lw2,maxbl0)
         i1 = 1
         i2 = idim
         i3 = 1
         j1 = 1
         j2 = jdim
         j3 = 1
         k1 = 1
         k2 = kdim
         k3 = 1
            jdw = (j2-j1)/j3 + 1
            kdw = (k2-k1)/k3 + 1
            idw = (i2-i1)/i3 + 1
            itemp = jdw*kdw*idw*14 + jdim*kdim*idim*2
         need  = max(need,itemp)
         ii = mblk2nd(nbl)
         mem_req_node(ii) = max(mem_req_node(ii),itemp)
         write(66,*)' 23(plot3davg)   itemp, need = ',itemp,need
         write(66,*)'       parallel: itemp, need = ',itemp,
     .   mem_req_node(ii)
      end do
      end if
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
            itemp = jdw*kdw*idw*17 + jdim*kdim*idim*2
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
     .                  itemp = itemp + maxdims*7*2
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
     .                  itemp = itemp + maxdims*7*2
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
     .                     itemp = itemp + jmax1*kmax1*7*2
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
     .                     itemp = itemp + jmax1*kmax1*7*2
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
     .                  itemp = itemp + maxdims*7*2
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
     .                  itemp = itemp + maxdims*7*2
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
     .                  itemp = itemp + maxdims*7*3
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
     .                  itemp = itemp + maxdims*7*3
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
      needi_trn =1550*lmaxbl       + 119      + 9*lmaxbl*lmxsegdg
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
         call cntsurf(ns2014,maxbl0,maxgr0,maxseg0,ngrid,nblg,
     .                nbci0,nbcj0,nbck0,nbcidim,nbcjdim,
     .                nbckdim,ibcinfo,jbcinfo,kbcinfo,2014)
         call cntsurf(ns2024,maxbl0,maxgr0,maxseg0,ngrid,nblg,
     .                nbci0,nbcj0,nbck0,nbcidim,nbcjdim,
     .                nbckdim,ibcinfo,jbcinfo,kbcinfo,2024)
         call cntsurf(ns2034,maxbl0,maxgr0,maxseg0,ngrid,nblg,
     .                nbci0,nbcj0,nbck0,nbcidim,nbcjdim,
     .                nbckdim,ibcinfo,jbcinfo,kbcinfo,2034)
         call cntsurf(ns1005,maxbl0,maxgr0,maxseg0,ngrid,nblg,
     .                nbci0,nbcj0,nbck0,nbcidim,nbcjdim,
     .                nbckdim,ibcinfo,jbcinfo,kbcinfo,1005)
         call cntsurf(ns1006,maxbl0,maxgr0,maxseg0,ngrid,nblg,
     .                nbci0,nbcj0,nbck0,nbcidim,nbcjdim,
     .                nbckdim,ibcinfo,jbcinfo,kbcinfo,1006)
         call cntsurf(ns2016,maxbl0,maxgr0,maxseg0,ngrid,nblg,
     .                nbci0,nbcj0,nbck0,nbcidim,nbcjdim,
     .                nbckdim,ibcinfo,jbcinfo,kbcinfo,2016)
         nsurf = ns2004 + ns2014 + ns2024 + ns2034 + ns1005 + ns1006 + 
     .           ns2016
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
      ivert = 0
      if (iunst.gt.1 .or. idef_ss.gt.0) then
         do nbl = 1,nblock
          if (levelg(nbl).ge.lglobal) then
            if(abs(isktyp).eq.1) then
              call lead(nbl,lw,lw2,maxbl0)
              iskp = iskip(nbl,1)
              jskp = jskip(nbl,1)
              kskp = kskip(nbl,1)
              do i=1,idim,iskp
                do j=1,jdim,jskp
                  do k=1,kdim,kskp
                    ivert = ivert + 1
                  end do
                end do
              end do
            else if(abs(isktyp).eq.2) then
              itot = 0
              jtot = 0
              ktot = 0
              do ii = 1,500
                if(iskip(nbl,ii).eq.0) goto 300
                itot = itot + 1
              enddo
300           continue
              do ii = 1,500
                if(jskip(nbl,ii).eq.0) goto 320
                jtot = jtot + 1
              enddo
320           continue
              do ii = 1,500
                if(kskip(nbl,ii).eq.0) goto 340
                ktot = ktot + 1
              enddo
340           continue
              ivert = ivert + itot*jtot*ktot
            end if
          end if
         end do
      end if
c
      lnslave = max(ivert,lnslave)
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
         deallocate(ltot)
         deallocate(jjmax1)
         deallocate(kkmax1)
         deallocate(iiint1)
         deallocate(iiint2)    
         deallocate(n14)
         deallocate(iwrk)
         deallocate(icsinfo)
         deallocate(iv)
         deallocate(memblock)
         deallocate(lwdat)
         deallocate(nblfine)
         deallocate(mem_req_node)
         deallocate(no_of_points)
         deallocate(lw)
         deallocate(lw2)
         deallocate(nblk)
         deallocate(limblk)
         deallocate(isva)
         deallocate(nblon)
         deallocate(lig)
         deallocate(lbg)
         deallocate(iovrlp)
         deallocate(ibpntsg)
         deallocate(iipntsg)
         deallocate(rkap0g)
         deallocate(levelg)
         deallocate(igridg)
         deallocate(iflimg)
         deallocate(ifdsg)
         deallocate(iviscg)
         deallocate(jdimg)
         deallocate(kdimg)
         deallocate(idimg)
         deallocate(idiagg)
         deallocate(nblcg)
         deallocate(idegg)
         deallocate(jsg)
         deallocate(ksg)
         deallocate(isg)
         deallocate(jeg)
         deallocate(keg)
         deallocate(ieg)
         deallocate(mit)
         deallocate(jlamlog)
         deallocate(klamlog)
         deallocate(ilamlog)
         deallocate(jlamhig)
         deallocate(klamhig)
         deallocate(ilamhig)
         deallocate(iwfg)
         deallocate(utrans)
         deallocate(vtrans)
         deallocate(wtrans)
         deallocate(omegax)
         deallocate(omegay)
         deallocate(omegaz)
         deallocate(xorig)
         deallocate(yorig)
         deallocate(zorig)
         deallocate(dxmx)
         deallocate(dymx)
         deallocate(dzmx)
         deallocate(dthxmx)
         deallocate(dthymx)
         deallocate(dthzmx)
         deallocate(thetax)
         deallocate(thetay)
         deallocate(thetaz)
         deallocate(rfreqt)
         deallocate(rfreqr)
         deallocate(xorig0)
         deallocate(yorig0)
         deallocate(zorig0)
         deallocate(time2)
         deallocate(thetaxl)
         deallocate(thetayl)
         deallocate(thetazl)
         deallocate(itrans)
         deallocate(irotat)
         deallocate(idefrm)
         deallocate(bcvali)
         deallocate(bcvalj)
         deallocate(bcvalk)
         deallocate(nbci0)
         deallocate(nbcj0)
         deallocate(nbck0)
         deallocate(nbcidim)
         deallocate(nbcjdim)
         deallocate(nbckdim)
         deallocate(ibcinfo)
         deallocate(jbcinfo)
         deallocate(kbcinfo)
         deallocate(bcfilei)
         deallocate(bcfilej)
         deallocate(bcfilek)
         deallocate(ncgg)
         deallocate(nblg)
         deallocate(iemg)
         deallocate(inewgg)
         deallocate(inpl3d)
         deallocate(inpr)
         deallocate(iadvance)
         deallocate(iforce)
         deallocate(iindex)
         deallocate(iindx)
         deallocate(llimit)
         deallocate(iitmax)
         deallocate(mmcxie)
         deallocate(mmceta)
         deallocate(ncheck)
         deallocate(iifit)
         deallocate(mblkpt)
         deallocate(iic0)
         deallocate(iiorph)
         deallocate(iitoss)
         deallocate(ifiner)
         deallocate(dx)
         deallocate(dy)
         deallocate(dz)
         deallocate(dthetx)
         deallocate(dthety)
         deallocate(dthetz)
         deallocate(isav_blk)
         deallocate(isav_prd)
         deallocate(isav_pat)
         deallocate(isav_pat_b)
         deallocate(isav_dpat)
         deallocate(isav_dpat_b)
         deallocate(isav_emb)
         deallocate(mblk2nd)
         deallocate(mglevg)
         deallocate(nemgl)
         deallocate(ipl3dtmp)
         deallocate(ifrom)
         deallocate(xif1)
         deallocate(etf1)
         deallocate(xif2)
         deallocate(etf2)
         deallocate(utrnsae)
         deallocate(vtrnsae)
         deallocate(wtrnsae)
         deallocate(omgxae)
         deallocate(omgyae)
         deallocate(omgzae)
         deallocate(xorgae)
         deallocate(yorgae)
         deallocate(zorgae)
         deallocate(xorgae0)
         deallocate(yorgae0)
         deallocate(zorgae0)
         deallocate(icouple)
         deallocate(thtxae)
         deallocate(thtyae)
         deallocate(thtzae)
         deallocate(rfrqtae)
         deallocate(rfrqrae)
         deallocate(icsi)
         deallocate(icsf)
         deallocate(jcsi)
         deallocate(jcsf)
         deallocate(kcsi)
         deallocate(kcsf)
         deallocate(idfrmseg)
         deallocate(iaesurf)
         deallocate(nsegdfrm)
         deallocate(freq)
         deallocate(gmass)
         deallocate(x0)
         deallocate(gf0)
         deallocate(damp)
         deallocate(perturb)
         deallocate(aesrfdat)
         deallocate(iskip)
         deallocate(jskip)
         deallocate(kskip)
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
        if (abs(jbcinfo(nbl,mseg,1,1)) .eq. 2004 .or.
     .      abs(jbcinfo(nbl,mseg,1,1)) .eq. 2014 .or.
     .      abs(jbcinfo(nbl,mseg,1,1)) .eq. 2024 .or.
     .      abs(jbcinfo(nbl,mseg,1,1)) .eq. 2034 .or.
     .      abs(jbcinfo(nbl,mseg,1,1)) .eq. 2016) ibcjmin=1 
      end do
      do mseg=1,nbcjdim(nbl)
        if (abs(jbcinfo(nbl,mseg,1,2)) .eq. 2004 .or.
     .      abs(jbcinfo(nbl,mseg,1,2)) .eq. 2014 .or.
     .      abs(jbcinfo(nbl,mseg,1,2)) .eq. 2024 .or.
     .      abs(jbcinfo(nbl,mseg,1,2)) .eq. 2034 .or.
     .      abs(jbcinfo(nbl,mseg,1,2)) .eq. 2016) ibcjmax=1 
      end do
c
      ibckmin=0
      ibckmax=0
      do mseg=1,nbck0(nbl)
        if (abs(kbcinfo(nbl,mseg,1,1)) .eq. 2004 .or.
     .      abs(kbcinfo(nbl,mseg,1,1)) .eq. 2014 .or.
     .      abs(kbcinfo(nbl,mseg,1,1)) .eq. 2024 .or.
     .      abs(kbcinfo(nbl,mseg,1,1)) .eq. 2034 .or.
     .      abs(kbcinfo(nbl,mseg,1,1)) .eq. 2016) ibckmin=1
      end do
      do mseg=1,nbckdim(nbl)
        if (abs(kbcinfo(nbl,mseg,1,2)) .eq. 2004 .or.
     .      abs(kbcinfo(nbl,mseg,1,2)) .eq. 2014 .or.
     .      abs(kbcinfo(nbl,mseg,1,2)) .eq. 2024 .or.
     .      abs(kbcinfo(nbl,mseg,1,2)) .eq. 2034 .or.
     .      abs(kbcinfo(nbl,mseg,1,2)) .eq. 2016) ibckmax=1
      end do
c
      ibcimin=0
      ibcimax=0
      do mseg=1,nbci0(nbl)
        if (abs(ibcinfo(nbl,mseg,1,1)) .eq. 2004 .or.
     .      abs(ibcinfo(nbl,mseg,1,1)) .eq. 2014 .or.
     .      abs(ibcinfo(nbl,mseg,1,1)) .eq. 2024 .or.
     .      abs(ibcinfo(nbl,mseg,1,1)) .eq. 2034 .or.
     .      abs(ibcinfo(nbl,mseg,1,1)) .eq. 2016) ibcimin=1
      end do
      do mseg=1,nbcidim(nbl)
        if (abs(ibcinfo(nbl,mseg,1,2)) .eq. 2004 .or.
     .      abs(ibcinfo(nbl,mseg,1,2)) .eq. 2014 .or.
     .      abs(ibcinfo(nbl,mseg,1,2)) .eq. 2024 .or.
     .      abs(ibcinfo(nbl,mseg,1,2)) .eq. 2034 .or.
     .      abs(ibcinfo(nbl,mseg,1,2)) .eq. 2016) ibcimax=1
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
