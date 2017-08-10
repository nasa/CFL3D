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
      subroutine mgblk(iseq,lw,lw2,w,mgwk,timesave,wk,nwork,iw,ireq_qb,
     .                 iwk,iwork,maxbl,maxgr,maxseg,iitot,intmax,nsub1,
     .                 maxxe,ncycmax,iovrlp,lig,lbg,ibpntsg,iipntsg,
     .                 iibg,kkbg,jjbg,ibcg,dxintg,dyintg,dzintg,iiig,
     .                 jjig,kkig,qb,lwdat,nblk,nbli,limblk,isva,nblon,
     .                 nblock,levelg,igridg,iviscg,idimg,jdimg,kdimg,
     .                 jsg,ksg,isg,jeg,keg,ieg,nblcg,mit,bcvali,bcvalj,
     .                 bcvalk,nbci0,nbcidim,nbcj0,nbcjdim,nbck0,
     .                 nbckdim,ibcinfo,jbcinfo,kbcinfo,bcfilei,bcfilej,
     .                 bcfilek,utrans,vtrans,wtrans,omegax,omegay,
     .                 omegaz,xorig,yorig,zorig,dxmx,dymx,dzmx,dthxmx,
     .                 dthymx,dthzmx,thetax,thetay,thetaz,rfreqt,
     .                 rfreqr,xorig0,yorig0,zorig0,time2,thetaxl,
     .                 thetayl,thetazl,itrans,irotat,idefrm,
     .                 rms,clw,cdw,cdpw,cdvw,cxw,cyw,czw,cmxw,
     .                 cmyw,cmzw,n_clcd,clcd,nblocks_clcd,blocks_clcd,
     $                 chdw,swetw,fmdotw,cfttotw,cftmomw,
     .                 cftpw,cftvw,rmstr,nneg,ntr,
     .                 ihstry,iadvance,iforce,lfgm,resmx,imx,jmx,kmx,
     .                 vormax,ivmax,jvmax,kvmax,sx,sy,sz,stot,pav,
     .                 ptav,tav,ttav,xmav,fmdot,cfxp,cfyp,cfzp,cfdp,
     .                 cflp,cftp,cfxv,cfyv,cfzv,cfdv,cflv,cftv,cfxmom,
     .                 cfymom,cfzmom,cfdmom,cflmom,cftmom,cfxtot,
     .                 cfytot,cfztot,cfdtot,cfltot,cfttot,icsinfo,ncs,
     .                 windex,ninter,iindex,nblkpt,windx,nintr,iindx,
     .                 llimit,iitmax,mmcxie,mmceta,ncheck,iifit,mblkpt,
     .                 iic0,iiorph,iitoss,ifiner,msub1,dthetxx,dthetyy,
     .                 dthetzz,swett,clt,cdt,cxt,cyt,czt,cmxt,cmyt,
     .                 cmzt,cdpt,cdvt,mblk2nd,geom_miss,epsc0,
     .                 period_miss,epsrot,isav_blk,isav_prd,lbcprd,
     .                 isav_pat,isav_pat_b,isav_emb,lbcemb,mxbli,
     .                 maxcs,intmx,mxxe,mptch,ncgg,nblg,iemg,ngrid,
     .                 dx,dy,dz,dthetx,dthety,dthetz,isav_dpat,
     .                 isav_dpat_b,lout,ifrom,xif1,xif2,etf1,
     .                 etf2,jjmax1,kkmax1,iiint1,iiint2,nblk1,
     .                 nblk2,jimage,kimage,jte,kte,jmm,kmm,
     .                 xte,yte,zte,xmi,ymi,zmi,xmie,ymie,
     .                 zmie,sxie,seta,sxie2,seta2,xie2s,
     .                 eta2s,temp,x2,y2,z2,x1,y1,z1,ip3dsurf,
     .                 factjlo,factjhi,factklo,factkhi,nplots,
     .                 nou,bou,nbuf,ibufdim,istat2_bl,
     .                 istat2_pa,istat2_pe,istat2_em,istat_size,
     .                 nplot3d,nprint,inpl3d,inpr,bcfiles,mxbcfil,
     .                 utrnsae,vtrnsae,wtrnsae,omgxae,omgyae,omgzae,
     .                 xorgae,yorgae,zorgae,thtxae,thtyae,thtzae,
     .                 rfrqtae,rfrqrae,icsi,icsf,jcsi,jcsf,kcsi,kcsf,
     .                 freq,gmass,damp,x0,gf0,nmds,maxaes,aesrfdat,
     .                 perturb,islavept,nslave,iskip,jskip,kskip,bmat,
     .                 stm,stmi,xs,xxn,gforcn,gforcnm,gforcs,nsegdfrm,
     .                 idfrmseg,iaesurf,maxsegdg,nmaster,aehist,
     .                 timekeep,xorgae0,yorgae0,zorgae0,icouple,
     .                 iprnsurf,nblelst,iskmax,jskmax,kskmax,ue,nummem)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Advance the solution in time using multigrid acceleration
c     At each level in the multigrid procedure, all the blocks are 
c     advanced before moving to the next level.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
#if defined DIST_MPI
#     include "mpif.h"
#   ifdef DBLE_PRECSN
#      ifdef CMPLX
#        define MY_MPI_REAL MPI_DOUBLE_COMPLEX
#      else
#        define MY_MPI_REAL MPI_DOUBLE_PRECISION
#      endif
#   else
#      ifdef CMPLX
#        define MY_MPI_REAL MPI_COMPLEX
#      else
#        define MY_MPI_REAL MPI_REAL
#      endif
#   endif
      dimension istat(MPI_STATUS_SIZE)
c
#endif
c     maxbl  - maximum number of blocks
c     maxgr  - maximum number of grids
c     ifamax - maximum for total number of starting and ending indices
c              for flux accumulation in the I-, J-, and K-directions
c
      logical stop_me
c
      character*80 grid,plt3dg,plt3dq,output,residual,turbres,blomx,
     .             output2,printout,pplunge,ovrlap,patch,restrt,
     .             subres,subtur,grdmov,alphahist,errfile,preout,
     .             aeinp,aeout,sdhist,avgg,avgq
      character*241 avgq2,avgq2pert,clcds,clcdp,output_dir
      character*120 bou(ibufdim,nbuf)
      character*80  bcfiles(mxbcfil)
      character*32 basedesired
c
      integer lout(msub1),xif1(msub1),xif2(msub1),etf1(msub1),
     .        etf2(msub1),ifrom(msub1)
      integer bcfilei,bcfilej,bcfilek
c
      dimension nou(nbuf)
      dimension istat2_bl(istat_size,mxbli*5),
     .          istat2_pa(istat_size,intmax*nsub1*3),
     .          istat2_em(istat_size,lbcemb*3),
     .          istat2_pe(istat_size,lbcprd*5)
      dimension w(mgwk),wk(nwork),lw(65,maxbl),lw2(43,maxbl),iwk(iwork)
      dimension iw(maxbl*7*3),timesave(maxbl),lwdat(maxbl,maxseg,6)
      dimension inpl3d(nplots,11),inpr(nplots,11)
      dimension deltat(5),nptsr(20)
      dimension mblk2nd(maxbl),ireq_qb(maxbl*13),mytag_qb(13)
      dimension iovrlp(maxbl),lig(maxbl),lbg(maxbl),ibpntsg(maxbl,4),
     .          iipntsg(maxbl),iibg(iitot),kkbg(iitot),jjbg(iitot),
     .          ibcg(iitot),dxintg(iitot),dyintg(iitot),dzintg(iitot),
     .          iiig(iitot),jjig(iitot),kkig(iitot),qb(iitot,5,3)
      dimension nblk(2,mxbli),limblk(2,6,mxbli),isva(2,2,mxbli),
     .          nblon(mxbli)
      dimension levelg(maxbl),igridg(maxbl),iviscg(maxbl,3),
     .          jdimg(maxbl),kdimg(maxbl),idimg(maxbl),nblcg(maxbl),
     .          jsg(maxbl),ksg(maxbl),isg(maxbl),jeg(maxbl),keg(maxbl),
     .          ieg(maxbl),mit(5,maxbl)
      dimension bcvali(maxbl,maxseg,12,2),
     .          bcvalj(maxbl,maxseg,12,2),bcvalk(maxbl,maxseg,12,2),
     .          nbci0(maxbl),nbcidim(maxbl),nbcj0(maxbl),nbcjdim(maxbl),
     .          nbck0(maxbl),nbckdim(maxbl),ibcinfo(maxbl,maxseg,7,2),
     .          jbcinfo(maxbl,maxseg,7,2),kbcinfo(maxbl,maxseg,7,2)
      dimension bcfilei(maxbl,maxseg,2),bcfilej(maxbl,maxseg,2),
     .          bcfilek(maxbl,maxseg,2)
      dimension utrans(maxbl),vtrans(maxbl),wtrans(maxbl),omegax(maxbl),
     .          omegay(maxbl),omegaz(maxbl),xorig(maxbl),yorig(maxbl),
     .          zorig(maxbl),dxmx(maxbl),dymx(maxbl),dzmx(maxbl),
     .          dthymx(maxbl),dthzmx(maxbl),dthxmx(maxbl),
     .          thetax(maxbl),thetay(maxbl),thetaz(maxbl),
     .          rfreqt(maxbl),rfreqr(maxbl),xorig0(maxbl),
     .          yorig0(maxbl),zorig0(maxbl),time2(maxbl),thetaxl(maxbl),
     .          thetayl(maxbl),thetazl(maxbl),itrans(maxbl),
     .          irotat(maxbl),idefrm(maxbl)
      dimension utrnsae(maxbl,maxsegdg),vtrnsae(maxbl,maxsegdg),
     .          wtrnsae(maxbl,maxsegdg),omgxae(maxbl,maxsegdg),
     .          omgyae(maxbl,maxsegdg),omgzae(maxbl,maxsegdg),
     .          xorgae(maxbl,maxsegdg),yorgae(maxbl,maxsegdg),
     .          zorgae(maxbl,maxsegdg),thtxae(maxbl,maxsegdg),
     .          thtyae(maxbl,maxsegdg),thtzae(maxbl,maxsegdg),
     .          rfrqtae(maxbl,maxsegdg),rfrqrae(maxbl,maxsegdg)
      dimension xorgae0(maxbl,maxsegdg),yorgae0(maxbl,maxsegdg),
     .          zorgae0(maxbl,maxsegdg),icouple(maxbl,maxsegdg)
      dimension icsi(maxbl,maxsegdg),icsf(maxbl,maxsegdg),
     .          jcsi(maxbl,maxsegdg),jcsf(maxbl,maxsegdg),
     .          kcsi(maxbl,maxsegdg),kcsf(maxbl,maxsegdg)
      dimension nsegdfrm(maxbl),idfrmseg(maxbl,maxsegdg),
     .          iaesurf(maxbl,maxsegdg)
      dimension freq(nmds,maxaes),gmass(nmds,maxaes),x0(2*nmds,maxaes),
     .          gf0(2*nmds,maxaes),damp(nmds,maxaes),
     .          perturb(nmds,maxaes,4),aehist(ncycmax,3,nmds,maxaes)
      dimension aesrfdat(5,maxaes),islavept(nslave,nmaster,5)
      dimension ue(3*nslave)
      dimension iskip(maxbl,500),jskip(maxbl,500),kskip(maxbl,500)
      dimension bmat(2*nmds,2*nmds,maxaes),gforcn(2*nmds,maxaes),
     .          gforcnm(2*nmds,maxaes),gforcs(2*nmds,maxaes),
     .          stm(2*nmds,2*nmds,maxaes),stmi(2*nmds,2*nmds,maxaes),
     .          xs(2*nmds,maxaes),xxn(2*nmds,maxaes)
      dimension ncgg(maxgr),nblg(maxgr),iemg(maxgr)
      integer blocks_clcd
      dimension rms(ncycmax),clw(ncycmax),cdw(ncycmax),cdpw(ncycmax),
     .          cdvw(ncycmax),cxw(ncycmax),cyw(ncycmax),czw(ncycmax),
     .          cmxw(ncycmax),cmyw(ncycmax),cmzw(ncycmax),
     .          clcd(2,n_clcd,ncycmax), blocks_clcd(2,nblocks_clcd),
     .          chdw(ncycmax),swetw(ncycmax),fmdotw(ncycmax),
     .          cfttotw(ncycmax),cftmomw(ncycmax),cftpw(ncycmax),
     .          cftvw(ncycmax),rmstr(ncycmax,nummem),
     .          nneg(ncycmax,nummem),timekeep(ncycmax)
      dimension iadvance(maxbl),iforce(maxbl)
      dimension resmx(maxbl),imx(maxbl),jmx(maxbl),kmx(maxbl)
      dimension vormax(maxbl),ivmax(maxbl),jvmax(maxbl),kvmax(maxbl)
      dimension sx(maxcs),sy(maxcs),sz(maxcs),stot(maxcs),
     .          pav(maxcs),ptav(maxcs),tav(maxcs),ttav(maxcs),
     .          xmav(maxcs),fmdot(maxcs),
     .          cfxp(maxcs),cfyp(maxcs),cfzp(maxcs),
     .          cfdp(maxcs),cflp(maxcs),cftp(maxcs),
     .          cfxv(maxcs),cfyv(maxcs),cfzv(maxcs),
     .          cfdv(maxcs),cflv(maxcs),cftv(maxcs),
     .          cfxmom(maxcs),cfymom(maxcs),cfzmom(maxcs),
     .          cfdmom(maxcs),cflmom(maxcs),cftmom(maxcs),
     .          cfxtot(maxcs),cfytot(maxcs),cfztot(maxcs),
     .          cfdtot(maxcs),cfltot(maxcs),cfttot(maxcs),
     .          icsinfo(maxcs,9)
      dimension windex(maxxe,2),iindex(intmax,6*nsub1+9),
     .          nblkpt(maxxe)
      dimension windx(mxxe,2),iindx(intmx,6*msub1+9),
     .          llimit(intmx),iitmax(intmx),mmcxie(intmx),
     .          mmceta(intmx),ncheck(maxbl),iifit(intmx),
     .          mblkpt(mxxe),iic0(intmx),iiorph(intmx),iitoss(intmx),
     .          ifiner(intmx)
      dimension dthetxx(intmax,nsub1),dthetyy(intmax,nsub1),
     .          dthetzz(intmax,nsub1)
      dimension swett(maxbl),clt(maxbl),cdt(maxbl),cxt(maxbl),
     .          cyt(maxbl),czt(maxbl),cmxt(maxbl),cmyt(maxbl),
     .          cmzt(maxbl),cdpt(maxbl),cdvt(maxbl)
      dimension geom_miss(2*mxbli),period_miss(lbcprd)
      dimension isav_blk(2*mxbli,17),isav_prd(lbcprd,12)
      dimension isav_pat(intmax,17),isav_pat_b(intmax,nsub1,6)
      dimension isav_dpat(intmx,17),isav_dpat_b(intmx,msub1,6)
      dimension isav_emb(lbcemb,12)
      dimension dx(intmx,msub1),dy(intmx,msub1),dz(intmx,msub1),
     .          dthetx(intmx,msub1),dthety(intmx,msub1),
     .          dthetz(intmx,msub1)
      dimension jjmax1(nsub1),kkmax1(nsub1),iiint1(nsub1),iiint2(nsub1)
      dimension jimage(msub1,mptch+2,mptch+2),
     .          kimage(msub1,mptch+2,mptch+2)
      dimension jte(msub1),kte(msub1)
      dimension jmm(mptch+2),kmm(mptch+2)
      dimension nblk1(mptch+2),nblk2(mptch+2)
      dimension xte(mptch+2,mptch+2,msub1),yte(mptch+2,mptch+2,msub1),
     .          zte(mptch+2,mptch+2,msub1)
      dimension xmi(mptch+2,mptch+2,msub1),ymi(mptch+2,mptch+2,msub1),
     .          zmi(mptch+2,mptch+2,msub1)
      dimension xmie(mptch+2,mptch+2,msub1),ymie(mptch+2,mptch+2,msub1),
     .          zmie(mptch+2,mptch+2,msub1)
      dimension sxie(mptch+2,mptch+2,msub1),seta(mptch+2,mptch+2,msub1),
     .          sxie2(mptch+2,mptch+2),seta2(mptch+2,mptch+2)
      dimension xie2s(mptch+2,mptch+2),eta2s(mptch+2,mptch+2)
      dimension temp((mptch+2)*(mptch+2))
      dimension x2(mptch+2,mptch+2),y2(mptch+2,mptch+2),
     .          z2(mptch+2,mptch+2)
      dimension x1(mptch+2,mptch+2),y1(mptch+2,mptch+2),
     .          z1(mptch+2,mptch+2)
      dimension factjlo(intmx,msub1),factjhi(intmx,msub1),
     .          factklo(intmx,msub1),factkhi(intmx,msub1)
      dimension istop_data(10)
      dimension nblelst(maxbl,2)
      dimension iskmax(maxbl)
      dimension jskmax(maxbl)
      dimension kskmax(maxbl)
      dimension idimt(maxbl),jdimt(maxbl),kdimt(maxbl)
      dimension nnegb(nummem),rmstb(nummem),sumn(nummem),negn(nummem)
c
      common /cfl/ dt0,dtold
      common /elastic/ ndefrm,naesrf
      common /cgns/ icgns,iccg,ibase,nzones,nsoluse,irind,jrind,krind
      common /unit5/ iunit5
      common /chk/ ichk
      common /account/ iaccnt,ioutsub
      common /ginfo/ jdim,kdim,idim,jj2,kk2,ii2,nblc,js,ks,is,je,ke,ie,
     .        lq,lqj0,lqk0,lqi0,lsj,lsk,lsi,lvol,ldtj,lx,ly,lz,lvis,
     .        lsnk0,lsni0,lq1,lqr,lblk,lxib,lsig,lsqtq,lg,
     .        ltj0,ltk0,lti0,lxkb,lnbl,lvj0,lvk0,lvi0,lbcj,lbck,lbci,
     .        lqc0,ldqc0,lxtbi,lxtbj,lxtbk,latbi,latbj,latbk,
     .        lbcdj,lbcdk,lbcdi,lxib2,lux,lcmuv,lvolj0,lvolk0,lvoli0,
     .        lxmdj,lxmdk,lxmdi,lvelg,ldeltj,ldeltk,ldelti,
     .        lxnm2,lynm2,lznm2,lxnm1,lynm1,lznm1,lqavg
      common /ginfo2/ lq2avg,iskip_blocks,inc_2d(3),inc_coarse(3)
      common /avgdata/ xnumavg,iteravg,xnumavg2,ipertavg,iclcd,isubit_r
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /mgv/ epsssc(3),epsssr(3),issc,issr
      common /nfablk/ nfajki(3)
      common /sklton/ isklton
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /wrbl/ nwrest
      common /moov/movie,nframes,icall1,lhdr,icoarsemovie,i2dmovie
      common /filenam/ grid,plt3dg,plt3dq,output,residual,turbres,blomx,
     .                 output2,printout,pplunge,ovrlap,patch,restrt,
     .                 subres,subtur,grdmov,alphahist,errfile,preout,
     .                 aeinp,aeout,sdhist,avgg,avgq
      common /filenam2/ avgq2,avgq2pert,clcds,clcdp,output_dir
      common /alphait/ ialphit,cltarg,rlxalph,dalim,dalpha,icycupdt
      common /precond/ cprec,uref,avn
      common /mydist2/ nnodes,myhost,myid,mycomm
      common /maxiv/ ivmx
      common /motionmc/ xmc0,ymc0,zmc0,utransmc,vtransmc,wtransmc,
     .                  omegaxmc,omegaymc,omegazmc,xorigmc,yorigmc,
     .                  zorigmc,xorig0mc,yorig0mc,zorig0mc,thetaxmc,
     .                  thetaymc,thetazmc,dxmxmc,dymxmc,dzmxmc,
     .                  dthxmxmc,dthymxmc,dthzmxmc,rfreqtmc,
     .                  rfreqrmc,itransmc,irotatmc,time2mc
      common /trim/ dmtrmn,dmtrmnm,dlcln,dlclnm,trtol,cmy,cnw,alf0,
     .              alf1,dzdt,thtd0,thtd1,zrg0,zrg1,dtrmsmx,dtrmsmn,
     .              dalfmx,ddtmx,ddtrm0,ddtrm1,itrmt,itrminc,fp(4,4),
     .              tp(4,4),zlfct,epstr,relax,ittrst 
      common /rigidbody/ irigb,irbtrim
      common /turbconv/ cflturb(7),edvislim,iturbprod,nsubturb,nfreeze,
     .                  iwarneddy,itime2read,itaturb,tur1cut,tur2cut,
     .                  iturbord,tur1cutlev,tur2cutlev
      common /deformz/ beta1,beta2,alpha1,alpha2,isktyp,negvol,meshdef,
     .                 nsprgit,ndgrd,ndwrt 
c
c      time step cycle:
c
c      nt       - time step counter
c
c      ntstep   - number of time steps to take
c                 = 1 for steady state
c
c      isklton  - flag to generate verbose output to main output
c                 file during first timestep/multigrid cycle
c                 = 1  verbose output
c                 = 0  minimal output
c
c
c
c      multigrid cycle:
c
c      mseq     - number of sequenced global solutions
c
c      iseq     - counter for mesh sequencing
c                 = 1 : coarsest global solution
c                   .
c                   .
c                   .
c                 = mseq : finest global solution
c
c      note: mesh embedding allowed only for finest global solution (iseq=mseq)
c
c      levelt   - starting level for cycle of sequence iseq
c
c      levelb   - ending level    "          "
c
c      mgflag   - multigrid flag
c                 = 0  no multigrid
c                 = 1  multigrid on global grids only
c                 = 2  multigrid on global and embedded grids
c
c      levt     - starting level for cycle for sequence iseq
c
c      levb     - ending level    "          "
c
c      lfgm     - level of finest global mesh
c
c      level    - counter for grid levels
c                 = 1 : coarsest grids
c                   .
c                   .
c                   .
c                 = mseq : finest global grids
c                   .
c                   .
c                   .
c                 = mseq + nembed : finest embedded grids
c
c      nembed   - number of embedded levels for sequence iseq
c
c      lglobal  - level of finest global grids for sequence iseq
c
c      kxpand   - flag to indicate direction of multigrid process
c                 =  1 : downward (residual/restriction) leg
c                 = -1 : upward (prolongation) leg
c
c      nit      - number of solution updates to do at each level
c                 of the multigrid cycle
c
c      mtt      - number of additional solution updates to do 
c                 on upward leg of the multigrid cycle
c
c      ioutsub  - flag to indicate at which point in the subiteration
c                 loop to monitor global residual and force/moment data
c                 (this is the residual and force/moment data that gets
c                 written out to the main output file - unit 11 - and the 
c                 main convergence history file - unit 12. ioutsub is only
c                 applicable to time accurate cases with subiteration)
c                 =    1 : monitor at start of subiteration loop
c                 = ncyc : monitor at end of subiteration loop
c
c      kode=1   mode=0   standard update
c                        res=r(q)
c                        l(dq)=res
c
c      kode=2   mode=1   first time on lower level multigrid
c                        q1=q
c                        res=r(q)
c                        qr=qr-res
c                        res=res+qr
c                        l(dq)=res
c
c      kode=1   mode=1   subsequent time on lower level multigrid
c                        res=r(q)
c                        res=res+qr
c                        l(dq)=res
c
      levt  = levelt(iseq)
      levb  = levelb(iseq)
      nitfo = nitfo1(iseq)
      ncyc  = ncyc1 (iseq)
c
      stop_me = .false.
c
c      check storage for nsm and mit arrays
c
      if (levt-levb+1.gt.5) then
         write(11,1000) 
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
c
      lglobal = lfgm-(mseq-iseq)
      nembed  = levt-lglobal
c
c     find first (starting) block at global level, nblstag
c
      do 5000 nbl=1,nblock
      if (mblk2nd(nbl).eq.myid .or. myid.eq.myhost) then
         if (lglobal.ne.levelg(nbl)) go to 5000
         nblstag = nbl
         go to 5001
      end if
 5000 continue
 5001 continue
c
c     find last (ending) block at global level, nblendg
c
      do 5010 nbl=1,nblock
      if (mblk2nd(nbl).eq.myid .or. myid.eq.myhost) then
         if (lglobal.ne.levelg(nbl)) go to 5010
         nblendg = nbl
      end if
 5010 continue
c
c     find first (starting) block at top level, nblstat
c
      do 5020 nbl=1,nblock
      if (mblk2nd(nbl).eq.myid .or. myid.eq.myhost) then
         if (levt.ne.levelg(nbl)) go to 5020
         nblstat = nbl
         go to 5021
      end if
 5020 continue
 5021 continue
c
#if defined DIST_MPI
c     all processors need to know the host's value of nblstat to 
c     get same order of ramped cfl/dt output as sequential code 
c
      do 5025 nbl=1,nblock
      if (levt.ne.levelg(nbl)) go to 5025
      nblstat_h = nbl
      go to 5026
 5025 continue
 5026 continue
#else
      nblstat_h = nblstat
#endif
c
c     find last (ending) block at top level, nblendt
c
      do 5030 nbl=1,nblock
      if (mblk2nd(nbl).eq.myid .or. myid.eq.myhost) then
         if (levt.ne.levelg(nbl)) go to 5030
         nblendt = nbl
      end if
 5030 continue
c
#if defined DIST_MPI
c     set baseline tag values
c
      ioffset = maxbl
      itag_qc = 1
      itag_q1 = itag_qc + ioffset
      itag_qv = itag_q1 + ioffset 
      itag_qt = itag_qv + ioffset
      itag_qr = itag_qt + ioffset
c
      mytag_qb(1) = 0
      do ll=2,13
        mytag_qb(ll) = mytag_qb(ll-1) + maxbl
      end do
#endif
c
c     initialize mismatch check arrays
c
      do nn=1,2*mxbli
         geom_miss(nn) = 0.
      end do
      do nn=1,lbcprd
         period_miss(nn) = 0.
      end do
c
c********************************************************************
c
c     begin time advancement
c
c********************************************************************
c
c     output plot3d file of initial flow field for animation
c
c
      ncall = 1
      if (real(dt).gt.0. and. movie.lt.0) then
         isklton = 0
         lhdr    = 0
         if (iunst .gt. 0) then
c
c           calculate grid speeds and metrics for use in boundary
c           conditions but do not change grid position (iupdat=0)
c           since time (but not the grid) is updated before the
c           restart file is written, must decrement block time counter
c           by dt in order to get the right temporal metrics
c
            do 108 nbl=1,nblock
            timesave(nbl) = time2(nbl)
            time2(nbl)    = time2(nbl) - dt
 108        continue
            iupdat = 0
            level  = levelt(iseq)
            irigb0   = 0
            irbtrim0 = 0
            call updateg(lw,lw2,w,mgwk,wk,nwork,iupdat,iseq,maxbl,
     .                   maxgr,maxseg,nbci0,nbcj0,nbck0,nbcidim,
     .                   nbcjdim,nbckdim,ibcinfo,jbcinfo,kbcinfo,
     .                   nblock,levelg,igridg,utrans,vtrans,wtrans,
     .                   omegax,omegay,omegaz,xorig,yorig,zorig,
     .                   thetax,thetay,thetaz,rfreqt,rfreqr,xorig0,
     .                   yorig0,zorig0,time2,thetaxl,thetayl,thetazl,
     .                   itrans,irotat,idefrm,ncgg,iadvance,nou,
     .                   bou,nbuf,ibufdim,myid,myhost,mycomm,mblk2nd,
     .                   irigb0,irbtrim0,nt)
            do 109 nbl=1,nblock
            time2(nbl) = timesave(nbl)
 109        continue
c
c           get interpolation data for dynamic boundaries (for qout)
c
            iwk1   = 1
            iwk2   = iwk1  + 65*maxbl
            iwk3   = iwk2  + 3*intmx
            iwork1 = iwork - iwk3
            if (iwork1.lt.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),'(''stopping...not enough integer '',
     .                    ''work space for subroutine dynptch'')')
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
            do iii = 1,iwk3
               iwk(iii) = 0
            end do
            call dynptch(lw,lw2,w,mgwk,wk,nwork,ncall,
     .                   maxgr,maxbl,msub1,intmx,mxxe,mptch,jdimg,
     .                   kdimg,idimg,xorig,yorig,zorig,nblock,ngrid,
     .                   levelg,ncgg,nblg,windex,ninter,iindex,nblkpt,
     .                   windx,nintr,iindx,mblkpt,llimit,iitmax,
     .                   mmcxie,mmceta,ncheck,iifit,iic0,
     .                   iiorph,iitoss,ifiner,factjlo,factjhi,
     .                   factklo,factkhi,dx,dy,dz,dthetx,dthety,
     .                   dthetz,dthetxx,dthetyy,dthetzz,
     .                   isav_dpat,isav_dpat_b,intmax,maxxe,nsub1,
     .                   iwk(iwk1),iwk(iwk2),lout,ifrom,xif1,xif2,etf1,
     .                   etf2,jjmax1,kkmax1,iiint1,iiint2,nblk1,
     .                   nblk2,jimage,kimage,jte,kte,jmm,kmm,
     .                   xte,yte,zte,xmi,ymi,zmi,xmie,ymie,
     .                   zmie,sxie,seta,sxie2,seta2,xie2s,
     .                   eta2s,temp,x2,y2,z2,x1,y1,z1,
     .                   myid,myhost,mycomm,mblk2nd,nou,bou,nbuf,
     .                   ibufdim,igridg,iemg)
#if defined DIST_MPI
c
c           rebroadcast ninter, since dynptch can change it (only needed 
c           after first call to dynptch)
c
            if (ncall .eq. 1) then
               call MPI_Bcast (ninter,1,MPI_INTEGER,myhost,mycomm,ierr)
            end if
#endif
            ncall = ncall + 1
         end if
         iwk1 = 1
         iwk2 = iwk1 + 65*maxbl
         iwk3 = iwk2 + nsub1
         iwk4 = iwk3 + nsub1
         iwk5 = iwk4 + nsub1
         iwk6 = iwk5 + nsub1
         iwork1 = iwork - iwk6
         if (iwork1.lt.0) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''stopping...not enough integer '',
     .                       ''work space for subroutine pre_bc'')')
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         do iii = 1,iwk6
            iwk(iii) = 0
         end do
         call pre_bc(lw,lw2,iwk(iwk1),maxbl,maxgr,maxseg,ninter,intmax,
     .               nsub1,iindex,isav_pat,iwk(iwk2),
     .               iwk(iwk3),iwk(iwk4),iwk(iwk5),mxbli,nbli,limblk,
     .               isva,nblon,nblk,lbcprd,isav_prd,bcvali,bcvalj,
     .               bcvalk,nblg,lbcemb,iemg,igridg,isav_emb,
     .               iviscg,jdimg,kdimg,idimg,nbci0,nbcj0,nbck0,
     .               nbcjdim,nbckdim,nbcidim,ibcinfo,jbcinfo,
     .               kbcinfo,iadvance,myid,myhost,mycomm,mblk2nd,
     .               nou,bou,nbuf,ibufdim,isav_pat_b,levelg,nblock,
     .               isav_blk)
         iwk1   = 1
         iwk2   = iwk1  + 3*nplots
         iwk3   = iwk2  + maxbl
         iwork1 = iwork - iwk3 + 1
         if (iwork1.le.0) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''stopping...not enough integer '',
     .                       ''work space for subroutine qout'')')
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         do iii = 1,iwk3
            iwk(iii) = 0
         end do
         do iii = 1,nwork
            wk(iii) = 0.0
         end do
         call qout(iseq,lw,lw2,w,mgwk,wk,nwork,
     .             nplots,iovrlp,iibg,kkbg,jjbg,ibcg,lbg,
     .             ibpntsg,qb,lwdat,nbci0,nbcj0,nbck0,
     .             nbcidim,nbcjdim,nbckdim,jbcinfo,kbcinfo,
     .             ibcinfo,bcfilei,bcfilej,bcfilek,itrans,
     .             irotat,idefrm,nblock,levelg,igridg,iviscg,jdimg,
     .             kdimg,idimg,nblg,clw,ncycmax,nplot3d,
     .             inpl3d,ip3dsurf,nprint,inpr,iadvance,mycomm,
     .             myid,myhost,mblk2nd,nou,bou,nbuf,ibufdim,maxbl,
     .             maxgr,maxseg,iitot,jsg,ksg,isg,jeg,keg,ieg,
     .             ninter,windex,iindex,nblkpt,intmax,nsub1,maxxe,
     .             nblk,nbli,limblk,isva,nblon,mxbli,thetay,
     .             iwk(iwk1),iwk(iwk2),iwk(iwk3),iwork1,xorig,yorig,
     .             zorig,period_miss,geom_miss,epsc0,epsrot,
     .             isav_blk,isav_pat,isav_pat_b,isav_emb,isav_prd,
     .             lbcprd,lbcemb,dthetxx,dthetyy,dthetzz,nblcg,lfgm,
     .             istat2_bl,istat2_pa,istat2_pe,istat2_em,
     .             istat_size,vormax,ivmax,jvmax,kvmax,bcfiles,
     .             mxbcfil,iprnsurf,nummem)
         do iii = 1,iwk3
            iwk(iii) = 0
         end do
         do iii = 1,nwork
            wk(iii) = 0.0
         end do
         nframes = nframes + 1
      end if
c
      if (myid.eq.myhost) then
         write(11,1010) iseq
         if (real(dt).lt.0.) then
            write(11,1020)
         else
            if (abs(ita).eq.1) then
               write(11,1030) ntstep
            end if
            if (abs(ita).eq.2) then
               write(11,1040) ntstep
            end if
            if (ita.lt.0) then
               write(11,1050)
            end if
         end if
         if (cprec.ne.0.) then
            write(11,1055) cprec,uref,avn
         end if 
      end if
c
c     isklton=1/0 for verbose/minumal output during first timestep/multigrid
c     cycle
c
      isklton = 1
c
c     ioutsub = 1/ncyc for global convergence history output at beginning/end 
c     of subiteration loop
c
      ioutsub = ncyc
c
c     loop over time steps
c
      do 6000 nt=1,ntstep
c
c     check for user stop file - will then stop on next time step
c
      if (myid.eq.myhost) then
        inquire (file='stop',exist=stop_me)
        if (stop_me) then
c         must set ntstep to current time step AND set number of cycles on all
c         higher levels to zero
          ntstep = nt
          do is  = iseq+1,mseq
            ncyc1(is) = 0
          end do
          write(11,*)
          write(11,'(a,i5)')
     .    'Stopping: user-invoked stop file detected at time step ',
     .    ntstep
          write(11,*)
        end if
      end if
#if defined DIST_MPI
      istop_data(1) = ntstep
      do is =1,mseq
        istop_data(is+1) = ncyc1(is)
      end do
      call MPI_Bcast (istop_data,mseq+1,MPI_INTEGER,myhost,mycomm,ierr)
      ntstep = istop_data(1)
      do is =1,mseq
        ncyc1(is) = istop_data(is+1)
      end do
#endif

c
      itatemp=ita
c     always do 1st order time 1st iteration if restarting and changing time step
      if (nt .eq. 1 .and. dt .ne. dtold) then
        ita=sign(1,ita)
      end if
c
      if (nt.gt.1) isklton = 0
c
c     start on top level
      level  = levelt(iseq)
c
c     aeroelastic predictor step
c
      if (iunst.gt.1 .and. naesrf.gt.0) then
#if defined DIST_MPI
         if (irigb .eq. 1 .or. irbtrim .eq. 1) then
            alf1   = thetay(nblstag)
            xorgrb = xorig(nblstag)
            yorgrb = yorig(nblstag)
            zorgrb = zorig(nblstag)
            call MPI_Bcast(clw(ntt),1,MY_MPI_REAL,myhost,
     .                     mycomm,ierr)
            call MPI_Bcast(cmyw(ntt),1,MY_MPI_REAL,myhost,
     .                     mycomm,ierr)
         end if
#endif
c
c    nblstat is the first block at the top level (levelt) on a
c    particular processor; all blocks at this level
c    are assumed to have the same value of thetay.
c
         call ae_pred(aesrfdat,stm,stmi,gforcn,gforcnm,xs,xxn,
     .                x0,perturb,cmyw(ntt),clw(ntt),xorgrb,yorgrb,
     .                zorgrb,nmds,maxaes,irbtrim,maxbl,myid)
c
c        rigid body trim predictor step
c
         if (irbtrim.gt.0) irigb = 0
         if (irigb.gt.0) then
            call rb_pred(nt)
         end if
         if (irigb .eq. 1 .or. irbtrim .eq. 1) then
            do nbl = 1,nblock
              thetayl(nbl)= thetay(nbl)
              thetay(nbl) = alf1
              omegay(nbl) = (thetay(nbl)-thetayl(nbl))/dt
            enddo
         end if
      end if
c
c     update grid position, speed, and dynamic patch interpolation
c     coefficients, if needed
c
      if (iunst .gt. 0) then
c
c        store off old grid position for 2nd order (in time) update
c        of deforming grid speeds
c
c
c      if meshdef = 1, deform mesh and calculate volumes and metrics only
c
       if(meshdef.eq.1) then
         if(myid.eq.myhost) then
           write(11,10229) nt
10229      format(' Deforming mesh and calculating volumes and metrics,'
     .             ,' time step ',i5)
         end if
       end if

         do nbl = 1,nblock
            if (myid.eq.mblk2nd(nbl) .and. (levelg(nbl).ge.lglobal .and.
     .         levelg(nbl).le.levelt(iseq))) then
               if (idefrm(nbl) .gt. 0) then
                  call lead(nbl,lw,lw2,maxbl)
                  mdim=jdim*kdim*idim
                  do lll=1,mdim
                     w(lxnm1+lll-1) = w(lx+lll-1)
                     w(lynm1+lll-1) = w(ly+lll-1)
                     w(lznm1+lll-1) = w(lz+lll-1)
                  end do
               end if
            end if
         end do
c
         iupdat = 1
         call updateg(lw,lw2,w,mgwk,wk,nwork,iupdat,iseq,maxbl,
     .                maxgr,maxseg,nbci0,nbcj0,nbck0,nbcidim,
     .                nbcjdim,nbckdim,ibcinfo,jbcinfo,kbcinfo,
     .                nblock,levelg,igridg,utrans,vtrans,wtrans,
     .                omegax,omegay,omegaz,xorig,yorig,zorig,
     .                thetax,thetay,thetaz,rfreqt,rfreqr,xorig0,
     .                yorig0,zorig0,time2,thetaxl,thetayl,thetazl,
     .                itrans,irotat,idefrm,ncgg,iadvance,nou,
     .                bou,nbuf,ibufdim,myid,myhost,mycomm,mblk2nd,
     .                irigb,irbtrim,nt)
c
c        reset grid position if grid has been translated/rotated
c        past limits
c
         iwk1 = maxbl
         iwork1 = iwork - iwk1
         if (iwork1.lt.0) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''stopping...not enough integer '',
     .                       ''work space for subroutine resetg'')')
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         do iii = 1,iwk1
            iwk(iii) = 0
         end do
         call resetg(lw,lw2,w,mgwk,wk,nwork,iupdat,iseq,maxbl,
     .               maxgr,maxseg,nbci0,nbcj0,nbck0,nbcidim,
     .               nbcjdim,nbckdim,ibcinfo,jbcinfo,kbcinfo,
     .               nblock,levelg,igridg,utrans,vtrans,wtrans,
     .               omegax,omegay,omegaz,xorig,yorig,zorig,
     .               thetax,thetay,thetaz,rfreqt,rfreqr,xorig0,
     .               yorig0,zorig0,time2,thetaxl,thetayl,thetazl,
     .               itrans,irotat,idefrm,ncgg,iadvance,
     .               dxmx,dymx,dzmx,dthymx,dthzmx,dthxmx,
     .               iitot,iovrlp,lig,lbg,iipntsg,ibpntsg,
     .               qb,iibg,kkbg,jjbg,ibcg,nou,bou,nbuf,ibufdim,
     .               myid,myhost,mycomm,mblk2nd,iwk,nt)
         if (iunst .gt. 1) then
            iwk1 = 2*maxbl
            iwork1 = iwork - iwk1
            if (iwork1.lt.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),'(''stopping...not enough'',
     .              ''integer work space for subroutine updatedg'')')
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
            do iii = 1,iwk1
               iwk(iii) = 0
            end do
            iupdat = 1
            call updatedg(lw,lw2,w,mgwk,wk,nwork,iupdat,iseq,maxbl,
     .                    maxgr,maxseg,nbci0,nbcj0,nbck0,nbcidim,
     .                    nbcjdim,nbckdim,ibcinfo,jbcinfo,kbcinfo,
     .                    nblock,levelg,igridg,idefrm,ncgg,iadvance,nou,
     .                    bou,nbuf,ibufdim,myid,myhost,mycomm,mblk2nd,
     .                    utrnsae,vtrnsae,wtrnsae,omgxae,omgyae,omgzae,
     .                    xorgae,yorgae,zorgae,thtxae,thtyae,thtzae,
     .                    rfrqtae,rfrqrae,icsi,icsf,jcsi,jcsf,
     .                    kcsi,kcsf,freq,gmass,damp,x0,gf0,nmds,maxaes,
     .                    aesrfdat,perturb,itrans,irotat,islavept,
     .                    nslave,iskip,jskip,kskip,xs,xxn,nsegdfrm,
     .                    idfrmseg,iaesurf,maxsegdg,iwk,nmaster,nt,
     .                    xorig,yorig,zorig,xorgae0,yorgae0,zorgae0,
     .                    icouple,iwk(maxbl+1),nnodes,nblelst,iskmax,
     .                    jskmax,kskmax,ue)
         end if
      end if
      if (iunst.gt.0) then
         iwk1     = 1
         iwk2     = iwk1  + 65*maxbl
         iwk3     = iwk2  + 3*intmx
         iwork1   = iwork - iwk3
         if (iwork1.lt.0) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''stopping...not enough integer '',
     .                 ''work space for subroutine dynptch'')')
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         do iii = 1,iwk3
            iwk(iii) = 0
         end do
         call dynptch(lw,lw2,w,mgwk,wk,nwork,ncall,
     .                maxgr,maxbl,msub1,intmx,mxxe,mptch,jdimg,
     .                kdimg,idimg,xorig,yorig,zorig,nblock,ngrid,
     .                levelg,ncgg,nblg,windex,ninter,iindex,nblkpt,
     .                windx,nintr,iindx,mblkpt,llimit,iitmax,
     .                mmcxie,mmceta,ncheck,iifit,iic0,
     .                iiorph,iitoss,ifiner,factjlo,factjhi,
     .                factklo,factkhi,dx,dy,dz,dthetx,dthety,
     .                dthetz,dthetxx,dthetyy,dthetzz,
     .                isav_dpat,isav_dpat_b,intmax,maxxe,nsub1,
     .                iwk(iwk1),iwk(iwk2),lout,ifrom,xif1,xif2,etf1,
     .                etf2,jjmax1,kkmax1,iiint1,iiint2,nblk1,
     .                nblk2,jimage,kimage,jte,kte,jmm,kmm,
     .                xte,yte,zte,xmi,ymi,zmi,xmie,ymie,
     .                zmie,sxie,seta,sxie2,seta2,xie2s,
     .                eta2s,temp,x2,y2,z2,x1,y1,z1,
     .                myid,myhost,mycomm,mblk2nd,nou,bou,nbuf,
     .                ibufdim,igridg,iemg)
#if defined DIST_MPI
c
c        rebroadcast ninter, since dynptch can change it (only needed 
c        after first call to dynptch)
c
         if (ncall .eq. 1) then
            call MPI_Bcast (ninter,1,MPI_INTEGER,myhost,mycomm,ierr)
         end if
#endif
         ncall = ncall + 1
      end if
c
c     initial setup for exchange of bc data between blocks/processors
c
      if (nt.eq.1) then
         iwk1   = 1
         iwk2   = iwk1 + 65*maxbl
         iwk3   = iwk2 + nsub1
         iwk4   = iwk3 + nsub1
         iwk5   = iwk4 + nsub1
         iwk6   = iwk5 + nsub1
         iwork1 = iwork - iwk6
         if (iwork1.lt.0) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''stopping...not enough integer '',
     .                       ''work space for subroutine pre_bc'')')
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         do iii = 1,iwk6
            iwk(iii) = 0
         end do
         call pre_bc(lw,lw2,iwk(iwk1),maxbl,maxgr,maxseg,ninter,intmax,
     .               nsub1,iindex,isav_pat,iwk(iwk2),
     .               iwk(iwk3),iwk(iwk4),iwk(iwk5),mxbli,nbli,limblk,
     .               isva,nblon,nblk,lbcprd,isav_prd,
     .               bcvali,bcvalj,bcvalk,nblg,lbcemb,
     .               iemg,igridg,isav_emb,
     .               iviscg,jdimg,kdimg,idimg,nbci0,nbcj0,nbck0,
     .               nbcjdim,nbckdim,nbcidim,ibcinfo,jbcinfo,
     .               kbcinfo,iadvance,myid,myhost,mycomm,mblk2nd,
     .               nou,bou,nbuf,ibufdim,isav_pat_b,levelg,nblock,
     .               isav_blk)
      end if
c
c***********************************************
c
c     begin multigrid cycles
c
c***********************************************
c
c     loop over multigrid cycles
c
      do 7000 icyc=1,ncyc
c
      if (icyc.gt.1) isklton = 0
c
c     check for user stop file - will then stop on next MG cycle
c     (if unsteady, we don't want to stop until the next time step,
c     so always to all MG cycles)
c
      if (myid .eq. myhost) then
        if (real(dt) .lt. 0.) then
          inquire (file='stop',exist=stop_me)
          if (stop_me) then
c            must set ncyc to current cycle AND set number of cycles on all
c            higher levels to zero
             ncyc = icyc
             do is  = iseq+1,mseq
               ncyc1(is) = 0
             end do
             write(11,*)
             write(11,'(a,i5)') 
     .       'Stopping: user-invoked stop file detected at cycle ',icyc
             write(11,*)
          end if
        end if
      end if
#if defined DIST_MPI
      istop_data(1) = ncyc
      do is =1,mseq
        istop_data(is+1) = ncyc1(is)
      end do
      call MPI_Bcast (istop_data,mseq+1,MPI_INTEGER,myhost,mycomm,ierr)
      ncyc = istop_data(1)
      do is =1,mseq
        ncyc1(is) = istop_data(is+1)
      end do
#endif
c
      if (isklton.eq.1 .and. myid.eq.myhost) then
         write(11,1060)
         write(11,1070)iseq
         write(11,1080)levt
         write(11,1090)levb
         if (nembed.gt.0) write(11,2000)nembed
         lggl = levt-levb-nembed+1
         write(11,2010)lggl
         write(11,2020)lglobal
      end if
c
      do 5070 i=1,level-levb+1
      nsm(i) = ngam
 5070 continue
c
c     the variable "iaccnt" is used to ensure that the residual
c     and force coefficients are accumulated only once per iteration,
c     on the initial restriction portion of the multigrid cycle. This
c     is important when using embedded grids AND the w-cycle in which
c     case the global grid is not the top of the cycle.
c
      iaccnt = 1
c
c     nit = number of solution updates on each multigrid level
c     note: this is not the same as the number of subiterations
c     used for time-accurate solutions
c
      nit = mit(level-levb+1,iseq)
      if (nit.gt.20) then
         write(11,2030)
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
c
c     multigrid cycle starts with downward leg
c
      kxpand = 1
      kode   = 1
      mode   = 0
c
 9000 continue
c
c***************************************************
c     begin/continue downward leg of multigrid cycle
c***************************************************
c
c     find last (ending) block on current level, nblend
c
      do 5090 nbl=1,nblock
      if (level.ne.levelg(nbl)) go to 5090
      nblend = nbl
 5090 continue
c
c     zero out residual and force coefficients for the current
c     multigrid cycle (or time step if time accurate)
c
      if (level.eq.levt) then
         ntt0     = ntt
         if ((real(dt).lt.0. .and. nt.eq.1) .or.
     .       (real(dt).ge.0. .and. icyc.eq.ioutsub)) then
            ntt         = ntt + 1 
            rms(ntt)    = 0.
            clw(ntt)    = 0.
            cdw(ntt)    = 0.
            cyw(ntt)    = 0.
            rmstr(ntt,:) = 0.
            nneg(ntt,:)  = 0
            chdw(ntt)   = 0.
            swetw(ntt)  = 0.
            cmxw(ntt)   = 0.
            cmyw(ntt)   = 0.
            cmzw(ntt)   = 0.
            clcd(:,:,ntt)= 0.
            cxw(ntt)    = 0.
            czw(ntt)    = 0.
            cdpw(ntt)   = 0.
            cdvw(ntt)   = 0.
            fmdotw(ntt) = 0.
            cfttotw(ntt)= 0. 
            cftmomw(ntt)= 0.
            cftpw(ntt)  = 0.
            cftvw(ntt)  = 0.
         end if
         do 6010 ntime=1,nit
         nptsr(ntime)   = 0
         nptsr(ntime+1) = 0
 6010    continue
         nptsrb = 0
         nnegb(:) = 0
         rmstb(:) = 0.
         rmsb   = 0.
         sumn(:)  = 0.
         negn(:)  = 0
      end if
c
c   meshdef = 1 if flow solution is to be bipassed, and mesh deformation only:
c
      if(meshdef.eq.1) goto 6998
c
c     nit1 = number of residual evaluations on each multigrid level;
c     on the downward leg, one extra residual evaluation is performed
c     so that the latest residuals (after solution is updated) from
c     the current level are restricted to the coarser level; the
c     exception is the bottom level - no further residual
c     restrictions are needed, so neither is the extra evaluation
c
      nit1 = nit + min(1,mgflag)
      if (level.eq.levb) nit1 = nit
      if (level.gt.lglobal) nit1 = nit + 1
c
      if (isklton.gt.0 .and. myid.eq.myhost) then
         write(11,2040)level,kode,mode
      end if
      ilc = 0
c
      do 6500 ntime=1,nit1
c
c     initialize density/pressure bcs for comparison in a later call
c     to bcchk, after all bc routines have been called, in order to
c     detect missing bc's
c
      if (isklton.gt.0 .and. ntime.eq.1) then
         do 6990 nbl=1,nblock
            if (iadvance(nbl).lt.0) go to 6990
            if (level.ne.levelg(nbl)) go to 6990
            if (mblk2nd(nbl).eq.myid) then
               call lead(nbl,lw,lw2,maxbl)
               ibcflg = 0
               istop  = 1
               call bcchk(idim,jdim,kdim,w(lq),w(lqi0),w(lqj0),w(lqk0),
     .         w(lblk),ibcflg,nbl,nou,bou,nbuf,ibufdim,myid,istop,
     .         igridg,maxbl)
            end if
 6990    continue
      end if
c
c     update physical boundary conditions
c
      if (iipv.gt.0) then
         if (level.ge.lglobal .and. ntime.eq.nit) then
            nttuse=max(ntt-1,1)
            clwuse = clw(nttuse)
#if defined(DIST_MPI)
            call MPI_Bcast (clwuse,1,MY_MPI_REAL,myhost,
     .                      mycomm,ierr)
#endif
         end if
      end if
      do 6510 nbl=1,nblock
      if (iadvance(nbl).lt.0) go to 6510
      if (level.ne.levelg(nbl)) go to 6510
      if (mblk2nd(nbl).eq.myid) then
         call lead(nbl,lw,lw2,maxbl)
         lres = 1
         nsafe = nwork-lres+1
         call bc(ntime,nbl,lw,lw2,w,mgwk,wk(lres),nsafe,clwuse,
     .           nou,bou,nbuf,ibufdim,maxbl,maxgr,maxseg,itrans,
     .           irotat,idefrm,igridg,nblg,nbci0,nbcj0,nbck0,
     .           nbcidim,nbcjdim,nbckdim,ibcinfo,jbcinfo,
     .           kbcinfo,bcfilei,bcfilej,bcfilek,lwdat,myid,
     .           idimg,jdimg,kdimg,bcfiles,mxbcfil,nummem)
      end if
 6510 continue
c
#if defined(DIST_MPI)
      if (myid.ne.myhost) then
#endif
c
c     update periodic boundary conditions
c
      lres   = 1
      nsafe  = nwork-lres+1
      mneed  = lbcprd*5
      iwk1   = 1
      iwk2   = iwk1 + mneed
      iwk3   = iwk2 + mneed
      iwk4   = iwk3 + mneed
      iwk5   = iwk4 + mneed*2
      iwk6   = iwk5 + mneed
      iwork1 = iwork - iwk6
      if (iwork1.lt.0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),'(''stopping...not enough integer work '',
     .                  ''space for subroutine bc_period'')')
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
      do iii = 1,iwk6
         iwk(iii) = 0
      end do
      call bc_period(ntime,nbl,lw,lw2,w,mgwk,wk(lres),nsafe,maxbl,
     .               maxgr,maxseg,iadvance,bcfilei,bcfilej,bcfilek,
     .               lwdat,xorig,yorig,zorig,jdimg,kdimg,idimg,lbcprd,
     .               isav_prd,period_miss,epsrot,iwk(iwk1),iwk(iwk2),
     .               iwk(iwk3),iwk(iwk4),iwk(iwk5),myid,myhost,
     .               mycomm,mblk2nd,nou,bou,nbuf,ibufdim,
     .               istat2_pe,istat_size,bcfiles,mxbcfil,nummem)
c
c     update embeded-grid boundary conditions
c
      lres   = 1
      nsafe  = nwork-lres+1
      mneed  = lbcemb*3
      iwk1   = 1
      iwk2   = iwk1 + mneed
      iwk3   = iwk2 + mneed
      iwk4   = iwk3 + mneed
      iwk5   = iwk4 + mneed*2
      iwk6   = iwk5 + mneed
      iwork1 = iwork - iwk6
      if (iwork1.lt.0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),'(''stopping...not enough integer work '',
     .                  ''space for subroutine bc_embed'')')
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
      do iii = 1,iwk6
         iwk(iii) = 0
      end do
      call bc_embed(ntime,nbl,lw,lw2,w,mgwk,wk(lres),nsafe,maxbl,
     .              maxgr,lbcemb,iadvance,idimg,jdimg,kdimg,
     .              isav_emb,iwk(iwk1),iwk(iwk2),
     .              iwk(iwk3),iwk(iwk4),iwk(iwk5),myid,myhost,
     .              mycomm,mblk2nd,nou,bou,nbuf,ibufdim,iviscg,
     .              istat2_em,istat_size,nummem)
c
c     update 1-1 block boundary conditions
c
      lres   = 1
      nsafe  = nwork-lres+1
      mneed  = mxbli*5
      iwk1   = 1
      iwk2   = iwk1 + mneed
      iwk3   = iwk2 + mneed
      iwk4   = iwk3 + mneed
      iwk5   = iwk4 + mneed*2
      iwk6   = iwk5 + mneed
      iwork1 = iwork - iwk6
      if (iwork1.lt.0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),'(''stopping...not enough integer work '',
     .                  ''space for subroutine bc_blkint'')')
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
      do iii = 1,iwk6
         iwk(iii) = 0
      end do
      call bc_blkint(ntime,nbl,lw,lw2,w,mgwk,wk(lres),nsafe,maxbl,
     .               maxgr,mxbli,iadvance,geom_miss,epsc0,nblk,nbli,
     .               limblk,isva,nblon,jdimg,kdimg,idimg,
     .               mblk2nd,isav_blk,iwk(iwk1),
     .               iwk(iwk2),iwk(iwk3),iwk(iwk4),iwk(iwk5),
     .               nou,bou,nbuf,ibufdim,myid,myhost,mycomm,
     .               istat2_bl,istat_size,nummem)
c
c     update patch-grid boundary conditions
c
      lres = 1
      nsafe = nwork-lres+1
      mneed  = intmax*nsub1*3
      iwk1   = 1
      iwk2   = iwk1 + mneed
      iwk3   = iwk2 + mneed
      iwk4   = iwk3 + mneed
      iwk5   = iwk4 + mneed*2
      iwk6   = iwk5 + mneed*2
      iwork1 = iwork - iwk6
      if (iwork1.lt.0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),'(''stopping...not enough integer work '',
     .                  ''space for subroutine bc_patch'')')
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
      do iii = 1,iwk6
         iwk(iii) = 0
      end do
      call bc_patch(ntime,nbl,lw,lw2,w,mgwk,wk(lres),nsafe,maxbl,
     .              maxgr,intmax,nsub1,maxxe,iadvance,jdimg,kdimg,
     .              idimg,ninter,windex,iindex,nblkpt,dthetxx,
     .              dthetyy,dthetzz,isav_pat,isav_pat_b,
     .              iwk(iwk1),iwk(iwk2),iwk(iwk3),
     .              iwk(iwk4),iwk(iwk5),myid,myhost,mycomm,
     .              mblk2nd,nou,bou,nbuf,ibufdim,
     .              istat2_pa,istat_size,nummem)
c
#if defined(DIST_MPI)
      end if
#endif
c
c     update chimera boundary conditions
c
      do 6997 nbl=1,nblock
      if (iadvance(nbl).lt.0) go to 6997
      if (level.ne.levelg(nbl)) go to 6997
      if (mblk2nd(nbl).eq.myid) then
         call lead(nbl,lw,lw2,maxbl)
         lres = 1
         nsafe = nwork-lres+1
         call bc_xmera(ntime,nbl,lw,lw2,w,mgwk,wk(lres),nsafe,maxbl,
     .                 iitot,iviscg,iovrlp,lbg,ibpntsg,qb,iibg,kkbg,
     .                 jjbg,ibcg,nou,bou,nbuf,ibufdim,1,nummem)
      end if
 6997 continue
c
c     output bc info to main output file and, the first time through:
c     fill endpoints for safety with multi-plane vectorization (ibcflg=1)
c     check to make sure all boundary conditions have been set (ibcflg=2)
c
      do 6988 nbl=1,nblock
      if (level.ne.levelg(nbl)) go to 6988
      if (mblk2nd(nbl).eq.myid) then
         call lead(nbl,lw,lw2,maxbl)
         lres = 1
         nsafe = nwork-lres+1
         call bc_info(ntime,nbl,w,mgwk,ibcinfo,jbcinfo,
     .                kbcinfo,nbci0,nbcj0,nbck0,nbcidim,nbcjdim,
     .                nbckdim,bcfilei,bcfilej,bcfilek,igridg,itrans,
     .                irotat,idefrm,idimg,jdimg,kdimg,nblg,iibg,
     .                kkbg,jjbg,ibcg,ibpntsg,iipntsg,lig,lbg,
     .                isav_blk,isav_prd,iemg,nbli,
     .                geom_miss,period_miss,ninter,iindex,nou,bou,nbuf,
     .                myid,mblk2nd,ibufdim,maxbl,maxseg,iitot,
     .                intmax,nsub1,mxbli,lbcprd,epsc0,epsrot,lwdat,
     .                maxgr,iovrlp,bcfiles,mxbcfil)
      end if
c
      if (isklton.gt.0) then
#   ifdef FASTIO
         call writ_buffast(nbl,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .                 mycomm,mblk2nd,maxbl,02)
#   else
         call writ_buf(nbl,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .                 mycomm,mblk2nd,maxbl)
#   endif
      end if
c
 6988 continue
c
      do 6989 nbl=1,nblock
      if (iadvance(nbl).lt.0) go to 6989
      if (level.ne.levelg(nbl)) go to 6989
      if (mblk2nd(nbl).eq.myid) then
         call lead(nbl,lw,lw2,maxbl)
         if (isklton.gt.0 .and. ntime.eq.1) then
            ibcflg = 1
            istop  = 1
            call bcchk(idim,jdim,kdim,w(lq),w(lqi0),w(lqj0),w(lqk0),
     .      w(lblk),ibcflg,nbl,nou,bou,nbuf,ibufdim,myid,istop,
     .      igridg,maxbl)
            ibcflg = 2
            istop  = 1
            if (ichk.eq.-1 .or. ichk.gt.1) istop = 0
            call bcchk(idim,jdim,kdim,w(lq),w(lqi0),w(lqj0),w(lqk0),
     .      w(lblk),ibcflg,nbl,nou,bou,nbuf,ibufdim,myid,istop,
     .      igridg,maxbl)
            if (idefrm(nbl).eq.1) then
               call chkdef(nbl,idim,jdim,kdim,w(lbci),w(lbcj),w(lbck),
     .                     icsi,icsf,jcsi,jcsf,kcsi,kcsf,nsegdfrm,
     .                     idfrmseg,maxbl,maxsegdg,nou,bou,nbuf,
     .                     ibufdim,myid)
            end if
         end if
      end if
      if (isklton.gt.0) then
#   ifdef FASTIO
         call writ_buffast(nbl,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .                 mycomm,mblk2nd,maxbl,03)
#   else
         call writ_buf(nbl,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .                 mycomm,mblk2nd,maxbl)
#   endif
      end if
 6989 continue
c
c     cycle through blocks
c
#if defined DIST_MPI
      nreq = 1
c
#endif
      if (iteravg .eq. 1 .or. iteravg .eq. 2) then
        if (level .ge. lglobal .and. ntime .eq. nit) then
          if (real(dt) .lt. 0. .or. (real(dt) .gt. 0. and. 
     +        icyc .eq. ncyc)) then
            xnumavg=xnumavg+1.
            xnumavg2=xnumavg2+1.
          end if
        end if
      end if
      do 8000 nbl=1,nblock
c
      if (level.ne.levelg(nbl)) go to 8000
c
      call lead(nbl,lw,lw2,maxbl)
c
#if defined DIST_MPI
c
c     ramp time step/cfl number on host processor
c
      if (myid.eq.myhost) then
         if (kode.eq.-1 .or. ntime.gt.1) go to 130
c
c        local time step: ramp initial CFL number to a value fmax times larger.
c        The time step distribution is recalculated for each cycle.
c
         dttol = 1.e-6
         if (iflagts.gt.0 .and. real(dt).lt.0 .and. icyc.gt.1) then 
           if (nt.eq.1 .and. abs(real(dt)) .lt. 
     .         real(fmax)*abs(real(dt0))-real(dttol)) then
c             ramp CFL on first block on top level
              if (level.eq.levt .and. nbl.eq.nblstat) then
                 t2   = 1.e0/float(iflagts)
                 fact = fmax**t2
                 dt   = fact*dt
c                nou(1) = min(nou(1)+1,ibufdim)
c                write(bou(nou(1),1),*)'cycle, cfl = ',icyc,real(dt)
              end if
           end if
         end if
c
c        global time step:ramp initial time step to a value fmax times larger.
c        The CFL number distribution is recalculated for each time step.
c
         dttol = 1.e-6
         if (iflagts.gt.0 .and. real(dt).gt.0 .and. nt.gt.1) then
            if (icyc.eq.1 .and. abs(real(dt)) .lt. 
     .         real(fmax)*abs(real(dt0))-real(dttol)) then
c              ramp time step on first block on top level
               if (level.eq.levt .and. nbl.eq.nblstat) then
                 t2   = 1.e0/float(iflagts)
                 fact = fmax**t2
                 dt   = fact*dt
c                nou(1) = min(nou(1)+1,ibufdim)
c                write(bou(nou(1),1),*)'time step, dt = ',nt,real(dt)
              end if
           end if
         end if
 130     continue
      end if
#endif
c
      if (mblk2nd(nbl) .eq. myid) then
c
c        for time accurate calculations, compute and store
c        dqc0 and qc0, where qc0 is the solution at the start
c        of the subiterations (i.e. the solution at the last
c        time step) and dqc0 is the change in solution between
c        the last time step and the next-to-last time step (needed
c        for temporally second-order accurate solutions).
c        also, if the grid is undergoing rotation, rotate initial 
c        solution (i.e. the solution from the last time step) to
c        correspond to the new grid position. this helps convergence
c        of the subiterations, particularly if large timesteps are used.
c
c        NOTE: this prerotation will not preserve freestream unless a
c        very large number of subiterations are done. if performing
c        freestream preservation tests, set iprerot = 0 to skip the
c        prerotation step for rotaing grids.
c
         iprerot = 1
c
         if (real(dt).gt.0.0 .and. iadvance(nbl).ge.0) then
            if (icyc.eq.1 .and. level.ge.lglobal .and. ntime.eq.1) then
               if (abs(ita).eq.2) then
                  call setdqc0(jdim,kdim,idim,w(lq),w(lqc0),w(ldqc0))
               end if
               call setqc0(jdim,kdim,idim,w(lq),w(lqc0))
               if (irotat(nbl).gt.0 .and. iprerot.gt.0) then
                  dthtx  = thetax(nbl) - thetaxl(nbl)
                  dthty  = thetay(nbl) - thetayl(nbl)
                  dthtz  = thetaz(nbl) - thetazl(nbl)
                  lqwk = 1
                  call rotateq(jdim,kdim,idim,w(lq),wk(lqwk),1,idim,
     .                         1,jdim,1,kdim,dthtx,dthty,dthtz)
                  do 40 l=1,jdim*kdim*idim*5
                  w(l+lq-1) = wk(l+lqwk-1)
   40             continue
               end if
            end if
         end if
c
c        accumulate fluxes at embedded grid boundaries
c
         ifluxa = 0
         if (level.ge.lglobal .and. level.ne.levt)ifluxa = min(1,iconsf)
         if (ifluxa.gt.0 .and. iadvance(nbl).ge.0) then
            ifamax = maxbl*7*3
            call fa(nbl,lw,w,mgwk,wk,nwork,iw,ifamax,nwfa,nifa,nfajki,
     .              maxbl,maxseg,jdimg,kdimg,idimg,jsg,ksg,isg,jeg,
     .              keg,ieg,jbcinfo,kbcinfo,ibcinfo,nblock,nblcg,
     .              nou,bou,nbuf,ibufdim)
         else
            nwfa      = 0
            nfajki(1) = 0
            nfajki(2) = 0
            nfajki(3) = 0
         end if
c
         lres      = nwfa   + 1
         lwj0      = lres   + jdim*kdim*(idim-1)*5
         lwk0      = lwj0   + kdim*idim*22
         lwi0      = lwk0   + jdim*idim*22
         lvmuk     = lwi0   + kdim*jdim*22
         lvmuj     = lvmuk  + 2*(jdim-1)*(idim-1)
         lvmui     = lvmuj  + 2*(kdim-1)*(idim-1)
         ltot      = lvmui  + 2*(kdim-1)*(jdim-1)
         lsafe     = nwork - ltot
         if (lsafe.lt.0) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2050) lsafe
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
c        compute residual = r(q)
c          residual is comprised of: 
c          1) flux contributions
c          2) subiteration terms and 2nd order time terms 
c             (qc0 and dqc0; time accurate only)
c          3) correction from finer level 
c             (coarser levels only)
c
c        flux contributions
c
         isf = 0
         if (level.gt.lglobal .and. ntime.eq.nit1) isf = min(1,iconsf)
         call resid(nbl,ntime,jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),
     .              w(lqi0),w(lsj),w(lsk),w(lsi),w(lvol),w(ldtj),w(lx),
     .              w(ly),w(lz),w(lvis),w(lsni0),w(lsnk0),w(lsni0),
     .              wk(lres),wk(lwj0),wk(lwk0),wk(lwi0),wk(lvmuk),
     .              wk(lvmuj),wk(lvmui),wk(ltot),lsafe,
     .              isf,iw,wk,deltat(level),w(lblk),iovrlp(nbl),
     .              nblendg,nblstat,nblstag,w(lxib),w(lsig),
     .              w(lsqtq),w(lg),w(ltj0),w(ltk0),w(lti0),
     .              w(lxkb),w(lnbl),w(lvj0),w(lvk0),w(lvi0),w(lbcj),
     .              w(lbck),w(lbci),nt,sumn,negn,w(lux),
     .              w(lxib2),w(lcmuv),w(lvolj0),w(lvolk0),w(lvoli0),
     .              nou,bou,nbuf,ibufdim,myid,mblk2nd,maxbl,maxseg,
     .              nbci0,nbcj0,nbck0,nbcidim,nbcjdim,nbckdim,ibcinfo,
     .              jbcinfo,kbcinfo,vormax,ivmax,jvmax,kvmax,idefrm,
     .              iadvance,w(lqavg),nummem)
c
c        add 2nd order time terms and subiteration terms
c
         if (real(dt).gt.0.e0 .and. iadvance(nbl).ge.0) then
            call resadd(jdim,kdim,idim,w(lq),w(lqc0),w(ldqc0),
     .                  wk(lres),w(lvol),iovrlp(nbl),w(lblk))
         end if
c
      end if
c
      if (isklton.gt.0) then
#   ifdef FASTIO
         call writ_buffast(nbl,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .                 mycomm,mblk2nd,maxbl,04)
#   else
         call writ_buf(nbl,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .                 mycomm,mblk2nd,maxbl)
#   endif
         if (ivmx.eq.2 .or. ivmx.eq.3) then
#   ifdef FASTIO
            call writ_buffast(nbl,14,nou,bou,nbuf,ibufdim,myhost,myid,
     .                    mycomm,mblk2nd,maxbl,05)
#   else
            call writ_buf(nbl,14,nou,bou,nbuf,ibufdim,myhost,myid,
     .                    mycomm,mblk2nd,maxbl)
#   endif
         end if
      else
c
c        insure ramped time step/cfl number info gets written 
c        to the output file in the correct order
c
         if (iflagts.gt.0 .and. abs(real(fmax)).gt.0. .and.
     .      real(dt).lt.0.) then
            if (((icyc.eq.iflagts+1) .or.
     .         (icyc.eq.ncyc.and.iflagts+1.gt.ncyc))  .and.
     .         level.eq.levt.and.nbl.eq.nblstat_h.and.nt.eq.1) then
#   ifdef FASTIO
               call writ_buffast(nbl,11,nou,bou,nbuf,ibufdim,myhost,
     .                       myid,mycomm,mblk2nd,maxbl,06)
#   else
               call writ_buf(nbl,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .                       mycomm,mblk2nd,maxbl)
#   endif
            end if
         end if
         if (iflagts.gt.0 .and. abs(real(fmax)) .gt. 0. .and.
     .      real(dt).gt.0.) then
            if (((nt.eq.iflagts+1) .or.
     .         (nt.eq.ntstep.and.iflagts+1.gt.ntstep)) .and.
     .         level.eq.levt.and.nbl.eq.nblstat_h.and.icyc.eq.1) then
#   ifdef FASTIO
               call writ_buffast(nbl,11,nou,bou,nbuf,ibufdim,myhost,
     .                       myid,mycomm,mblk2nd,maxbl,07)
#   else
               call writ_buf(nbl,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .                       mycomm,mblk2nd,maxbl)
#   endif
            end if
         end if
      end if
c
c     add correction from finer levels
c
      if (mgflag.gt.0) then
         if (level.ne.levt) then
            if (mblk2nd(nbl) .eq. myid) then
               if (iadvance(nbl).ge.0) then
               call tau(mgflag,nbl,jdim,kdim,idim,w(lq),wk(lres),
     .                  w(lq1),w(lqr),lw,w,nou,bou,nbuf,ibufdim,maxbl,
     .                  maxgr,nblock,igridg,nblcg,jsg,ksg,isg,
     .                  jeg,keg,ieg,iemg)
               end if
            end if
            if (nbl.eq.nblend .and. kode.eq.2) kode = 1
         end if
      end if
c
      if (isklton.gt.0) then
#   ifdef FASTIO
         call writ_buffast(nbl,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .                 mycomm,mblk2nd,maxbl,08)
#   else
         call writ_buf(nbl,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .                 mycomm,mblk2nd,maxbl)
#   endif
      end if
c
c     compute force coefficients and L2 norm of residual
c
      if (level.ge.lglobal .and. iaccnt.gt.0 .and. ntime.eq.nit) then
         if (myid.eq.mblk2nd(nbl) .or. myid.eq.myhost) then
#if defined DIST_MPI
            if (myid.eq.myhost) then
c              the following were defined earlier for the nodes,
c              but not for the host
               nwfa      = 0
               nfajki(1) = 0
               nfajki(2) = 0
               nfajki(3) = 0
               lres      = nwfa   + 1
               lwj0      = lres   + jdim*kdim*(idim-1)*5
               lwk0      = lwj0   + kdim*idim*22
               lwi0      = lwk0   + jdim*idim*22
               lvmuk     = lwi0   + kdim*jdim*22
               lvmuj     = lvmuk  + 2*(jdim-1)*(idim-1)
               lvmui     = lvmuj  + 2*(kdim-1)*(idim-1)
               ltot      = lvmui  + 2*(kdim-1)*(jdim-1)
               lsafe     = nwork - ltot
               if (lsafe.lt.0) then
                  write(11,2050) lsafe
                  call termn8(myid,-1,ibufdim,nbuf,bou,nou)
               end if
             end if
#endif
             call resp(nbl,ntime,ntt0,jdim,kdim,idim,w(lq),w(lsj),
     .                 w(lsk),w(lsi),w(lvol),w(lx),w(ly),w(lz),
     .                 wk(lres),wk(lwk0),wk(lvmuk),wk(lvmuj),
     .                 wk(lvmui),wk(ltot),lsafe,nblstag,nblendg,
     .                 nblstat,nblendt,nptsr(ntime),nptsrb,
     .                 nnegb,w(lbcj),w(lbck),w(lbci),
     .                 w(lblk),nt,sumn,negn,rmst,
     .                 w(lxtbj),w(lxtbk),w(lxtbi),w(lqc0),w(ldqc0),
     .                 iseq,w(lqj0),w(lqk0),w(lqi0),maxbl,maxgr,
     .                 maxseg,rms,rmsb,rmstb,clw,cdw,cdpw,cdvw,
     .                 cxw,cyw,czw,cmxw,cmyw,cmzw,
     $                 n_clcd,clcd,nblocks_clcd,blocks_clcd,
     $                 chdw,swetw,fmdotw,
     .                 cfttotw,cftmomw,cftpw,cftvw,rmstr,
     .                 nneg,ihstry,ngrid,nblg,iemg,levelg,
     .                 iviscg,itrans,irotat,iforce,swett,clt,cdt,
     .                 cxt,cyt,czt,cmxt,cmyt,cmzt,cdpt,cdvt,nbci0,
     .                 nbcj0,nbck0,nbcidim,nbcjdim,nbckdim,ibcinfo,
     .                 jbcinfo,kbcinfo,resmx,imx,jmx,kmx,vormax,
     .                 ivmax,jvmax,kvmax,sx,sy,sz,stot,pav,ptav,
     .                 tav,ttav,xmav,fmdot,cfxp,cfyp,cfzp,cfdp,
     .                 cflp,cftp,cfxv,cfyv,cfzv,cfdv,cflv,cftv,
     .                 cfxmom,cfymom,cfzmom,cfdmom,cflmom,cftmom,
     .                 cfxtot,cfytot,cfztot,cfdtot,cfltot,cfttot,
     .                 ncs,icsinfo,myid,myhost,mycomm,mblk2nd,nou,
     .                 bou,nbuf,ibufdim,ncycmax,maxcs,thetay,iadvance,
     .                 idefrm,igridg,nummem)
         end if
c
c        after last block on global level is accounted for, set
c        iaccnt = 0 to insure forces are not accumulated again
c        in this multigrid cycle
c
         if(nbl.eq.nblendg) iaccnt = 0
c
      end if
c
      if (iadvance(nbl).lt.0) go to 7999
c
c     smooth residual, if desired
c
      if (mblk2nd(nbl) .eq. myid) then
         if (issr.gt.0) then
            icall = idim-1
            call rsmooth(epsssr,idim,jdim,kdim,icall,wk(lres),wk(ltot),
     .                   nou,bou,nbuf,ibufdim)
         end if
      end if
c
c     update solution
c
      if (ntime.le.nit ) then
c
         if (mblk2nd(nbl) .eq. myid) then
            call update(jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),w(lqi0),
     .                  w(lsj),w(lsk),w(lsi),w(lvol),w(ldtj),w(lvis),
     .                  w(lblk),w(lx),w(ly),w(lz), wk(lres),wk(lwk0),
     .                  wk(lvmuk),wk(lvmuj),wk(lvmui),wk(ltot),lsafe,
     .                  nbl,iovrlp(nbl),w(lvk0),w(lbcj),w(lbck),w(lbci),
     .                  nou,bou,nbuf,ibufdim,myid,mblk2nd,maxbl,
     .                  w(lvolk0),w(lxib),w(ltk0),w(lcmuv),
     .                  iadvance,nummem,w(lux))
c 
c           update the overlapped values if chimera scheme is used
c
            if (iovrlp(nbl).eq.1) then
                ldim = 5
                if (isklton.eq.1) then
                   nou(1) = min(nou(1)+1,ibufdim)
                   write(bou(nou(1),1),2060) nbl
                end if
                need  = (idim+1)*(jdim+1)*(kdim+1)*5
                if (need .gt. nwork) then
                   nou(1) = min(nou(1)+1,ibufdim)
                   write(bou(nou(1),1),'(''stopping...not enough work'',
     .             '' space for subroutine intrbc (q)'')')
                   nou(1) = min(nou(1)+1,ibufdim)
                   write(bou(nou(1),1),'(''need, have '',2i11)')
     .             need,nwork
                   call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                end if
                call intrbc(w(lq),jdim,kdim,idim,nbl,ldim,maxbl,
     .                      iitot,lig,iipntsg,dxintg,dyintg,dzintg,
     .                      iiig,jjig,kkig,qb,w(lqj0),w(lqk0),
     .                      w(lqi0),wk,w(lbcj),w(lbck),w(lbci),nou,
     .                      bou,nbuf,ibufdim,1,1)
c               turbulence quantities
                if(iviscg(nbl,1).ge.2 .or. iviscg(nbl,2).ge.2 .or.
     .             iviscg(nbl,3).ge.2) then
                   ldim = 1
                   call intrbc(w(lvis),jdim,kdim,idim,nbl,ldim,maxbl,
     .                         iitot,lig,iipntsg,dxintg,dyintg,dzintg,
     .                         iiig,jjig,kkig,qb,w(lvj0),w(lvk0),
     .                         w(lvi0),wk,w(lbcj),w(lbck),w(lbci),nou,
     .                         bou,nbuf,ibufdim,1,2)
                end if
                if(iviscg(nbl,1).ge.4 .or. iviscg(nbl,2).ge.4 .or.
     .             iviscg(nbl,3).ge.4) then
                   ldim = nummem
                   call intrbc(w(lxib),jdim,kdim,idim,nbl,ldim,maxbl,
     .                         iitot,lig,iipntsg,dxintg,dyintg,dzintg,
     .                         iiig,jjig,kkig,qb,w(ltj0),w(ltk0),
     .                         w(lti0),wk,w(lbcj),w(lbck),w(lbci),nou,
     .                         bou,nbuf,ibufdim,1,3)
                end if
            end if
c
         end if
c
#if defined DIST_MPI
c        share qb information with the processors
c
         if (iovrlp(nbl).eq.1) then
c
            nvalqb  = iipntsg(nbl)
c
            if (myid.ne.mblk2nd(nbl) .and. myid.ne.myhost) then

c              recieve the data on non-local nodes
c
               nd_srce = mblk2nd(nbl)
c              receive q data 
               do ll=1,5
                  mytag  = mytag_qb(ll) + nbl
                  call MPI_Irecv (qb(lig(nbl),ll,1), nvalqb,
     .                            MY_MPI_REAL,
     .                            nd_srce, mytag, mycomm,
     .                            ireq_qb(nreq), ierr)
                  nreq = nreq + 1
               end do
               if (ivmx.ge.2) then
c                 receive vist3d data
                  mytag = mytag_qb(6) + nbl
                  call MPI_Irecv (qb(lig(nbl),1,2), nvalqb,
     .                            MY_MPI_REAL,
     .                            nd_srce, mytag, mycomm,
     .                           ireq_qb(nreq), ierr)
                  nreq = nreq + 1
               end if
               if (ivmx.ge.4) then
c                 receive turb. data            
                    do ll=1,nummem
                       mytag = mytag_qb(6+ll) + nbl
                       call MPI_Irecv (qb(lig(nbl),ll,3), nvalqb,
     .                                 MY_MPI_REAL,
     .                                 nd_srce, mytag, mycomm,
     .                                 ireq_qb(nreq), ierr)
                       nreq = nreq + 1
                    end do
               end if
c
            else if (myid.eq.mblk2nd(nbl)) then
c
c              send the data to non-local nodes
c
               do inode = 1,nnodes
                  if (inode .ne. myid) then
                     nd_dest = inode
c                    send q data
                     do ll=1,5
                        mytag  = mytag_qb(ll) + nbl
                        call MPI_Send (qb(lig(nbl),ll,1), nvalqb,
     .                                 MY_MPI_REAL,
     .                                 nd_dest, mytag,
     .                                 mycomm, ierr)
                     end do
                     if (ivmx.ge.2) then
c                       send vist3d data
                        mytag = mytag_qb(6) + nbl
                        call MPI_Send (qb(lig(nbl),1,2), nvalqb,
     .                                 MY_MPI_REAL,
     .                                 nd_dest, mytag,
     .                                 mycomm, ierr)
                     end if
                     if (ivmx.ge.4) then
c                       send turb. data
                        do ll=1,nummem
                           mytag = mytag_qb(6+ll) + nbl
                           call MPI_Send (qb(lig(nbl),ll,3), nvalqb,
     .                                    MY_MPI_REAL,
     .                                    nd_dest, mytag,
     .                                    mycomm, ierr)
                        end do
                     end if
                  end if
               end do
c
            end if
c
         end if
c
#endif
      end if
c
c     restrict q and r to coarser grid levels
c
      if (mgflag.gt.0 .and. ntime.eq.nit1 .and. level.ne.levb) then
c
         if (mblk2nd(nbl) .eq. myid) then
c
            if (level.le.lglobal) then
c
c              global grids
c
               lqc   = lw( 1,nbl+1)
               lvolc = lw( 8,nbl+1)
               lvtc  = lw(13,nbl+1)
               lqrc  = lw(17,nbl+1)
               lxibc = lw(19,nbl+1)
               lqc0c  = lw(34,nbl+1)
               ldqc0c = lw(35,nbl+1)
               nroom = nwork-lwj0
               mdim  = jdim*kdim*idim*5
               if (nroom.lt.mdim) then
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),2070)
                  call termn8(myid,-1,ibufdim,nbuf,bou,nou)
               end if
               call collq(w(lq),w(lqc),w(lvol),w(lvolc),jdim,kdim,idim,
     .                    jj2,kk2,ii2,wk(lres),w(lqrc),wk(lwj0),
     .                    w(lvis),w(lvtc),w(lxib),w(lxibc),
     .                    nbl,nou,bou,nbuf,ibufdim,nummem)
c
c              restrict qc0 and dqc0 for time-accurate multigrid
c
               if (real(dt).gt.0.e0) then
               call collqc0(w(lqc0),w(lqc0c),w(lvol),w(lvolc),
     .                      jdim,kdim,idim,jj2,kk2,ii2,
     .                      w(ldqc0),w(ldqc0c),wk(lwj0),nbl,
     .                      nou,bou,nbuf,ibufdim)
               end if
c
            else
c
c              embedded grids
c
               jc    = jdimg(nblc)
               kc    = kdimg(nblc)
               ic    = idimg(nblc)
               lqc   = lw( 1,nblc)
               lvolc = lw( 8,nblc)
               lqrc  = lw(17,nbl)
               lvtc  = lw(13,nblc)
               lxtc  = lw(19,nblc)
               nroom = nwork-lwj0
               mdim  = jdim*kdim*idim*5
               if (nroom.lt.mdim) then
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),2080)
                  call termn8(myid,-1,ibufdim,nbuf,bou,nou)
               end if
               call coll2q(w(lq),w(lqc),w(lvol),w(lvolc),jdim,kdim,idim,
     .                     jc,kc,ic,wk(lres),w(lqrc),wk(lwj0),
     .                     js,ks,is,je,ke,ie,nbl,nblc,
     .                     w(lvis),w(lvtc),w(lxib),w(lxtc),
     .                     nou,bou,nbuf,ibufdim,nummem)
            end if
c
         end if
c
#if defined DIST_MPI
         if (level.gt.lglobal) then
c
            nd_dest = mblk2nd(nblc)
            nd_srce = mblk2nd(nbl)
c
            if (nd_dest .ne. nd_srce) then
c
c              embedded and global blocks reside on different
c              processors, so need to pass qc data
c
               jc    = jdimg(nblc)
               kc    = kdimg(nblc)
               ic    = idimg(nblc)
               lqc   = lw( 1,nblc)
               lvolc = lw( 8,nblc)
               lqrc  = lw(17,nbl)
               lvtc  = lw(13,nblc)
               lxtc  = lw(19,nblc)
c
               mytag = itag_qc + nbl
               jki5  =  (je-js+1)*(ke-ks+1)*(ie-is+1)*5
               if (myid.eq.nd_srce) then
                  call ld_qc(w(lqc),wk,jc,kc,ic,is,ie,js,je,
     .                       ks,ke,5)
                  call MPI_Send(wk,jki5,MY_MPI_REAL,
     .                          nd_dest,mytag,mycomm,ierr)
               else if (myid .eq. nd_dest) then
                  call MPI_Recv(wk,jki5,MY_MPI_REAL,
     .                          nd_srce,mytag,mycomm,istat,ierr)
                  call unld_qc(w(lqc),wk,jc,kc,ic,is,ie,js,je,
     .                         ks,ke,5)
               end if
               if (ivmx.ge.2) then
                  mytag = itag_qv + nbl
                  jki  =  (je-js+1)*(ke-ks+1)*(ie-is+1)
                  if (myid.eq.nd_srce) then
                     call ld_qc(w(lvtc),wk,jc,kc,ic,is,ie,js,je,
     .                          ks,ke,1)
                     call MPI_Send(wk,jki,MY_MPI_REAL,
     .                             nd_dest,mytag,mycomm,ierr)
                  else if (myid .eq. nd_dest) then
                     call MPI_Recv(wk,jki,MY_MPI_REAL,
     .                             nd_srce,mytag,mycomm,istat,ierr)
                     call unld_qc(w(lvtc),wk,jc,kc,ic,is,ie,js,je,
     .                            ks,ke,1)
                  end if
               end if
               if (ivmx.ge.4) then
                  mytag = itag_qt + nbl
                  jki2  =  (je-js+1)*(ke-ks+1)*(ie-is+1)*nummem
                  if (myid.eq.nd_srce) then
                     call ld_qc(w(lxtc),wk,jc,kc,ic,is,ie,js,je,
     .                          ks,ke,nummem)
                     call MPI_Send(wk,jki2,MY_MPI_REAL,
     .                             nd_dest,mytag,mycomm,ierr)
                  else if (myid .eq. nd_dest) then
                     call MPI_Recv(wk,jki2,MY_MPI_REAL,
     .                             nd_srce,mytag,mycomm,istat,ierr)
                     call unld_qc(w(lxtc),wk,jc,kc,ic,is,ie,js,je,
     .                            ks,ke,nummem)
                  end if
               end if
               mytag = itag_qr + nbl
               jkim5 = (je-js+1)*(ke-ks+1)*(ie-is)*5
               if (myid.eq.nd_srce) then
                  call MPI_Send(w(lqrc),jkim5,MY_MPI_REAL,
     .                          nd_dest,mytag,mycomm,ierr)
               else if (myid .eq. nd_dest) then
                  call MPI_Recv(w(lqrc),jkim5,MY_MPI_REAL,
     .                          nd_srce,mytag,mycomm,istat,ierr)
               end if
c
            end if
c
         end if
#endif
      end if
c
 7999 continue
c
      if (isklton.gt.0) then
#   ifdef FASTIO
         call writ_buffast(nbl,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .                 mycomm,mblk2nd,maxbl,09)
#   else
         call writ_buf(nbl,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .                 mycomm,mblk2nd,maxbl)
#   endif
      end if
c
 8000 continue
c
#if defined DIST_MPI
c
c     wait for chimera data here
c
      if (nreq-1 .gt. 0) then
         do nnn=1,nreq-1
            call MPI_Wait(ireq_qb(nnn), istat, ierr)
         end do
      end if
c
#endif
 6500 continue
c     if using fixed Cl option, update angle of attack, then perform
c     instantaneous velocity rotation for all points in flowfield;
c     rotation angle is change in alpha, velocity magnitude remains unchanged
c
      if (icycupdt .gt. 0) then
      if (mod(ntt,icycupdt) .eq. 0) then
      if (ialphit .ne. 0 .and. level.eq.levt) then
         call newalpha(ncycmax,rms,clw,myid,myhost,mycomm)
         do nbl = 1, nblock
            if (mblk2nd(nbl) .eq. myid) then
               call lead(nbl,lw,lw2,maxbl)
               angx = 0
               angy = -dalpha
               angz = 0
               call rotateq(jdim,kdim,idim,w(lq),w(lq),1,idim,
     .                      1,jdim,1,kdim,angx,angy,angz)
            end if
         end do
      end if
      end if
      end if
c
      if (ilc.eq.0) then
         nsm(level-levb+1) = nsm(level-levb+1)-1
         if (nsm(level-levb+1).lt.0) nsm(level-levb+1) = ngam-1
         ilc = 1
      end if
c
      if (isklton.gt.0 .and. myid.eq.myhost)write(11,2090)
c
      kode = 1
c
      if (kxpand.eq.-1) then
c        finished with downward leg; begin upward leg
         go to 7500
      end if
c
      if (level.eq.levb .and. levb.eq.levt) then
c        single level case; start new cycle
         if(isklton.eq.1 .and. myid.eq.myhost) write(11,3000)
         go to 6999
      end if
c
      kode   = 2
      mode   = 1
      level  = level - 1
      ntime  = 0     
      nit    = mit(level-levb+1,iseq)
c
      if (level.gt.levb) then
c        not yet at bottom of cycle; continue downward leg
         go to 9000 
      else
c        at bottom of cycle; complete downward leg and set
c        kxpand so that upward leg is done after last level
c        on downward leg has been completed
         nit    = mit(1,iseq) + mtt
         kxpand = -1
         go to 9000
      end if
c
c*************************************************
c     end downward leg of multigrid cycle
c*************************************************
c
 7500 continue
c
c*************************************************
c     begin/continue upward leg of multigrid cycle
c*************************************************
c
      if (mgflag.ne.0) then
c
         if (myid.eq.myhost) then
            if (level.lt.lglobal) then
               if (isklton.gt.0)write(11,3010)level,kode,mode
            else if (mgflag.gt.1) then
               if (isklton.gt.0)write(11,3010)level,kode,mode
            end if
         end if
c
#if defined DIST_MPI
         nreq = 1
c
#endif
         do 7020 nbl=1,nblock
c
         if (iadvance(nbl).lt.0) go to 7020
         if (level.ne.levelg(nbl)) go to 7020
         lqc  = lw( 1,nbl)
         lq1c = lw(16,nbl)
         if (level.lt.lglobal) then
c
c           add corrections to global grids
c
            nblf = nbl-1
c
            if (mblk2nd(nbl) .eq. myid) then
c
               call lead(nblf,lw,lw2,maxbl)
               iwk1 = 1
               iwk2 = iwk1 + jdim*kdim*idim*5
               iwk3 = iwk2 + jj2*kk2*ii2*5
               nroom = nwork-iwk3
               mdim  = jdim*kk2*ii2*5
               if (nroom.lt.mdim) then
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),3020)
                  call termn8(myid,-1,ibufdim,nbuf,bou,nou)
               end if
c
              if (isklton.gt.0) then
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),3030)nbl,nblf,igridg(nblf)
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),3040)jdim,kdim,idim
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),3050)jj2,kk2,ii2
               end if
c
               call addx(w(lq),w(lqc),jdim,kdim,idim,jj2,kk2,ii2,
     .                   w(lq1c),wk(iwk1),wk(iwk2),wk(iwk3),nbl,w(lblk),
     .                   nou,bou,nbuf,ibufdim,5,myid)
c
c              add corrections to chimera grids
c
               if (iovrlp(nblf).eq.1) then
                  if (isklton.eq.1) then
                     nou(1) = min(nou(1)+1,ibufdim)
                     write(bou(nou(1),1),2060) nblf
                  end if
                  ldim=5
                  call intrbc(w(lq),jdim,kdim,idim,nblf,ldim,maxbl,
     .                        iitot,lig,iipntsg,dxintg,dyintg,dzintg,
     .                        iiig,jjig,kkig,qb,w(lqj0),w(lqk0),
     .                        w(lqi0),wk,w(lbcj),w(lbck),w(lbci),nou,
     .                        bou,nbuf,ibufdim,1,1)
c                 turbulence quantities
                  if(iviscg(nblf,1).ge.2 .or. iviscg(nblf,2).ge.2 .or.
     .               iviscg(nblf,3).ge.2) then
                     ldim = 1
                     call intrbc(w(lvis),jdim,kdim,idim,nblf,ldim,maxbl,
     .                           iitot,lig,iipntsg,dxintg,dyintg,dzintg,
     .                           iiig,jjig,kkig,qb,w(lvj0),w(lvk0),
     .                           w(lvi0),wk,w(lbcj),w(lbck),w(lbci),nou,
     .                           bou,nbuf,ibufdim,1,2)
                  end if
                  if(iviscg(nblf,1).ge.4 .or. iviscg(nblf,2).ge.4 .or.
     .               iviscg(nblf,3).ge.4) then
                     ldim = nummem
                     call intrbc(w(lxib),jdim,kdim,idim,nblf,ldim,maxbl,
     .                           iitot,lig,iipntsg,dxintg,dyintg,dzintg,
     .                           iiig,jjig,kkig,qb,w(ltj0),w(ltk0),
     .                           w(lti0),wk,w(lbcj),w(lbck),w(lbci),nou,
     .                           bou,nbuf,ibufdim,1,3)
                  end if
               end if
c
            end if
c
#if defined DIST_MPI
c
c           share corrected qb information with the processors
c
            if (iovrlp(nblf).eq.1) then
c
               nvalqb  = iipntsg(nblf)
c
               if (myid.ne.mblk2nd(nblf) .and. myid.ne.myhost) then

c                 recieve the data on non-local nodes
c
                  nd_srce = mblk2nd(nblf)
c                 receive q data 
                  do ll=1,5
                     mytag  = mytag_qb(ll) + nblf
                     call MPI_Irecv (qb(lig(nblf),ll,1), nvalqb,
     .                               MY_MPI_REAL,
     .                               nd_srce, mytag, mycomm,
     .                               ireq_qb(nreq), ierr)
                     nreq = nreq + 1
                  end do
                  if (ivmx.ge.2) then
c                    receive vist3d data
                     mytag = mytag_qb(6) + nblf
                     call MPI_Irecv (qb(lig(nblf),1,2), nvalqb,
     .                               MY_MPI_REAL,
     .                               nd_srce, mytag, mycomm,
     .                              ireq_qb(nreq), ierr)
                     nreq = nreq + 1
                  end if
                  if (ivmx.ge.4) then
c                    receive turb. data            
                       do ll=1,nummem
                          mytag = mytag_qb(6+ll) + nblf
                          call MPI_Irecv (qb(lig(nblf),ll,3), nvalqb,
     .                                    MY_MPI_REAL,
     .                                    nd_srce, mytag, mycomm,
     .                                    ireq_qb(nreq), ierr)
                          nreq = nreq + 1
                       end do
                  end if
c
               else if (myid.eq.mblk2nd(nblf)) then
c
c                 send the data to non-local nodes
c
                  do inode = 1,nnodes
                     if (inode .ne. myid) then
                        nd_dest = inode
c                       send q data
                        do ll=1,5
                           mytag  = mytag_qb(ll) + nblf
                           call MPI_Send (qb(lig(nblf),ll,1), nvalqb,
     .                                    MY_MPI_REAL,
     .                                    nd_dest, mytag,
     .                                    mycomm, ierr)
                        end do
                        if (ivmx.ge.2) then
c                          send vist3d data
                           mytag = mytag_qb(6) + nblf
                           call MPI_Send (qb(lig(nblf),1,2), nvalqb,
     .                                    MY_MPI_REAL,
     .                                    nd_dest, mytag,
     .                                    mycomm, ierr)
                        end if
                        if (ivmx.ge.4) then
c                          send turb. data
                           do ll=1,nummem
                              mytag = mytag_qb(6+ll) + nblf
                              call MPI_Send (qb(lig(nblf),ll,3), nvalqb,
     .                                       MY_MPI_REAL,
     .                                       nd_dest, mytag,
     .                                       mycomm, ierr)
                           end do
                        end if
                     end if
                  end do
c
               end if
c
            end if
c
#endif
         else
c
c           add corrections to embedded grids
c           (nbl is global, nblz is embedded)
c
            if (mgflag.gt.1) then
               ipass = 1
               jc    = jdimg(nbl)
               kc    = kdimg(nbl)
               ic    = idimg(nbl)
               do 6030 nblz=1,nblock
               if (nblz.eq.nbl) go to 6030
               nblcc = nblcg(nblz)
               if (nblcc.eq.nbl) then
                  igridc=igridg(nblz)
                  if (iemg(igridc).gt.0) then
                     call lead(nblz,lw,lw2,maxbl)
                     iwk1  = 1
                     iwk2  = iwk1 + jdim*kdim*idim*5
                     iwk3  = iwk2 + jc*kc*ic*5
                     iwk4  = iwk3 + jdim*kc*ic
                     nroom = nwork-iwk4
                     mdim  = jdim*kdim*ic
                     if (nroom.lt.mdim) then
                        nou(1) = min(nou(1)+1,ibufdim)
                        write(bou(nou(1),1),3060)
                        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                     end if
#if defined DIST_MPI
                     nd_dest = mblk2nd(nblz)
                     nd_srce = mblk2nd(nbl)
c
                     if (nd_dest.ne.nd_srce) then
c
c                    embedded and coarser blocks on different procs
c
                     jki5 = jc*kc*ic*5
                     mytag   = itag_q1 + nbl
                     if (myid.eq.nd_srce) then
                        call MPI_Send(w(lq1c),jki5,
     .                       MY_MPI_REAL,nd_dest,
     .                       mytag,mycomm,ierr)
                     else if (myid .eq. nd_dest) then
                        iwk5  = iwk4 + jdim*kdim*ic
                        nroom = nwork-iwk5
                        if (nroom.lt.jki5) then
                           nou(1) = min(nou(1)+1,ibufdim)
                           write(bou(nou(1),1),3060)
                           call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                        end if
                        call MPI_Recv(wk(iwk5),jki5,
     .                  MY_MPI_REAL,nd_srce,
     .                  mytag,mycomm,istat,ierr)
                     end if
                     mytag   = itag_qc + nbl
                     if (myid.eq.nd_srce) then
                        call MPI_Send(w(lqc),jki5,
     .                       MY_MPI_REAL,nd_dest,
     .                       mytag,mycomm,ierr)
                     else if (myid .eq. nd_dest) then
                        iwk6  = iwk5 + jki5
                        nroom = nwork-iwk6
                        if (nroom.lt.jki5) then
                           nou(1) = min(nou(1)+1,ibufdim)
                           write(bou(nou(1),1),3060)
                           call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                        end if
                        call MPI_Recv(wk(iwk6),jki5,
     .                  MY_MPI_REAL,nd_srce,
     .                  mytag,mycomm,istat,ierr)
                     end if
                     if (myid .eq. mblk2nd(nblz)) then
                        if (isklton.gt.0) then
                           nou(1) = min(nou(1)+1,ibufdim)
                           write(bou(nou(1),1),3070)nbl,nblz,
     .                     igridg(nblz)
                           nou(1) = min(nou(1)+1,ibufdim)
                           write(bou(nou(1),1),3040)jdim,kdim,idim
                           nou(1) = min(nou(1)+1,ibufdim)
                           write(bou(nou(1),1),3050)jj2,kk2,ii2
                        end if
                        call add2x(w(lq),wk(iwk6),jdim,kdim,
     .                  idim,jc,kc,ic,wk(iwk5),wk(iwk1),
     .                  wk(iwk2),wk(iwk3),wk(iwk4),
     .                  js,ks,is,je,ke,ie,ipass,nbl,nblz,
     .                  nou,bou,nbuf,ibufdim,5,myid)
                     end if
c
                     else
c
c                    embedded and coarser block reside on same proc
c
                     if (myid .eq. mblk2nd(nblz)) then
#endif
                     if (isklton.gt.0) then
                        nou(1) = min(nou(1)+1,ibufdim)
                        write(bou(nou(1),1),3070)nbl,nblz,
     .                  igridg(nblz)
                        nou(1) = min(nou(1)+1,ibufdim)
                        write(bou(nou(1),1),3040)jdim,kdim,idim
                        nou(1) = min(nou(1)+1,ibufdim)
                        write(bou(nou(1),1),3050)jj2,kk2,ii2
                     end if
                     call add2x(w(lq),w(lqc),jdim,kdim,idim,
     .                          jc,kc,ic,w(lq1c),wk(iwk1),wk(iwk2),
     .                          wk(iwk3),wk(iwk4),js,ks,is,
     .                          je,ke,ie,ipass,nbl,nblz,
     .                          nou,bou,nbuf,ibufdim,5,myid)
#if defined DIST_MPI
                     end if
c
                     end if
#endif
                  end if
               end if
c
               if (isklton.gt.0) then
#   ifdef FASTIO
                  call writ_buffast(nblz,11,nou,bou,nbuf,ibufdim,myhost,
     .                          myid,mycomm,mblk2nd,maxbl,10)
#   else
                  call writ_buf(nblz,11,nou,bou,nbuf,ibufdim,myhost,
     .                          myid,mycomm,mblk2nd,maxbl)
#   endif
               end if
c
 6030          continue
            end if
         end if
c
         if (isklton.gt.0) then
#   ifdef FASTIO
            call writ_buffast(nbl,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .                    mycomm,mblk2nd,maxbl,11)
#   else
            call writ_buf(nbl,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .                    mycomm,mblk2nd,maxbl)
#   endif
         end if
c
 7020    continue
c
#if defined DIST_MPI
         if (nreq-1 .gt. 0) then
            do nnn=1,nreq-1
               call MPI_Wait(ireq_qb(nnn), istat, ierr)
            end do
         end if
c
#endif
         if (myid.eq.myhost) then
            if (level.lt.lglobal) then
               if (isklton.gt.0)write(11,3080)
            else if (mgflag.gt.1) then
               if (isklton.gt.0)write(11,3080)
            end if
         end if
c
         level  = level + 1
c
         if (level.eq.levt) then
c           finished with multigrid cycle
            if(isklton.eq.1 .and. myid.eq.myhost) write(11,3000)
            go to 6999
         end if
c
         nit    = mtt
         if (nsm(level-levb+1).gt.0) then
            nit    = mit(level-levb+1,iseq) + mtt
            kxpand = 1
         end if
         kode   = 1
         ntime  = 0
c
         if (nit.eq.0) then
c           continue with upward leg
            go to 7500
         else
c           done with upward leg; begin downward leg
            go to 9000
         end if
c
      else
c
c        embedded grids without multigrid; set level=levt
c        and start new cycle
c
         level = levt
c
      end if
c
c***********************************************
c     end upward leg of multigrid cycle
c***********************************************
c
 6999 continue
c
c     write restart file if non time-accurate
c
      if (real(dt).lt.0.0) then
         if (level.eq.levt) then
            if (icyc.ge.ncyc .or. icyc/nwrest*nwrest.eq.icyc) then
               if (myid .eq. myhost) then
c                 initialize eddy viscosity limit flag
                  if (icyc.ge.ncyc) then
                     iwarneddy  = 0
                  else
                     iwarneddy  = -1
                  end if
               if (icgns .ne. 1) then
                  open(unit=2,file=restrt,form='unformatted',
     .                 status='unknown')
                  rewind( 2)
                  if (ndwrt .ne. 0 .and. iunst .gt. 1) then
                    if(mseq .gt. 1) then
                      nou(1) = min(nou(1)+1,ibufdim)
                      write(bou(nou(1),1)
     .                   ,'(''stopping...cannot write dynamic '',
     .                       ''grid file when mseq > 1'')')
                      call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                    end if
                    open(98,file='dgplot3d.bin',form='formatted'
     .                  ,status='unknown')
                    write(98,'(i8)') ngrid
                    do igrid=1,ngrid
                      nbl = nblg(igrid)
                      iem = iemg(igrid)
                      if (iem.eq.0) nbl = nbl + (mseq-iseq)
                      call lead(nbl,lw,lw2,maxbl)
                      idimt(igrid) = idim
                      jdimt(igrid) = jdim
                      kdimt(igrid) = kdim
                    enddo
                    write(98,'(3i8)') (idimt(igrid),jdimt(igrid),
     .                           kdimt(igrid),igrid=1,ngrid)
                  end if
                  if( iclcd .eq. 1 .or. iclcd .eq. 2 ) then
                     open(unit=102,file=clcds,form='unformatted',
     .                    status='unknown')
                     rewind(102)
                  end if
               else
#if defined CGNS
c   Open database
               basedesired='Base'
               idimdesired=3
               call wopencgns(grid,basedesired,idimdesired,iccg,
     .              ibase,nzones)
               iflagg=0
               call winfot(iccg,ibase,dt)
#endif
               end if
               end if
               if (iteravg .eq. 1 .or. iteravg .eq. 2) then
                 if (myid.eq.myhost) then
                   write(11,2222)
                   write(11,2223)
                   write(11,2225) avgg
                   write(11,2225) avgq
                   rewind(96)
                   rewind(97)
                   write(96) ngrid
                   write(97) ngrid

                   if (ipertavg .eq. 1 .or. ipertavg .eq. 2)then
                      write(11,2225) avgq2 
                      write(11,2225) avgq2pert 
                      rewind(95) 
                      rewind(98) 
                      write(95) ngrid 
                      write(98) ngrid 
                      write(95) (idimg(nblg(igrid)+mseq-iseq), 
     .                     jdimg(nblg(igrid)+mseq-iseq), 
     .                     kdimg(nblg(igrid)+mseq-iseq),
     $                     igrid=1,ngrid) 
                      write(96) (idimg(nblg(igrid)+mseq-iseq), 
     .                     jdimg(nblg(igrid)+mseq-iseq), 
     .                     kdimg(nblg(igrid)+mseq-iseq),
     $                     igrid=1,ngrid) 
                      write(97) (idimg(nblg(igrid)+mseq-iseq), 
     .                     jdimg(nblg(igrid)+mseq-iseq), 
     .                     kdimg(nblg(igrid)+mseq-iseq),
     $                     igrid=1,ngrid) 
                      write(98) (idimg(nblg(igrid)+mseq-iseq), 
     .                     jdimg(nblg(igrid)+mseq-iseq), 
     .                     kdimg(nblg(igrid)+mseq-iseq),
     $                     igrid=1,ngrid) 
                   else
                      write(96) (idimg(nblg(igrid)+mseq-iseq)-1,
     .                     jdimg(nblg(igrid)+mseq-iseq)-1,
     .                     kdimg(nblg(igrid)+mseq-iseq)-1,
     $                     igrid=1,ngrid)
                      write(97) (idimg(nblg(igrid)+mseq-iseq)-1,
     .                     jdimg(nblg(igrid)+mseq-iseq)-1,
     .                     kdimg(nblg(igrid)+mseq-iseq)-1,
     $                     igrid=1,ngrid)
                   end if
c
                 end if
               end if
               do 6040 igrid=1,ngrid
               iskipz = 0
               if (igrid.eq.1) iskipz = 1
               nbl = nblg(igrid)
               iem = iemg(igrid)
               if (iseq.ne.mseq .and. iem.gt.0) go to 6040
               if (icgns .eq. 1) then
               call lead(nbl,lw,lw2,maxbl)
               idima=idim
               jdima=jdim
               kdima=kdim
               end if
               if (iem.eq.0) nbl = nbl + (mseq-iseq)
               if (mblk2nd(nbl) .eq. myid .or. myid .eq. myhost) then
                  call lead(nbl,lw,lw2,maxbl)
               if (icgns .eq. 1) then
                 if(nwork .lt. (idima+1)*(jdima+1)*(kdima+1)) then
                 write(901,'('' not enough memory for cgns Q write.'')')
                 write(901,'('' nwork in wk='',i6,''.  Needed = '',i6)')
     +              nwork,(idima+1)*(jdima+1)*(kdima+1)
                 call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                 end if
               else
                 idima=idim
                 jdima=jdim
                 kdima=kdim
               end if
                  call wrest(nbl,jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),
     .            w(lqi0),ncycmax,rms,clw,cdw,cdpw,cdvw,cxw,cyw,czw,
     .            cmxw,cmyw,cmzw,
     $            n_clcd,clcd,nblocks_clcd,blocks_clcd,
     $            fmdotw,cftmomw,cftpw,cftvw,cfttotw,
     .            rmstr,nneg,iskipz,w(lvis),w(lxib),
     .            w(lsnk0),w(lsni0),w(lxkb),w(lnbl),w(lcmuv),
     .            thetay,maxbl,myid,myhost,mycomm,mblk2nd,igrid,wk,
     .            idima,jdima,kdima,w(lbcj),w(lbck),w(lbci),w(lvj0),
     .            w(lvk0),w(lvi0),w(ltj0),w(ltk0),w(lti0),w(lblk),
     .            iwk,iwork,iovrlp(nbl),nou,bou,nbuf,ibufdim,
     .            w(lqavg),w(lq2avg),w(lx),w(ly),w(lz),nummem)
               end if
 6040          continue
               if (myid.eq.myhost) then
               if (icgns .ne. 1) then
                  close(2)
                  close(102)
               else
#if defined CGNS
c   Close data base
                 write(901,'('' Closing CGNS database from mgblk'')')
                 call cg_close_f(iccg, ier)
                 if (ier .ne. 0) then
                   call cg_error_print_f
                   call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                 end if

#endif
               end if
               end if
            end if
         end if
      end if

c
c     the following line is added to allow exit from
c     loop 7000 if the user has altered ntstep via
c     the "stop" file
c
6998  continue
c
      if (icyc.ge.ncyc) go to 7011
c

c
 7000 continue
c
7011  continue
c     recover original ita
      ita=itatemp
c
c***********************************************
c
c     end multigrid cycles
c
c***********************************************
c
c     advance time if time-accurate
c
      if (real(dt).gt.0.) then
         if (level.eq.levt) then
            time    = time  + dt
            time2mc = time2mc + dt
            do 800 nbl=1,nblock
            time2(nbl) = time2(nbl) + dt
  800       continue
            timekeep(ntt) = time
         end if
      end if
c
c     aeroelastic corrector step
c
      if (naesrf.gt.0) then
         call ae_corr(stm,stmi,xs,xxn,gforcn,gforcs,gforcnm,
     .                gf0,lw,lw2,w,mgwk,maxbl,maxseg,
     .                aesrfdat,nmds,maxaes,nt,mblk2nd,iseq,
     .                levelg,iadvance,nblock,icsi,icsf,jcsi,jcsf,
     .                kcsi,kcsf,myid,nsegdfrm,idfrmseg,iaesurf,
     .                perturb,aehist,ncycmax,maxsegdg,myhost,mycomm)
      end if
c
c     rigid body mode corrector step
c
      if (irbtrim.ne.0) irigb = 0
      if (irigb.eq.1) then
#if defined DIST_MPI
         call MPI_Bcast(clw(ntt),1,MY_MPI_REAL,myhost,
     .                  mycomm,ierr)
         call MPI_Bcast(cmyw(ntt),1,MY_MPI_REAL,myhost,
     .                  mycomm,ierr)
#endif
         call rb_corr(aesrfdat,clw(ntt),cmyw(ntt),maxaes)
      end if
c
c
c   meshdef = 1 if flow solution is to be bipassed, and mesh deformation only:
c
      if(meshdef.eq.1) goto 5990
c
c     write restart file if time-accurate
c
      if (real(dt).gt.0.0) then
         if (level.eq.levt) then
            if (nt.ge.ntstep .or. nt/nwrest*nwrest.eq.nt) then
               if (myid .eq. myhost) then
c                 initialize eddy viscosity limit flag
                  if (nt.ge.ntstep) then
                     iwarneddy  = 0
                  else
                     iwarneddy  = -1 
                  end if
               if (icgns .ne. 1) then
                  open(unit=2,file=restrt,form='unformatted',
     .                 status='unknown')
                  rewind( 2)
                  if( iclcd .eq. 1 .or. iclcd .eq. 2 ) then
                     open(unit=102,file=clcds,form='unformatted',
     .                    status='unknown')
                     rewind(102)
                  end if
               else
#if defined CGNS
c   Open database
               basedesired='Base'
               idimdesired=3
               call wopencgns(grid,basedesired,idimdesired,iccg,
     .              ibase,nzones)
#endif
               end if
               end if
               if (iteravg .eq. 1 .or. iteravg .eq. 2) then
                 if (myid.eq.myhost) then
                   write(11,2222)
                   write(11,2223)
                   write(11,2224)
                   write(11,2225) avgg
                   write(11,2225) avgq
                   rewind(96)
                   rewind(97)
                   write(96) ngrid
                   write(97) ngrid

                   if (ipertavg .eq. 1 .or. ipertavg .eq. 2) then
                      write(11,2225) avgq2 
                      write(11,2225) avgq2pert 
                      rewind(95) 
                      rewind(98) 
                      write(95) ngrid 
                      write(98) ngrid 
                      write(95) (idimg(nblg(igrid)+mseq-iseq), 
     .                     jdimg(nblg(igrid)+mseq-iseq), 
     .                     kdimg(nblg(igrid)+mseq-iseq),
     $                     igrid=1,ngrid) 
                      write(96) (idimg(nblg(igrid)+mseq-iseq), 
     .                     jdimg(nblg(igrid)+mseq-iseq), 
     .                     kdimg(nblg(igrid)+mseq-iseq),
     $                     igrid=1,ngrid) 
                      write(97) (idimg(nblg(igrid)+mseq-iseq), 
     .                     jdimg(nblg(igrid)+mseq-iseq), 
     .                     kdimg(nblg(igrid)+mseq-iseq),
     $                     igrid=1,ngrid) 
                      write(98) (idimg(nblg(igrid)+mseq-iseq), 
     .                     jdimg(nblg(igrid)+mseq-iseq), 
     .                     kdimg(nblg(igrid)+mseq-iseq),
     $                     igrid=1,ngrid) 
                   else
                      write(96) (idimg(nblg(igrid)+mseq-iseq)-1,
     .                     jdimg(nblg(igrid)+mseq-iseq)-1,
     .                     kdimg(nblg(igrid)+mseq-iseq)-1,
     $                     igrid=1,ngrid)
                      write(97) (idimg(nblg(igrid)+mseq-iseq)-1,
     .                     jdimg(nblg(igrid)+mseq-iseq)-1,
     .                     kdimg(nblg(igrid)+mseq-iseq)-1,
     $                     igrid=1,ngrid)
                   end if
                 end if
               end if
               do 6050 igrid=1,ngrid
               iskipz = 0
               if (igrid.eq.1) iskipz = 1
               nbl = nblg(igrid)
               iem = iemg(igrid)
               if (iseq.ne.mseq .and. iem.gt.0) go to 6050
               if (icgns .eq. 1) then
               call lead(nbl,lw,lw2,maxbl)
               idima=idim
               jdima=jdim
               kdima=kdim
               end if
               if (iem.eq.0) nbl = nbl + (mseq-iseq)
               if (mblk2nd(nbl) .eq. myid .or. myid .eq. myhost) then
                  call lead(nbl,lw,lw2,maxbl)
               if (icgns .eq. 1) then
                 if(nwork .lt. (idima+1)*(jdima+1)*(kdima+1)) then
                 write(901,'('' not enough memory for cgns Q write.'')')
                 write(901,'('' nwork in wk='',i6,''.  Needed = '',i6)')
     +              nwork,(idima+1)*(jdima+1)*(kdima+1)
                 call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                 end if
               else
                 idima=idim
                 jdima=jdim
                 kdima=kdim
               end if
                  call wrest(nbl,jdim,kdim,idim,w(lq),w(lqj0),w(lqk0),
     .            w(lqi0),ncycmax,rms,clw,cdw,cdpw,cdvw,cxw,cyw,czw,
     .            cmxw,cmyw,cmzw,
     $            n_clcd,clcd,nblocks_clcd,blocks_clcd,
     $            fmdotw,cftmomw,cftpw,cftvw,cfttotw,
     .            rmstr,nneg,iskipz,w(lvis),w(lxib),
     .            w(lsnk0),w(lsni0),w(lxkb),w(lnbl),w(lcmuv),
     .            thetay,maxbl,myid,myhost,mycomm,mblk2nd,igrid,wk,
     .            idima,jdima,kdima,w(lbcj),w(lbck),w(lbci),w(lvj0),
     .            w(lvk0),w(lvi0),w(ltj0),w(ltk0),w(lti0),w(lblk),
     .            iwk,iwork,iovrlp(nbl),nou,bou,nbuf,ibufdim,
     .            w(lqavg),w(lq2avg),w(lx),w(ly),w(lz),nummem)
               end if
 6050          continue
c
c              set flag for writing additional data 2nd order (in time)
c              and/or dynamic mesh cases
c
               iflagg=0
               if (abs(ita) .eq. 2 .and. iunst .eq. 0) iflagg=1
               if (abs(ita) .ne. 2 .and. iunst .gt. 0) iflagg=2
               if (abs(ita) .eq. 2 .and. iunst .gt. 0) iflagg=3
               if (myid.eq.myhost) then
               if (icgns .ne. 1) then
                  write(2) iflagg
               else
#if defined CGNS
                  call winfot(iccg,ibase,dt)
#endif
               end if
               end if
c
c              write 2nd order data
c
               if (iflagg .eq. 1 .or. iflagg .eq. 3) then
               do 6060 igrid=1,ngrid
               nbl = nblg(igrid)
               iem = iemg(igrid)
               if (iseq.ne.mseq .and. iem.gt.0) go to 6060
               if (icgns .eq. 1) then
               call lead(nbl,lw,lw2,maxbl)
               idima=idim
               jdima=jdim
               kdima=kdim
               end if
               if (iem.eq.0) nbl = nbl + (mseq-iseq)
               if (mblk2nd(nbl) .eq. myid .or. myid .eq. myhost) then
                  call lead(nbl,lw,lw2,maxbl)
               if (icgns .eq. 1) then
                 if(nwork .lt. (idima+1)*(jdima+1)*(kdima+1)) then
                 write(901,'('' not enough memory for cgns'',
     +              '' unsteady Q write.'')')
                 write(901,'('' nwork in wk='',i6,''.  Needed = '',i6)')
     +              nwork,(idima+1)*(jdima+1)*(kdima+1)
                 call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                 end if
               else
                 idima=idim
                 jdima=jdim
                 kdima=kdim
               end if
                  call wrestg(nbl,jdim,kdim,idim,w(lx),w(ly),w(lz),
     .                 w(lxnm2),w(lynm2),w(lznm2),w(ldeltj),w(ldeltk),
     .                 w(ldelti),w(lqc0),0,0,utrans,vtrans,wtrans,
     .                 omegax,omegay,omegaz,xorig,yorig,zorig,
     .                 dxmx,dymx,dzmx,dthxmx,dthymx,
     .                 dthzmx,thetax,thetay,thetaz,rfreqt,
     .                 rfreqr,xorig0,yorig0,zorig0,time2,
     .                 thetaxl,thetayl,thetazl,itrans,irotat,idefrm,
     .                 utrnsae,vtrnsae,wtrnsae,omgxae,omgyae,omgzae,
     .                 xorgae,yorgae,zorgae,thtxae,thtyae,thtzae,
     .                 rfrqtae,rfrqrae,icsi,icsf,jcsi,jcsf,
     .                 kcsi,kcsf,freq,gmass,damp,x0,gf0,nmds,maxaes,
     .                 aesrfdat,perturb,myhost,myid,mycomm,mblk2nd,
     .                 maxbl,nsegdfrm,idfrmseg,iaesurf,maxsegdg,wk,
     .                 nwork,idima,jdima,kdima,igrid,w(lxib2),nummem)
               end if
 6060          continue
               end if
c
c              write dynamic mesh data
c
               if (iflagg .eq. 2 .or. iflagg .eq. 3) then
                  do 6070 igrid=1,ngrid
                  nbl = nblg(igrid)
                  iem = iemg(igrid)
                  if (iseq.ne.mseq .and. iem.gt.0) go to 6070
               if (icgns .eq. 1) then
               call lead(nbl,lw,lw2,maxbl)
               idima=idim
               jdima=jdim
               kdima=kdim
               end if
                  if (iem.eq.0) nbl = nbl + (mseq-iseq)
                  if (mblk2nd(nbl) .eq. myid .or. myid .eq. myhost) then
                     call lead(nbl,lw,lw2,maxbl)
               if (icgns .eq. 1) then
                 if(nwork .lt. (idima+1)*(jdima+1)*(kdima+1)) then
                 write(901,'('' not enough memory for cgns'',
     +              '' unsteady Q write.'')')
                 write(901,'('' nwork in wk='',i6,''.  Needed = '',i6)')
     +              nwork,(idima+1)*(jdima+1)*(kdima+1)
                 call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                 end if
               else
                 idima=idim
                 jdima=jdim
                 kdima=kdim
               end if
                     iuns = max(itrans(nbl),irotat(nbl),idefrm(nbl))
                     call wrestg(nbl,jdim,kdim,idim,w(lx),w(ly),w(lz),
     .                    w(lxnm2),w(lynm2),w(lznm2),w(ldeltj),
     .                    w(ldeltk),w(ldelti),w(lqc0),1,iuns,utrans,
     .                    vtrans,wtrans,omegax,omegay,omegaz,xorig,
     .                    yorig,zorig,dxmx,dymx,dzmx,dthxmx,dthymx,
     .                    dthzmx,thetax,thetay,thetaz,rfreqt,
     .                    rfreqr,xorig0,yorig0,zorig0,time2,
     .                    thetaxl,thetayl,thetazl,itrans,irotat,idefrm,
     .                    utrnsae,vtrnsae,wtrnsae,omgxae,omgyae,omgzae,
     .                    xorgae,yorgae,zorgae,thtxae,thtyae,thtzae,
     .                    rfrqtae,rfrqrae,icsi,icsf,jcsi,jcsf,
     .                    kcsi,kcsf,freq,gmass,damp,x0,gf0,nmds,maxaes,
     .                    aesrfdat,perturb,myhost,myid,mycomm,mblk2nd,
     .                    maxbl,nsegdfrm,idfrmseg,iaesurf,maxsegdg,wk,
     .                    nwork,idima,jdima,kdima,igrid,w(lxib2),nummem)
                  end if
 6070             continue
#if defined CGNS
               else
                  if (icgns .eq. 1 .and. myid .eq. myhost) then
                  do igrid=1,ngrid
c       need to delete certain nodes if they exist
                    call cg_goto_f(iccg,ibase,ier,'Zone_t',igrid,'end')
                    call cg_delete_node_f('RigidGridMotion',ier)
                    call cg_delete_node_f('ArbitraryGridMotion',ier)
                    call cg_delete_node_f('MovedGrid',ier)
                    call cg_delete_node_f('MovedGridLastDT',ier)
                  enddo
                  end if
#endif
               end if
c
c              write modal data if aeroelastic
c
               if (myid .eq. myhost) then
                   if (naesrf.gt.0) then
                      if (icgns .ne. 1) then
                      write(2) (timekeep(nn),nn=1,ntt)
                      do iaes=1,naesrf
                         nmodes = aesrfdat(5,iaes)
                         write(2) (xxn(n,iaes),n=1,2*nmodes)
                         write(2) (gforcn(n,iaes),n=1,2*nmodes)
                         write(2) (gforcnm(n,iaes),n=1,2*nmodes)
                         write(2) (((aehist(nn,ll,n,iaes),nn=1,ntt),
     .                               ll=1,3),n=1,nmodes)
                      end do
                      else
#if defined CGNS
                         maxnum=2*nmds*naesrf
                         maxnum2=ntt*3*nmds*naesrf
                         if ((maxnum+maxnum2) .gt. nwork) then
                           write(901,'('' not enough room in wk for'',
     .                       '' writing data in waeromode'')')
                           call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                         end if
                         call waeromode(iccg,ibase,timekeep,xxn,gforcn,
     .                    gforcnm,aehist,ncycmax,nmds,maxaes,
     .                    ntt,naesrf,wk,maxnum,wk(maxnum+1),maxnum2)
#endif
                      end if
#if defined CGNS
                   else
                  if (icgns .eq. 1) then
c       need to delete CFL3DAeroModeData node if it exists
                     call cg_goto_f(iccg,ibase,ier,'end')
                     call cg_delete_node_f('CFL3DAeroModeData',ier)
                  end if
#endif
                   end if
               if (icgns .ne. 1) then
                   close(2)
                   if (ndwrt .ne. 0) then
                      close(98)
                   end if
               else
#if defined CGNS
c   Close data base
                 write(901,'('' Closing CGNS database from mgblk'')')
                 call cg_close_f(iccg, ier)
                 if (ier .ne. 0) then
                   call cg_error_print_f
                   call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                 end if
#endif
               end if
               end if
c
            end if
         end if
      end if
c
c     output new plot3d file every movie time steps
c
      movabs = abs(movie)
      if (ipertavg .eq. 0) then
         if (real(dt).gt.0. .and. movabs.gt.0) then
            if ((movabs.gt.0 .and. (nt.eq.nt/movabs*movabs)) .or.
     .           (iteravg.gt.0)) then 
               if (nt.eq.nt/movabs*movabs) then
                  lhdr   = 0
                  iwk1   = 1
                  iwk2   = iwk1  + 3*nplots
                  iwk3   = iwk2  + maxbl
                  iwork1 = iwork - iwk3 + 1
                  if (iwork1.le.0) then
                     nou(1) = min(nou(1)+1,ibufdim)
                     write(bou(nou(1),1),'(''stopping...not '',
     $                    ''enough integer '',
     .                    ''work space for subroutine qout'')')
                     call termn8(myid,-1,ibufdim,nbuf,bou,nou)
                  end if
                  do iii = 1,iwk3
                     iwk(iii) = 0
                  end do
                  do iii = 1,nwork
                     wk(iii) = 0.0
                  end do
                  call qout(iseq,lw,lw2,w,mgwk,wk,nwork,
     .                 nplots,iovrlp,iibg,kkbg,jjbg,ibcg,lbg,
     .                 ibpntsg,qb,lwdat,nbci0,nbcj0,nbck0,
     .                 nbcidim,nbcjdim,nbckdim,jbcinfo,kbcinfo,
     .                 ibcinfo,bcfilei,bcfilej,bcfilek,itrans,
     .                 irotat,idefrm,nblock,levelg,igridg,iviscg,jdimg,
     .                 kdimg,idimg,nblg,clw,ncycmax,nplot3d,
     .                 inpl3d,ip3dsurf,nprint,inpr,iadvance,mycomm,
     .                 myid,myhost,mblk2nd,nou,bou,nbuf,ibufdim,maxbl,
     .                 maxgr,maxseg,iitot,jsg,ksg,isg,jeg,keg,ieg,
     .                 ninter,windex,iindex,nblkpt,intmax,nsub1,maxxe,
     .                 nblk,nbli,limblk,isva,nblon,mxbli,thetay,
     .                 iwk(iwk1),iwk(iwk2),iwk(iwk3),iwork1,xorig,yorig,
     .                 zorig,period_miss,geom_miss,epsc0,epsrot,
     .                 isav_blk,isav_pat,isav_pat_b,isav_emb,isav_prd,
     .                 lbcprd,lbcemb,dthetxx,dthetyy,dthetzz,nblcg,lfgm,
     .                 istat2_bl,istat2_pa,istat2_pe,istat2_em,
     .                 istat_size,vormax,ivmax,jvmax,kvmax,bcfiles,
     .                 mxbcfil,iprnsurf,nummem)

                  do iii = 1,iwk3
                     iwk(iii) = 0
                  end do
                  do iii = 1,nwork
                     wk(iii) = 0.0
                  end do
                  nframes = nframes + 1
               end if
            end if
         end if
      end if
      if (ipertavg .eq. 1 .or. ipertavg .eq. 2) then
         if ((real(dt).gt.0. .and. movabs.gt.0).or.
     $        (iteravg.gt.0)) then
c     
            if ((movabs.gt.0 .and. (nt.eq.nt/movabs*movabs)) .or.
     .           (iteravg.gt.0)) then
               lhdr   = 0
               iwk1   = 1
               iwk2   = iwk1  + 3*nplots
               iwk3   = iwk2  + maxbl
               iwork1 = iwork - iwk3 + 1
               if (iwork1.le.0) then
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),
     $                 '(''stopping...not enough integer '',
     .                    ''work space for subroutine qoutavg'')')
                  call termn8(myid,-1,ibufdim,nbuf,bou,nou)
               end if
               do iii = 1,iwk3
                  iwk(iii) = 0
               end do
               do iii = 1,nwork
                  wk(iii) = 0.0
               end do
c     
               call qoutavg(iseq,lw,lw2,w,mgwk,wk,nwork, 
     .              nplots,iovrlp,iibg,kkbg,jjbg,ibcg,lbg,
     .              ibpntsg,qb,lwdat,nbci0,nbcj0,nbck0,
     .              nbcidim,nbcjdim,nbckdim,jbcinfo,kbcinfo,
     .              ibcinfo,bcfilei,bcfilej,bcfilek,itrans,
     .              irotat,idefrm,nblock,levelg,igridg,iviscg,jdimg,
     .              kdimg,idimg,nblg,clw,ncycmax,nplot3d,
     .              inpl3d,ip3dsurf,nprint,inpr,iadvance,mycomm,
     .              myid,myhost,mblk2nd,nou,bou,nbuf,ibufdim,maxbl,
     .              maxgr,maxseg,iitot,jsg,ksg,isg,jeg,keg,ieg,
     .              ninter,windex,iindex,nblkpt,intmax,nsub1,maxxe,
     .              nblk,nbli,limblk,isva,nblon,mxbli,thetay,
     .              iwk(iwk1),iwk(iwk2),iwk(iwk3),iwork1,xorig,yorig,
     .              zorig,period_miss,geom_miss,epsc0,epsrot,
     .              isav_blk,isav_pat,isav_pat_b,isav_emb,isav_prd,
     .              lbcprd,lbcemb,dthetxx,dthetyy,dthetzz,nblcg,lfgm,
     .              istat2_bl,istat2_pa,istat2_pe,istat2_em,
     .              istat_size,vormax,ivmax,jvmax,kvmax,bcfiles,
     .              mxbcfil,iprnsurf,nt,movabs,nummem)
               do iii = 1,iwk3
                  iwk(iii) = 0
               end do
               do iii = 1,nwork
                  wk(iii) = 0.0
               end do
               nframes = nframes + 1
            end if
         end if
      end if
c     
      if (icoarsemovie .ne. 0) then
c     
         iflag_coarse = icoarsemovie
         mov_coarse = icoarsemovie
         iinc_coarse = inc_coarse(1)
         jinc_coarse = inc_coarse(2)
         kinc_coarse = inc_coarse(3)
         
         if (real(dt).gt.0. .and. iflag_coarse.gt.0) then                  
            if (ntt.eq.ntt/mov_coarse*mov_coarse) then                      
               lhdr   = 0
               iwk1   = 1
               iwk2   = iwk1  + 3*nblocks/lglobal
               iwk3   = iwk2  + maxbl
               iwork1 = iwork - iwk3 + 1
               if (iwork1.le.0) then
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),
     $                 '(''stopping...not enough integer '',
     .                 ''work space for subroutine qout_coarse'')')
                  call termn8(myid,-1,ibufdim,nbuf,bou,nou)
               end if
               do iii = 1,iwk3
                  iwk(iii) = 0
               end do
               do iii = 1,nwork
                  wk(iii) = 0.0
               end do
               
               if (nt.eq.1) then
                  icallcrs=0
                  icall2d=0
               endif
               call qout_coarse(iseq,lw,lw2,w,mgwk,wk,nwork,
     .              nplots,iovrlp,iibg,kkbg,jjbg,ibcg,lbg,
     .              ibpntsg,qb,lwdat,nbci0,nbcj0,nbck0,
     .              nbcidim,nbcjdim,nbckdim,jbcinfo,kbcinfo,
     .              ibcinfo,bcfilei,bcfilej,bcfilek,itrans,
     .              irotat,idefrm,nblock,levelg,igridg,iviscg,jdimg,
     .              kdimg,idimg,nblg,clw,ncycmax,nplot3d,
     .              inpl3d,ip3dsurf,nprint,inpr,iadvance,mycomm,
     .              myid,myhost,mblk2nd,nou,bou,nbuf,ibufdim,maxbl,
     .              maxgr,maxseg,iitot,jsg,ksg,isg,jeg,keg,ieg,
     .              ninter,windex,iindex,nblkpt,intmax,nsub1,maxxe,
     .              nblk,nbli,limblk,isva,nblon,mxbli,thetay,
     .              iwk(iwk1),iwk(iwk2),iwk(iwk3),iwork1,xorig,yorig,
     .              zorig,period_miss,geom_miss,epsc0,epsrot,
     .              isav_blk,isav_pat,isav_pat_b,isav_emb,isav_prd,
     .              lbcprd,lbcemb,dthetxx,dthetyy,dthetzz,nblcg,lfgm,
     .              istat2_bl,istat2_pa,istat2_pe,istat2_em,
     .              istat_size,vormax,ivmax,jvmax,kvmax,bcfiles,
     .              mxbcfil,iprnsurf,nt,
     .              mov_coarse,iinc_coarse,jinc_coarse,kinc_coarse,
     .              nummem)
               do iii = 1,iwk3
                  iwk(iii) = 0
               end do
               do iii = 1,nwork
                  wk(iii) = 0.0
               end do
            end if
         end if
c     
      end if
c     
      if (i2dmovie .ne. 0) then
c
         iflag_2d = i2dmovie
c     ! output snapshot every mov_2d time steps  
         mov_2d = i2dmovie
         
         iinc_2d = inc_2d(1)
         jinc_2d = inc_2d(2)
         kinc_2d = inc_2d(3)
         
         if (real(dt).gt.0. .and. iflag_2d.gt.0) then                  
            if (ntt.eq.ntt/mov_2d*mov_2d) then                        
               lhdr   = 0
               iwk1   = 1
               iwk2   = iwk1  + 3*nblocks/lglobal                         
               iwk3   = iwk2  + maxbl
               iwork1 = iwork - iwk3 + 1
               if (iwork1.le.0) then
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),
     $                 '(''stopping...not enough integer '',
     .                 ''work space for subroutine qout_2d'')')
                  call termn8(myid,-1,ibufdim,nbuf,bou,nou)
               end if
               do iii = 1,iwk3
                  iwk(iii) = 0
               end do
               do iii = 1,nwork
                  wk(iii) = 0.0
               end do
               call qout_2d(iseq,lw,lw2,w,mgwk,wk,nwork,              
     .              nplots,iovrlp,iibg,kkbg,jjbg,ibcg,lbg,
     .              ibpntsg,qb,lwdat,nbci0,nbcj0,nbck0,
     .              nbcidim,nbcjdim,nbckdim,jbcinfo,kbcinfo,
     .              ibcinfo,bcfilei,bcfilej,bcfilek,itrans,
     .              irotat,idefrm,nblock,levelg,igridg,iviscg,jdimg,
     .              kdimg,idimg,nblg,clw,ncycmax,nplot3d,
     .              inpl3d,ip3dsurf,nprint,inpr,iadvance,mycomm,
     .              myid,myhost,mblk2nd,nou,bou,nbuf,ibufdim,maxbl,
     .              maxgr,maxseg,iitot,jsg,ksg,isg,jeg,keg,ieg,
     .              ninter,windex,iindex,nblkpt,intmax,nsub1,maxxe,
     .              nblk,nbli,limblk,isva,nblon,mxbli,thetay,
     .              iwk(iwk1),iwk(iwk2),iwk(iwk3),iwork1,xorig,yorig,
     .              zorig,period_miss,geom_miss,epsc0,epsrot,
     .              isav_blk,isav_pat,isav_pat_b,isav_emb,isav_prd,
     .              lbcprd,lbcemb,dthetxx,dthetyy,dthetzz,nblcg,lfgm,
     .              istat2_bl,istat2_pa,istat2_pe,istat2_em,
     .              istat_size,vormax,ivmax,jvmax,kvmax,bcfiles,
     .              mxbcfil,iprnsurf,nt,
     .              mov_2d,iinc_2d,jinc_2d,kinc_2d,nummem)
               
               do iii = 1,iwk3
                  iwk(iii) = 0
               end do
               do iii = 1,nwork
                  wk(iii) = 0.0
               end do
            end if
         end if
      end if
c
5990  continue
c
c     the following line is added to allow exit from
c     loop 6000 if the user has altered ntstep via
c     the "stop" file
c
      if (nt.ge.ntstep) go to 6011
c

 6000 continue
c
 6011 continue
c
      if (myid.eq.myhost) then
         write(11,3090) iseq
      end if
c
c********************************************************************
c
c     end time advancement
c
c********************************************************************
c
c
 1000 format(2x,49hstopping - too many levels in multigrid - max = 5)
 1010 format(/40h***** BEGINNING TIME ADVANCEMENT, iseq =,
     .       i2,6h *****/)
 1020 format(1x,25hsteady-state computations)
 1030 format(1x,47htime-accurate computations - 1st order in time,,
     .       1x,22hnumber of time steps: ,i4)
 1040 format(1x,47htime-accurate computations - 2nd order in time,,
     .       1x,22hnumber of time steps: ,i4)
 1050 format(1x,46hpseudo-time term included on LHS for stability,
     .       1x,13hwith large dt)
 1055 format(/,1x,43hlow Mach number preconditioning used, with:,/,
     .       1x,21h  cprec, uref, avn = ,1x,f7.4,1h,,1x,f7.4,1h,,
     .       1x,f7.4)
 1060 format(/37h***** BEGINNING MULTIGRID CYCLE *****/)
 1070 format(6h iseq=,i5)
 1080 format(12h level top =,i3)
 1090 format(15h level bottom =,i3)
 2000 format(33h number of embedded grid levels =,i3)
 2010 format(31h number of global grid levels =,i3)
 2020 format(9h lglobal=,i3)
 2030 format(39h stopping...increase dimension of nptsr)
 2040 format(/49h***** BEGINNING RESIDUAL/RESTRICTION LOOP, level,,
     .       13h kode, mode =,3i3)
 2050 format(1x,48h insufficient memory in mgblk *stopping*  lsafe=,i7)
 2060 format(1x,43hinterpolating for cells overlapped to block,i4)
 2070 format(43h not enough work space for subroutine collq)
 2080 format(44h not enough work space for subroutine coll2q)
 2090 format(/38h***** ENDING RESIDUAL/RESTRICTION LOOP,3i3/)
 2222 format(/,50h writing unformatted running-average cell-centered,
     .       30h PLOT3D files (MG,3D,iblanked))
 2223 format(38h   in same precision that code was run)
 2224 format(54h   grid is LATEST cell-center grid (not time-averaged))
 2225 format(8x,a80)
 3000 format(/34h***** ENDING MULTIGRID CYCLE *****)
 3010 format(/47h***** BEGINNING PROLONGATION LOOP, level, kode,,
     .       6h mode=,3i3/)
 3020 format(42h not enough work space for subroutine addx)
 3030 format(44h interpolating correction from coarser block,i4,
     .       17h to   finer block,i4,6h (grid,i4,1h))
 3040 format(31h   jdim,kdim,idim (finer grid)=,3i5)
 3050 format(31h   jj2,kk2,ii2  (coarser grid)=,3i5)
 3060 format(26h not enough work space for,17h subroutine add2x)
 3070 format(44h interpolating correction from coarser block,i4,
     .       17h to embeded block,i4,6h (grid,i4,1h))
 3080 format(/30h***** ENDING PROLONGATION LOOP)
 3090 format(/37h***** ENDING TIME ADVANCEMENT, iseq =,i2,6h *****/)
c
      return
      end
