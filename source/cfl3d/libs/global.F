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
      subroutine global(myid,maxbl,maxgr,maxseg,maxcs,nplots,mxbli,
     .                  bcvali,bcvalj,bcvalk,nbci0,nbcj0,nbck0,
     .                  nbcidim,nbcjdim,nbckdim,ibcinfo,jbcinfo,
     .                  kbcinfo,bcfilei,bcfilej,bcfilek,nblk,nbli,
     .                  limblk,isva,nblon,rkap0g,nblock,levelg,
     .                  igridg,iflimg,ifdsg,iviscg,jdimg,kdimg,
     .                  idimg,idiagg,nblcg,idegg,jsg,ksg,isg,jeg,
     .                  keg,ieg,mit,ilamlog,ilamhig,jlamlog,
     .                  jlamhig,klamlog,klamhig,iwfg,utrans,vtrans,
     .                  wtrans,omegax,omegay,omegaz,xorig,yorig,
     .                  zorig,dxmx,dymx,dzmx,dthxmx,dthymx,
     .                  dthzmx,thetax,thetay,thetaz,rfreqt,rfreqr,
     .                  xorig0,yorig0,zorig0,itrans,irotat,idefrm,
     .                  ngrid,ncgg,nblg,iemg,inewgg,iovrlp,ninter,
     .                  nplot3d,inpl3d,ip3dsurf,nprint,inpr,
     .                  iadvance,iforce,lfgm,ncs,icsinfo,ihstry,
     .                  ncycmax,iv,time2,thetaxl,thetayl,thetazl,
     .                  intmax,nsub1,iindex,lig,lbg,ibpntsg,
     .                  iipntsg,icall,iunit11,nou,bou,ibufdim,nbuf,
     .                  mglevg,nemgl,ipl3dtmp,ntr,bcfiles,mxbcfil,
     .                  utrnsae,vtrnsae,wtrnsae,omgxae,omgyae,
     .                  omgzae,xorgae,yorgae,zorgae,thtxae,thtyae,
     .                  thtzae,rfrqtae,rfrqrae,icsi,icsf,jcsi,jcsf,
     .                  kcsi,kcsf,freq,gmass,damp,x0,gf0,nmds,maxaes,
     .                  aesrfdat,perturb,iskip,jskip,kskip,nsegdfrm,
     .                  idfrmseg,iaesurf,maxsegdg,xorgae0,yorgae0,
     .                  zorgae0,icouple,iprnsurf)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Read in the case input data.
c***********************************************************************
c
c     maxgr    - maximum number of grids
c     maxbl    - maximum number of blocks
c     nplots   - maximum number of data sets to output via PLOT3D or 
c                print options
c     mxbli    - maximum number of 1-1 interfaces (incl. coarser blocks!)
c     maxseg   - maximum number of segments on a boundary
c     maxsegdg - maximum number of segments on a deforming 
c     ncycmax  - maximum number of time-steps/cycles
c     intmax   - maximum number of block interpolations
c     nsub1    - maximum number of blocks a single patch face 
c                may be interpolated from
c     mxbcfil  - maximum number (less one) of auxiliary files that may 
c                be specified for 2000 series bc's 
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*80  grid,plt3dg,plt3dq,output,residual,turbres,blomx,
     .              output2,printout,pplunge,ovrlap,patch,restrt,
     .              subres,subtur,grdmov,alphahist,errfile
      character*80 bcfiles(mxbcfil)
      character*10 datahdr(10)
      character*120 bou(ibufdim,nbuf)
c
      integer bcfilei,bcfilej,bcfilek
c
#   ifdef CMPLX
      complex lreftra,lrefrot,lrefdef
#   else
      real lreftra,lrefrot,lrefdef
#   endif
c
      real realval(20)
c
      dimension nou(nbuf)
      dimension mglevg(maxbl),nemgl(maxbl)
      dimension im(11),iv(maxbl),iindex(intmax,2*nsub1+9)
      dimension lig(maxbl),lbg(maxbl),ibpntsg(maxbl,4),iipntsg(maxbl)
      dimension inpl3d(nplots,11),inpr(nplots,11),ipl3dtmp(11,nplots)
      dimension bcvali(maxbl,maxseg,12,2),bcvalj(maxbl,maxseg,12,2),
     .          bcvalk(maxbl,maxseg,12,2),nbci0(maxbl),nbcidim(maxbl),
     .          nbcj0(maxbl),nbcjdim(maxbl),nbck0(maxbl),nbckdim(maxbl),
     .          jbcinfo(maxbl,maxseg,7,2),ibcinfo(maxbl,maxseg,7,2),
     .          kbcinfo(maxbl,maxseg,7,2),bcfilei(maxbl,maxseg,2),
     .          bcfilej(maxbl,maxseg,2),bcfilek(maxbl,maxseg,2)
      dimension nblk(2,mxbli),limblk(2,6,mxbli),
     .          isva(2,2,mxbli),nblon(mxbli),iovrlp(maxbl)
      dimension rkap0g(maxbl,3),levelg(maxbl),igridg(maxbl),
     .          iflimg(maxbl,3),ifdsg(maxbl,3),iviscg(maxbl,3),
     .          jdimg(maxbl),kdimg(maxbl),idimg(maxbl),idiagg(maxbl,3),
     .          nblcg(maxbl),idegg(maxbl,3),
     .          jsg(maxbl),ksg(maxbl),isg(maxbl),jeg(maxbl),keg(maxbl),
     .          ieg(maxbl),mit(5,maxbl),ilamlog(maxbl),ilamhig(maxbl),
     .          jlamlog(maxbl),jlamhig(maxbl),klamlog(maxbl),
     .          klamhig(maxbl),iwfg(maxbl,3)
      dimension utrans(maxbl),vtrans(maxbl),wtrans(maxbl),omegax(maxbl),
     .          omegay(maxbl),omegaz(maxbl),xorig(maxbl),yorig(maxbl),
     .          zorig(maxbl),dxmx(maxbl),dymx(maxbl),dzmx(maxbl),
     .          dthxmx(maxbl),dthymx(maxbl),dthzmx(maxbl),
     .          thetax(maxbl),thetay(maxbl),thetaz(maxbl),rfreqt(maxbl),
     .          rfreqr(maxbl),xorig0(maxbl),yorig0(maxbl),zorig0(maxbl),
     .          time2(maxbl),thetaxl(maxbl),thetayl(maxbl),
     .          thetazl(maxbl),itrans(maxbl),irotat(maxbl),idefrm(maxbl)
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
      dimension freq(nmds,maxaes),gmass(nmds,maxaes),x0(2*nmds,maxaes),
     .          gf0(2*nmds,maxaes),damp(nmds,maxaes),
     .          perturb(nmds,maxaes,4)
      dimension aesrfdat(5,maxaes),iskip(maxbl,500),jskip(maxbl,500),
     .          kskip(maxbl,500)
      dimension nskpi1(maxbl),nskpj1(maxbl),nskpk1(maxbl)
      dimension nsegdfrm(maxbl),idfrmseg(maxbl,maxsegdg),
     .          iaesurf(maxbl,maxsegdg)
      dimension ncgg(maxgr),nblg(maxgr),iemg(maxgr),
     .          inewgg(maxgr)
      dimension iadvance(maxbl),iforce(maxbl),icsinfo(maxcs,9)
c
      integer stats
      real, dimension(:,:), allocatable :: riskp,rjskp,rkskp
      integer, dimension(:,:), allocatable :: iskipt,jskipt,kskipt
c
      common /avgdata/ xnumavg,iteravg,xnumavg2,ipertavg,iclcd,isubit_r
      common /cfl/ dt0,dtold
      common /complx/ xmach_img,alpha_img,beta_img,reue_img,tinf_img,
     .                geom_img,surf_img,xrotrate_img,yrotrate_img,
     .                zrotrate_img
      common /des/ cdes,ides,cddes
      common /elastic/ ndefrm,naesrf
      common /elastic_ss/ idef_ss
      common /params/ lmaxgr,lmaxbl,lmxseg,lmaxcs,lnplts,lmxbli,lmaxxe,
     .                lnsub1,lintmx,lmxxe,liitot,isum,lncycm,
     .                isum_n,lminnode,isumi,isumi_n,lmptch,
     .                lmsub1,lintmax,libufdim,lnbuf,llbcprd,
     .                llbcemb,llbcrad,lnmds,lmaxaes,lnslave,lmxsegdg,
     .                lnmaster
      common /precond/ cprec,uref,avn
      common /cpurate/ rate(5),ratesub(5),ncell(20)
      common /chk/ ichk
      common /degshf/ ideg(3)
      common /fsum/ sref,cref,bref,xmc,ymc,zmc
      common /fvfds/ rkap0(3),ifds(3)
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /lam/ ilamlo,ilamhi,jlamlo,jlamhi,klamlo,klamhi,
     .        i_lam_forcezero
      common /maxiv/ ivmx
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /mgv/ epsssc(3),epsssr(3),issc,issr
      common /ncyct/ ncyctot
      common /alphait/ ialphit,cltarg,rlxalph,dalim,dalpha,icycupdt
      common /reyue/ reue,tinf,ivisc(3)
      common /sminn/ isminc,ismincforce
      common /twod/ i2d
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /wallfun/ iwf(3)
      common /wrbl/ nwrest
      common /wrestq/ irest,irest2
      common /moov/movie,nframes,icall1,lhdr
      common /igrdtyp/ ip3dgrd,ialph
      common /conversion/ radtodeg
      common /unit5/ iunit5
      common /motionmc/ xmc0,ymc0,zmc0,utransmc,vtransmc,wtransmc,
     .                  omegaxmc,omegaymc,omegazmc,xorigmc,yorigmc,
     .                  zorigmc,xorig0mc,yorig0mc,zorig0mc,thetaxmc,
     .                  thetaymc,thetazmc,dxmxmc,dymxmc,dzmxmc,
     .                  dthxmxmc,dthymxmc,dthzmxmc,rfreqtmc,
     .                  rfreqrmc,itransmc,irotatmc,time2mc
      common /deformz/ beta1,beta2,alpha1,alpha2,isktyp,negvol,meshdef,
     .                 nsprgit,ndgrd,ndwrt 
      common /turbconv/ cflturb(7),edvislim,iturbprod,nsubturb,nfreeze,
     .                  iwarneddy,itime2read,itaturb,tur1cut,tur2cut,
     .                  iturbord,tur1cutlev,tur2cutlev
      common /konew/ ikoprod,isstdenom,pklimterm,ibeta8kzeta,i_bsl,
     .        keepambient,re_thetat0,i_wilcox06,i_wilcox06_chiw,
     .        i_turbprod_kterm,i_catris_kw,prod2d3dtrace,
     .        i_compress_correct,isstsf,i_wilcox98,i_wilcox98_chiw,
     .        isst2003
      common /axisym/ iaxi2plane,iaxi2planeturb,istrongturbdis,iforcev0
      common /fullns/ ifullns
      common /curvat/ isarc2d,sarccr3,ieasmcc2d,isstrc,sstrc_crc,
     .        isar,crot,isarc3d
      common /mms/ iexact_trunc,iexact_disc,iexact_ring
      common /constit/ i_nonlin,c_nonlin,snonlin_lim,i_tauijs,i_qcr2000,
     .                 i_qcr2013,i_qcr2013v
      common /easmlim/ cmulim
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .        tur10(7)
      common /cgns/ icgns,iccg,ibase,nzones,nsoluse,irind,jrind,krind
      common /writestuff/ ifort50write, j_ifort50write, i_ifort50write
      common /sa_options/ i_saneg,i_sanoft2,sa_cw2,sa_cw3,sa_cv1,sa_ct3,
     .                    sa_ct4,sa_cb1,sa_cb2,sa_sigma,sa_karman
      common /specialtop_kmax1001/ i_specialtop_kmax1001,
     .        a_specialtop_kmax1001,xc_specialtop_kmax1001,
     .        sig_specialtop_kmax1001,vtp_specialtop_kmax1001,
     .        wc_specialtop_kmax1001,fac_specialtop_kmax1001,
     .        cc_specialtop_kmax1001,xerf_specialtop_kmax1001,
     .        sigerf_specialtop_kmax1001
      common /reystressmodel/ issglrrw2012,i_sas_rsm,i_yapterm
      common /iupdate/ iupdatemean
c
      pi = 4.0*atan(1.0)
c
      itoti  = 0
      itotb  = 0
      nfiles = 1
c
c     set error flag to -99 if this routine is called during
c     array sizing step 
c
      if (icall .eq. 0) then
         ierrflg = -99
      else
         ierrflg = -1
      end if
c
c     read keyword-driven input, if any
c
      call readkey(ititr,myid,ibufdim,nbuf,bou,nou,iunit11,
     .             ierrflg)
c
c
c     free-form alphanumeric title

      if (ititr.eq.0) then
         read(iunit5,10)(realval(i),i=1,20)
         do i=1,20
            title(i) = realval(i)
         end do
      end if
   10 format(20a4)   
      write(iunit11,11)(real(title(i)),i=1,20)
   11 format(/1h ,20a4)
c
c     xmach  - freestream Mach number
c     alpha  - angle of attack
c     beta   - side-slip angle
c     reue   - freestream Reynolds number per unit length (millions)
c     tinf   - freestream temperature (degrees Rankine)
c     ialph  - determines how angle of attack is measured in plot3d 
c              type grids (ngrid <0) :
c            = 0 alpha measured in x-z plane (z "up"), i.e. alpha = 90 deg.
c              would give a freestream velocity vector in the +z direction
c            > 0 alpha measured in x-y plane (y "up"), i.e. alpha = 90 deg.
c              would give a freestream velocity vector in the +y direction
c              (note: for cfl3d type grids, alpha is always measured in
c              the x-z plane, z "up")
c     ihstry - determines variables tracked for convergence history:
c            = 0 standard convergence history: residual and cl,cd,cy,cm
c            > 0 control surface history: residual and mass flow,
c              pressure force, viscous force,thrust (momentum) force
c              (forces are resultant forces, i.e. if fx, fy,fz are the force
c              components in the x,y,and z directions, then the resultant
c              force is sqrt(fx**2+fy**2+fz**2); must have ncs > 0
c              and specify which control surfaces are to be tracked -
c              only those surfaces with inorm .ne. 0 are included in the
c              sum)
c
c     ***Note:  isnd and c2spe no longer used - they have been incorporated
c               into solid wall BC2004
c
      read(iunit5,10)
      read(iunit5,*) (realval(i),i=1,5),ialph,ihstry
#   ifdef CMPLX
c
c     complex perturbations are set by keyword input, with defaults = 0 
c
      xmach = cmplx(realval(1),real(xmach_img))
      alpha = cmplx(realval(2),real(alpha_img))
      beta  = cmplx(realval(3),real(beta_img))
      reue  = cmplx(realval(4),real(reue_img))
      tinf  = cmplx(realval(5),real(tinf_img))
#   else
      xmach = realval(1)
      alpha = realval(2)
      beta  = realval(3)
      reue  = realval(4)
      tinf  = realval(5)
#   endif
      write(iunit11,24)
   24 format(6x,4hMach,5x,5halpha,6x,4hbeta,6x,4hReUe,3x,7hTinf,dR,
     .5x,5hialph,4x,6hihstry)
c
      if (real(reue).le.0.e0) reue = 1.0e0
      reue = reue*1.0e+06
      if (real(tinf).le.0.e0) tinf = 460.e0
c
      write(iunit11,20) real(xmach),real(alpha),real(beta),real(reue),
     .                  real(tinf),ialph,ihstry
   20 format(3f10.5,e10.3,f10.5,i10,2i10)      
      if (cltarg.ne.99999.0) then
         write(iunit11,'('' Note: will adjust alpha to reach cltarg ='',
     .                  f12.7)') cltarg
      end if
#   ifdef CMPLX
c
      write(iunit11,'(''    Complex Perturbations:'')')
      if (real(geom_img).eq.0.) then
         write(iunit11,25)
   25    format(4x,6hMach_i,3x,7halpha_i,4x,6hbeta_i,4x,6hReUe_i,1x,
     .   9hTinf_i,dR)
         write(iunit11,'(10e10.3)') real(xmach_img),real(alpha_img),
     .                              real(beta_img),real(reue_img),
     .                              real(tinf_img) 
      else
         write(iunit11,'(4x,6hgeom_i)')
         write(iunit11,'(e10.3)') real(geom_img)
      end if
#   endif
c
      alpha = alpha/radtodeg
      beta  = beta/radtodeg
c
      write(iunit11,16)
   16 format(/,38h note: isnd and c2spe no longer used -,
     .         20h see notes on bc2004,/)
c
c     sref - reference area
c     cref - reference length
c     bref - reference span
c     xmc  - moment center in x-direction
c     ymc  - moment center in y-direction
c     zmc  - moment center in z-direction
c
      read(iunit5,10)
      read(iunit5,*) (realval(i),i=1,6)
      sref = realval(1)
      cref = realval(2)
      bref = realval(3)
      xmc  = realval(4)
      ymc  = realval(5)
      zmc  = realval(6) 
c
      if (real(sref).eq.0.) sref = 1.
      if (real(cref).le.0.) cref = 1.
      if (real(bref).le.0.) bref = 1.
c
      write(iunit11,23)
   23 format(6x,4hsref,6x,4hcref,6x,4hbref,7x,3hxmc,7x,3hymc,
     .7x,3hzmc)
      write(iunit11,22) real(sref),real(cref),real(bref),real(xmc),
     .                  real(ymc),real(zmc)
      if (real(sref).lt.0.) then
         write (iunit11,'(''WARNING: input sref<0, taking abs(sref)'')')
         write (iunit11,'(''  Note: fixed Cl option is no longer '',
     .                    ''triggered with sref<0...use keyword '',
     .                    ''input'')')
         sref = ccabs(sref)
      end if
c
   22 format(f10.3,5f10.5)
c
c     dt      - time step
c             < 0  local time stepping, CFL=abs(dt)
c             > 0  constant time step  (=dt)        
c     irest   = 0   start from freestream conditions
c             > 0   restart from previous solution 
c             < 0   restart from previous solution, but do not save
c                   previous history information
c     iflagts = 0  constant dt
c             > 0  dt ramped over iflagts steps to dt*fmax
c     fmax    - maximum increase in dt
c     iunst   = 0  steady grid (defaults to 0 if dt<0)
c             > 0  dynamic rigid/deforming grid
c               1  rigid mesh
c               2  deforming mesh
c               3  combination rigid and deforming
c     cfltau  - cfl number for tau-ts time step scheme
c
      read(iunit5,10)
      read(iunit5,*) realval(1),irest,iflagts,realval(2),iunst,
     .               realval(3)
      dt     = realval(1)
      dtold  = dt
      fmax   = realval(2)
      cfltau = realval(3)
      if (real(dt) .lt. 0.0)  iunst = 0
      if (real(dt) .lt. 0.0 .and. ides .ge. 1) then
        write(iunit11,'('' WARNING: DES/DDES SHOULD be run time-'',
     .   ''accurately (dt>0)'')')
      end if
      if (real(dt) .gt. 0.0 .and. cltarg.ne.99999.0) then
         write(iunit11,'('' Cannot use cltarg with time-accurate'')')
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
      write(iunit11,31)
   31 format(8x,2hdt,5x,5hirest,3x,7hiflagts,6x,4hfmax,5x,
     .       5hiunst,3x,7hcfl_tau)
      write(iunit11,30) real(dt),irest,iflagts,real(fmax),iunst,
     .                  real(cfltau)
   30 format(f10.5,2i10,f10.5,i10,f10.5)
c
      if (irest .lt. 0) then
        irest2=1
        irest=-irest
      else
        irest2=0
      end if
      if (irest .eq. 0 .and. iteravg .eq. 2) then
        write(iunit11,98)
   98   format(/,47h when irest=0, iteravg=2 starts running-average,
     .   31h from scratch (as if iteravg=1),/)
      end if
c 
      fmax   = ccabs(fmax)
c
c     ngrid   - number of grids input
c             > 0  grid is cfl3d type
c             < 0  grid is plot3d type  
c     nplot3d - number of flowfield data sets to be written in
c               plot3d format
c               if nplot3d = any number < 0, then the plot3d files
c               are automatically set to include all solid surfaces
c               (no field points). if nplot3d < 0, then abs(nplot3d)
c               plot3d data lines MUST appear in the plot3d input
c               plot3d input section below. However, since the only
c               piece of data used from the abs(nplot3d) data lines
c               is iptype, and only the last value is used, then the
c               recommended trigger for this option is nplot3d = -1,
c               with just one data line at the bottom, with the iptype
c               desired (see iptype options below). Note that all other
c               data will be ignored, and so any values may be used
c     nprint  - number of data sets to be sent to an output file
c               if nprint = any number < 0, then the prinout file
c               is automatically set to include all solid surfaces
c               (no field points). if nprint < 0, the same comments
c               apply as made for nplot3d <0 above.
c     nwrest  - number of iterations between updates of the binary
c                restart file
c     ichk    - checks for negative densities and/or pressures in 
c                subroutines gfluxr, hfluxr, ffluxr, and update if 
c                set equal to 1
c     i2d     - 2-d flag
c             = 0  3-d (default)
c             = 1  2-d
c             =-1  2-d with farfield point-vortex correction, 
c                  when bc 1003 used
c                  (Uses xmc & zmc as location of point vortex)
c     ntstep    - number of time steps (defaults to 1 if dt<0)
c     ita     - temporal accuracy
c             = 1  1st order
c             = 2  2nd order
c             =-1  1st order + pseudo-time term for large dt stability
c             =-2  2nd order + pseudo-time term for large dt stability
c                  note: for second-order accuracy in time, subiterations
c                  MUST be used to drive the first order errors from the
c                  approximate factorization (and diagonalization, if used)
c                  to zero.
c
      read(iunit5,10)
      read(iunit5,*) ngrid,nplot3d,nprint,nwrest,ichk,i2d,ntstep,ita
      write(iunit11,1639)
 1639 format(5x,5hngrid,3x,7hnplot3d,4x,6hnprint,4x,6hnwrest,
     .       6x,4hichk,7x,3hi2d,3x,7h ntstep,4x,6h   ita)
      if (real(dt) .lt. 0.0) ita = 1
      if (ita.eq.0 .and. real(dt).gt.0.0) then
         write(iunit11,*)' stopping...must chose abs(ita) > 0 for ',
     .   'time-accurate computations'
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (i2d.lt.-1 .or. i2d.gt.1) i2d = 0
      if (real(dt).lt.0.0) ntstep = 1
      if (nwrest .le. 0) then
        write(iunit11,
     .  '('' NWREST set incorrectly.  Resetting to 2000'')')
        nwrest=2000
      end if
      if (nwrest .lt. 20) then
        write(iunit11,
     .  '('' NWREST is very small... be sure you want this'')')
      end if
      write(iunit11,36) ngrid,nplot3d,nprint,nwrest,ichk,i2d,ntstep,ita
   36 format(8i10)
      iipv=0
      if (i2d .eq. -1) then
        iipv=1
        i2d=1
      end if
c   DES must be run in 3-D
      if (ides .ge. 1 .and. i2d .eq. 1) then
        write(iunit11,'('' WARNING: DES/DDES SHOULD be run 3-D'',
     +    '' (i2d=0)'')')
c        call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c   Check for SARC & EASMCC models
      if ((isarc2d .eq. 1 .or. ieasmcc2d .eq. 1) .and. i2d .eq. 0) then
        write(iunit11,'('' WARNING: SARC & EASMCC currently only'',
     .  '' coded for 2-D curvature in x-z or x-y planes;'')')
        write(iunit11,'(''          (i-index MUST be in spanwise'',
     .  '' direction; i-index derivatives are not accounted for)...'')')
        write(iunit11,'(''          usage in 3-D cases is'',
     .  '' discouraged'')')
      end if
      if (isarc2d .eq. 1 .and. isar .eq. 1) then
        write(iunit11,'('' ERROR: cannot have isarc2d and isar=1'')')
        call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (isarc3d .eq. 1 .and. isar .eq. 1) then
        write(iunit11,'('' ERROR: cannot have isarc3d and isar=1'')')
        call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (isarc3d .eq. 1 .and. isarc2d .eq. 1) then
        write(iunit11,'('' ERROR: cannot have isarc3d and isarc2d=1'')')
        call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
c     for cfl3d type grids, set ialph = 0 regardless of input value
c     note: if it is desired to allow ialph to be used with cfl3d 
c     type grids as well as plot3d type grids then comment out the
c     following 5 lines (which used to read: if (ngrid.gt.0) ialph=0)
c
      if (ngrid .gt. 0 .and. ialph .ne. 0) then
        ialph = 0
        write(iunit11,'(/,'' WARNING: IALPH currently hardwired to 0'',
     .   '' for NGRID>0 - type grids.  See manual.'')')
      end if
c
      ip3dgrd = 0
      if (ngrid.lt.0) then
         ip3dgrd = 1
         ngrid =  abs(ngrid)
      end if
c
         if (ialph .eq. 0 ) then
            uinf = xmach*cos(alpha)*cos(beta)
            winf = xmach*sin(alpha)*cos(beta)
            vinf = -xmach*sin(beta)
         else
            uinf = xmach*cos(alpha)*cos(beta)
            winf = xmach*sin(beta)
            vinf = xmach*sin(alpha)*cos(beta)
            temp=ymc
            ymc=-zmc
            zmc=temp
         end if
         write(iunit11,18) real(uinf),real(vinf),real(winf)
   18    format(/,48h the input values of Mach, alpha, beta and ialph,
     .            28h yield a freestream velocity,/,
     .            23h with x,y,z components:,
     .             6h u0 = ,f8.4,2h, ,6h v0 = ,f8.4,2h, ,
     .             6h w0 = ,f8.4,/,31h           ** make sure this is,
     .            29h consistent with your grid **,/)
c
      nchk = maxgr-ngrid
c
c     check maximum number of grids
c
      if (nchk.lt.0) then
         write(iunit11,1492)
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
 1492 format(55h stopping - insufficient maximum number of grids(maxgr))
c
c     check maximum number of output sets
c
      if (nplot3d.gt.nplots .or. nprint.gt.nplots) then
         write(iunit11,1524)
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
 1524 format(1x,16hnplots too small)
c
c
c     ncg  - number of coarser grids to construct for multigrid/mesh 
c             sequencing  (=0 for embedded mesh)
c     iem  - embedded mesh flag
c          = 0 for global grid
c          = l level of this embedded grid above global grid level
c     iadvance - flag to skip any residual/update calculations
c              >=0  proceed as usual
c              < 0  skip residual/update calculations 
c     iforce  -  flag indicates where to compute forces
c             = 3-digit integer (n1)(n2)(n3)
c                 n1 indicates i-direction
c                 n2 indicates j-direction
c                 n3 indicates k-direction
c               Set n1, n2, and n3 each to 0, 1, 2, or 3
c                 0 = do not compute forces/moments on any faces in this
c                     direction
c                 1 = compute forces/moments on n=1 face in this 
c                     direction
c                 2 = compute forces/moments on n=ndim face in this 
c                     direction
c                 3 = compute forces/moments on both n=1 and n=ndim faces 
c                     in this direction
c                 (EXAMPLES:  000 indicates no force computed on this block
c                             001 indicated force computed at k=1 face 
c                                 on this block, over those cells that have
c                                 solid wall boundary conditions applied
c                                 to them)
c     ivisc(m) - viscous/inviscid interaction flag   m= 1 : I-direction
c              = 0  inviscid                            2 : J-direction
c              = 1  laminar                             3 : K-direction
c              = 2  turbulent - Baldwin-Lomax
c              = 3  turbulent - Baldwin-Lomax with Degani-Schiff option
c              = 4  turbulent - Baldwin-Barth
c              = 5  turbulent - Spalart-Allmaras
c              = 6  turbulent - k-omega
c              = 7  turbulent - SST (Menter's k-omega version)
c              = 8  turbulent - k-omega Explicit Algebraic Stress Model
c                               (EASM-variable-g) in LINEAR formulation
c                               (linear version of #14)
c              = 9  turbulent - k-epsilon EASM-variable-g in LINEAR
c                               formulation
c                               (linear version of #13)
c              =10  turbulent - k-epsilon (Abid version)
c              =11  turbulent - k-epsilon EASM-Gatski/Speziale NONLINEAR
c              =12  turbulent - k-omega EASM-Gatski/Speziale NONLINEAR
c              =13  turbulent - k-epsilon EASM-variable-g NONLINEAR
c              =14  turbulent - k-omega EASM-variable-g NONLINEAR
c              =15  turbulent - k-enstrophy model
c              =16  turbulent - k-kL-MEAH2015 model
c       if ivisc(m) < 0 on input, a wall function is employed.  Should only
c                       (officially) be used for attached flow, when 1st 
c                       gridpoint off wall is not in the sublayer (typically
c                       y+ > 10 or more). The wall function with the
c                       B-L model (ivisc=-2 or -3) can be particularly
c                       non-robust and, although it may sometimes
c                       work well, is not recommended in general.
c    NOTE:  The thin layer viscous terms can be included in either the
c           J-, K-, or I-directions, separately.  The viscous (laminar)
c           terms can be included simultaneously in all three directions.
c
c           For ivisc=2 or 3 (B-L), the viscous (turbulent) terms can be included 
c           simultaneously in, at most, two directions, either J-K or I-K, 
c           for any particular grid.  It is preferable to let K be the primary
c           viscous direction and J be the secondary viscous direction.
c           If I is used as a primary viscous direction, it may be
c           necessary to switch the order of inversions in af3f.f
c           (e.g. from j,k,i to j,i,k).  For multiple-block cases, ivisc
c           should only = 2 on blocks with a body at the lower boundary
c           (K=1, J=1, or I=1).  On other blocks (hopefully relatively far
c           away from bodies), ivisc should be set =0.
c
c           For ivisc>=4 (1- and 2-equation field models), the viscous
c           (turbulent) terms can be included simultaneously in 1, 2, or
c           all 3 directions.  These models can all be used for multiple
c           blocks.  All blocks must have the same ivisc number in at least
c           one of the directions.
c
c***** repeat above for each grid, 1 through ngrid *****
c
      read(iunit5,10)
      nbl = 0
      nblock = 0
      write(iunit11,61)
   61 format(7x,3hncg,7x,3hiem,2x,8hiadvance,4x,6hiforce,
     .       2x,8hivisc(i),2x,8hivisc(j),2x,8hivisc(k))
      iemtot = 0
      do 71 igrid=1,ngrid
      nbl    = nbl+1
      read(iunit5,*) ncg,iem,iad,ifor,ivisc
      write(iunit11,36) ncg,iem,iad,ifor,ivisc
      nblock = nblock+ncg+1
      if (iem .gt. 0 .and. ncg .gt. 0) then
        write(iunit11,
     .  '('' Embedded grid must have ncg=0 (it uses grid'')')
        write(iunit11,
     .  '('' in which it is embedded for coarser levels)'')')
        call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (ncg .gt. 4) then
        write(iunit11,
     .  '('' known bug: ncg must not exceed 4'')')
        call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
      if(ivisc(1).eq.-2 .or. ivisc(2).eq.-2 .or. ivisc(3).eq.-2) then
        write(iunit11,'('' NOTICE:  ivisc=-2 is no longer B-L with '',
     .  ''Degani-Schiff'')')
        write(iunit11,
     .  '(''     ivisc < 0 now indicates Wall Function'')')
        write(iunit11,
     .  '(''     Use ivisc=3 for B-L with Degani-Schiff'')')
      end if
c
      if (ivisc(1).eq.-16 .or. ivisc(2).eq.-16 .or. ivisc(3).eq.-16 .or.
     .    ivisc(1).eq.-30 .or. ivisc(2).eq.-30 .or. ivisc(3).eq.-30 .or.
     .    ivisc(1).eq.-40 .or. ivisc(2).eq.-40 .or. ivisc(3).eq.-40)
     .    then
        write(iunit11,'('' cannot use wall fns with ivisc 16,30,40'')')
        call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c   Set wall iwf, wall function parameter:
      iwf(1)=0
      iwf(2)=0
      iwf(3)=0
      if (ivisc(1) .lt. 0) then
        iwf(1)=1
        ivisc(1) = abs(ivisc(1))
      end if
      if (ivisc(2) .lt. 0) then
        iwf(2)=1
        ivisc(2) = abs(ivisc(2))
      end if
      if (ivisc(3) .lt. 0) then
        iwf(3)=1
        ivisc(3) = abs(ivisc(3))
      end if
c   DES currently only works with models 5,6,7
      if (ides .eq. 1 .or. ides .eq. 2 .or. ides .eq. 3) then
        if (i2d .ne. 1) then
          if ((ivisc(1).eq.5 .or. ivisc(1).eq.6 .or. ivisc(1).eq.7 .or.
     .         ivisc(1).eq.40) .and.
     .        (ivisc(2).eq.5 .or. ivisc(2).eq.6 .or. ivisc(2).eq.7 .or.
     .         ivisc(2).eq.40) .and.
     .        (ivisc(3).eq.5 .or. ivisc(3).eq.6 .or. ivisc(3).eq.7 .or.
     .         ivisc(3).eq.40)) then
            continue
          else
            write(iunit11,'('' ides=1,2,3 currently works with ivisc'',
     .       '' 5,6,7,40 only, and must be on in all directions'')')
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
          end if
        else
          if ((ivisc(2).eq.5 .or. ivisc(2).eq.6 .or. ivisc(2).eq.7 .or.
     .         ivisc(2).eq.40) .and.
     .        (ivisc(3).eq.5 .or. ivisc(3).eq.6 .or. ivisc(3).eq.7 .or.
     .         ivisc(3).eq.40)) then
            continue
          else
            write(iunit11,'('' ides=1,2,3 currently works with ivisc'',
     .       '' 5,6,7,40 only, and must be on in all directions'')')
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
          end if
        end if
      end if
      if (ides .gt. 3) then
        if (i2d .ne. 1) then
          if ((ivisc(1).lt.5 .or. ivisc(1).gt.5) .and.
     .        (ivisc(2).lt.5 .or. ivisc(2).gt.5) .and.
     .        (ivisc(3).lt.5 .or. ivisc(3).gt.5)) then
            write(iunit11,'('' ides>3 currently works with ivisc'',
     .       '' 5 only, and must be on in all directions'')')
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
          end if
        else
          if ((ivisc(2).lt.5 .or. ivisc(2).gt.5) .and.
     .        (ivisc(3).lt.5 .or. ivisc(3).gt.5)) then
            write(iunit11,'('' ides>3 currently works with ivisc'',
     .       '' 5 only, and must be on in all directions'')')
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
          end if
        end if
      end if
c   for Degani schiff, set ideg=1 and reset ivisc=2 within code:
      ideg(1) = 0
      ideg(2) = 0
      ideg(3) = 0
      if (ivisc(1).eq.3) then
        ideg(1)=1
        ivisc(1)=2
      end if
      if (ivisc(2).eq.3) then
        ideg(2)=1
        ivisc(2)=2
      end if
      if (ivisc(3).eq.3) then
        ideg(3)=1
        ivisc(3)=2
      end if
      iv(igrid) = max(ivisc(1),ivisc(2))
      iv(igrid) = max(ivisc(3),iv(igrid))
      if (ivisc(1).gt.16.or. ivisc(2).gt.16.or. ivisc(3).gt.16) then
        if (ivisc(1).ne.25.and.ivisc(2).ne.25.and.ivisc(3).ne.25.and.
     .      ivisc(1).ne.30.and.ivisc(2).ne.30.and.ivisc(3).ne.30.and.
     .      ivisc(1).ne.40.and.ivisc(2).ne.40.and.ivisc(3).ne.40.and.
     .      ivisc(1).ne.72.and.ivisc(2).ne.72.and.ivisc(3).ne.72) then
          write(iunit11,'('' Unknown turb model type chosen'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
        end if
      end if
      if (ivisc(1).ge.4 .or. ivisc(2).ge.4 .or. ivisc(3).ge.4) then
        isminc=1
      else
        isminc=0
      end if
      if (ivisc(1).gt.1 .and. ivisc(2).gt.1 .and. 
     .    ivisc(1).ne.ivisc(2)) then
        write(iunit11,'('' cannot mix turbulence model types'')')
        call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (ivisc(1).gt.1 .and. ivisc(3).gt.1 .and. 
     .    ivisc(1).ne.ivisc(3)) then
        write(iunit11,'('' cannot mix turbulence model types'')')
        call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (ivisc(2).gt.1 .and. ivisc(3).gt.1 .and. 
     .    ivisc(2).ne.ivisc(3)) then
        write(iunit11,'('' cannot mix turbulence model types'')')
        call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
      if (ivisc(2).eq.2 .and. ivisc(1).eq.2) then
         write(iunit11,91)
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
   91 format(53h B-L progr. only for turbulent flow in two directions,
     .34h either K-J or K-I directions only)
c
      iemtot      = iemtot+iem
      ncgg(igrid) = ncg
      if (igrid.eq.1) then
         ncgmax   = ncg
      else
         ncgmax   = max(ncgmax,ncg)
      end if      
c
      iemg(igrid)   = iem
      iadvance(nbl) = iad
      iforce(nbl)   = ifor
      iviscg(nbl,1) = ivisc(1)
      iviscg(nbl,2) = ivisc(2)
      iviscg(nbl,3) = ivisc(3)
      idegg(nbl,1)  = ideg(1)
      idegg(nbl,2)  = ideg(2)
      idegg(nbl,3)  = ideg(3)
      iwfg(nbl,1)   = iwf(1)
      iwfg(nbl,2)   = iwf(2)
      iwfg(nbl,3)   = iwf(3)
      if (ncg.gt.0) then
         do 68 n=1,ncg
         nbl           = nbl+1
         iadvance(nbl) = iadvance(nbl-1)
         iforce(nbl)   = iforce(nbl-1)
         iviscg(nbl,1) = ivisc(1)
         iviscg(nbl,2) = ivisc(2)
         iviscg(nbl,3) = ivisc(3)
         idegg(nbl,1)  = ideg(1)
         idegg(nbl,2)  = ideg(2)
         idegg(nbl,3)  = ideg(3)
         iwfg(nbl,1)   = iwf(1)
         iwfg(nbl,2)   = iwf(2)
         iwfg(nbl,3)   = iwf(3)
   68    continue
      end if
   71 continue
c
c     check maximum number of blocks
c
      nchk   = maxbl-nblock
      if (nchk.lt.0) then
         write(iunit11,1649)
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
 1649 format(42h stopping - insufficient maximum number of,
     .14h blocks(maxbl))
c
      ivmx=iv(1)
      if (ngrid .gt. 1) then
        do 73 igrid=2,ngrid
          if (iv(igrid).gt.ivmx) ivmx=iv(igrid)
   73   continue
        if (ivmx.gt.3) then
          do 72 igrid=2,ngrid
            if (iv(igrid).ne.iv(igrid-1)) then
              write(iunit11,
     .        '('' for turb models > 3, need same number'',
     .        '' on ALL blocks'')')
              call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
            end if
   72     continue
        end if
      end if
c
      if ((isarc2d .eq. 1 .or. isarc3d .eq. 1) .and. ivmx .ne. 5) then
        write(iunit11,'('' WARNING: isarc2d/3d not active'',
     .  '' unless SA model (ivisc=5) used'')')
      end if
      if (ieasmcc2d .eq. 1 .and. (ivmx.ne.8 .and. ivmx.ne.9 .and.
     .  ivmx.ne.11 .and. ivmx.ne.12 .and. ivmx.ne.13 .and.
     .  ivmx.ne.14)) then
        write(iunit11,'('' WARNING: ieasmcc2d not active'',
     .  '' unless EASM model (ivisc=8,9,11,12,13, or 14) used'')')
      end if
c
      if (i_wilcox06 .eq. 1 .and. ivmx .ne. 6) then
        write(iunit11,'('' WARNING: i_wilcox06 not active'',
     .  '' unless Wilcox k-o model (ivisc=6) used'')')
      end if
      if (i_wilcox98 .eq. 1 .and. ivmx .ne. 6) then
        write(iunit11,'('' WARNING: i_wilcox98 not active'',
     .  '' unless Wilcox k-o model (ivisc=6) used'')')
      end if
      if (i_wilcox98 .eq. 1 .and. i_wilcox06 .eq. 1) then
        write(iunit11,'('' Error: cannot have both i_wilcox98=1'',
     .    '' and i_wilcox06=1'')')
        call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
      if (i_turbprod_kterm .eq. 1 .and. (ivmx .ne. 6 .and.
     .  ivmx .ne. 7)) then
        if (ikoprod .ne. 1) then
        write(iunit11,'('' WARNING: i_turbprod_kterm not active'',
     .  '' unless ivisc=6 or 7 used AND ikoprod=1'')')
        end if
      end if
c
      if (i_catris_kw .eq. 1 .and. (ivmx.ne.6 .and. ivmx.ne.7 .and.
     .  ivmx.ne.8 .and. ivmx.ne.12 .and. ivmx.ne.14)) then
        write(iunit11,'('' WARNING: i_catris_kw not active'',
     .  '' unless ivisc=6,7,8,12, or 14 used'')')
      end if
c
      if (isstrc .gt. 0 .and. (ivmx.ne.6 .and. ivmx.ne.7)) then
        write(iunit11,'('' WARNING: isstrc not active'',
     .  '' unless ivisc=6 or 7'')')
      end if
      if (isstsf .eq. 1 .and. (ivmx.ne.6 .and. ivmx.ne.7)) then
        write(iunit11,'('' WARNING: isstsf not active'',
     .  '' unless ivisc=6 or 7'')')
      end if
c
      if (isst2003 .eq. 1 .and. ivmx .ne. 7) then
        write(iunit11,'('' ERROR: isst2003 must have ivisc=7'')')
        call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
      if (abs(prod2d3dtrace-0.5) .lt. 0.01) then
        prod2d3dtrace=0.5
      else if (abs(prod2d3dtrace-0.33333333) .lt. 0.01) then
        prod2d3dtrace=0.33333333
      else
        prod2d3dtrace=0.
      end if
c
      if (ikoprod .gt. 0 .and. (ivmx.ne.6 .and. ivmx.ne.7 .and.
     .  ivmx.ne.10 .and. ivmx.ne.15 .and. ivmx.ne.16 .and.
     .  ivmx.ne.30 .and. ivmx.ne.40)) then
        write(iunit11,'('' WARNING: ikoprod not active'',
     .  '' unless ivisc=6,7,10,15,16,30, or 40 used'')')
      end if
c
      if (iturbprod .gt. 0 .and. (ivmx.ne.8 .and. ivmx.ne.9 .and.
     .  ivmx.ne.13 .and. ivmx.ne.14)) then
        write(iunit11,'('' WARNING: iturbprod not active'',
     .  '' unless ivisc=8,9,13, or 14 used'')')
      end if
c
      if (isstdenom .ne. 0 .and. ivmx .ne. 7) then
        write(iunit11,'('' WARNING: isstdenom not active'',
     .  '' unless ivisc=7 used'')')
      else if (isstdenom .ne. 0 .and. i_bsl .ne. 0) then
        write(iunit11,'('' WARNING: isstdenom not active'',
     .  '' when i_bsl=1 (BSL option in effect)'')')
      end if
c
      if (ibeta8kzeta .ne. 0 .and. ivmx .ne. 15) then
        write(iunit11,'('' WARNING: ibeta8kzeta not active'',
     .  '' unless ivisc=15 used'')')
      end if
c
      if (i_bsl .ne. 0 .and. ivmx .ne. 7) then
        write(iunit11,'('' WARNING: i_bsl not active'',
     .  '' unless ivisc=7 used'')')
      end if
c
      if (keepambient .ne. 0 .and. (ivmx.ne.6 .and. ivmx.ne.7 .and.
     .  ivmx.ne.8 .and. ivmx.ne.9 .and. ivmx.ne.10 .and.
     .  ivmx.ne.11 .and. ivmx.ne.12 .and. ivmx.ne.13 .and.
     .  ivmx.ne.14 .and. ivmx.ne.30 .and. ivmx.ne.40)) then
        write(iunit11,'('' WARNING: keepambient not active'',
     .  '' unless ivisc=6,7,8,9,10,11,12,13,14,30, or 40 used'')')
      end if
c
      if (ivmx .eq. 25 .or. ivmx .eq. 30 .or. ivmx .eq. 40 .or.
     .    ivmx .eq. 72 .or. ivmx .eq. 16) then
        write(iunit11,'('' WARNING: ivisc='',i4,'' still under'',
     .    '' development... use at your own risk!'')') ivmx
      end if
c
c     Need to check for keywords that do not work for ivmx=72:
      if (ivmx .eq. 72) then
        write(iunit11,'('' WARNING: keyword ifullns has no effect for'',
     .    '' ivisc=72 - it is ALWAYS full NS'')')
        write(iunit11,'(''          (but diffusion'',
     .    '' terms in turb eqns are still thin-layer)'')')
        if (ifullns .ne. 0) then
          write(iunit11,'('' ERROR: ifullns must be 0 for ivisc=72'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
        end if
        if (edvislim .ne. 1.e10) then
          write(iunit11,'('' WARNING: keyword edvislim has influence'',
     .      '' only in turb heat flux in energy eqn for ivisc=72'')')
        end if
        if (icgns .ne. 0) then
          write(iunit11,'('' ERROR: CGNS not working yet'',
     .      '' for ivisc=72'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
        end if
        if (itaturb .ne. 1) then
          write(iunit11,'('' ERROR: itaturb not implemented'',
     .      '' for ivisc=72'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
        end if
        if (iturbord .ne. 1) then
          write(iunit11,'('' ERROR: iturbord not implemented'',
     .      '' for ivisc=72 (default is 1st order)'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
        end if
c       if (nfreeze .ne. 0) then
c         write(iunit11,'('' ERROR: nfreeze not implemented'',
c    .      '' for ivisc=72'')')
c         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
c       end if
c       if (iaxi2planeturb .ne. 0) then
c         write(iunit11,'('' ERROR: iaxi2planeturb not implemented'',
c    .      '' for ivisc=72'')')
c         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
c       end if
        if (istrongturbdis .ne. 0) then
          write(iunit11,'('' ERROR: istrongturbdis not implemented'',
     .      '' for ivisc=72'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
        end if
c       if (cflturb(2) .ne. 0. .or. cflturb(3) .ne. 0. .or.
c    .      cflturb(4) .ne. 0. .or. cflturb(5) .ne. 0. .or.
c    .      cflturb(6) .ne. 0. .or. cflturb(7) .ne. 0.) then
c         write(iunit11,'('' ERROR: cflturb2...cflturb7 not'',
c    .      '' implemented for ivisc=72'')')
c         write(iunit11,'(''  ... use cflturb1 to set CFL for all'',
c    .      '' seven turb equations'')')
c         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
c       end if
        if (tur1cut .ne. 1.e-20 .or. tur2cut .ne. 1.e-20) then
          write(iunit11,'('' ERROR: tur1cut, tur2cut not implemented'',
     .      '' for ivisc=72'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
        end if
        if (tur1cutlev .ne. 0. .or. tur2cutlev .ne. 0.) then
          write(iunit11,'('' ERROR: tur1cutlev, tur2cutlev not'',
     .      '' implemented for ivisc=72'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
        end if
        if (i_lam_forcezero .ne. 0) then
          write(iunit11,'('' ERROR: i_lam_forcezero not'',
     .      '' implemented for ivisc=72'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
        end if
      end if
c
#   ifdef CMPLX
      if (ivmx==72) then
        write(iunit11,'('' Error... ivmx=72 not allowed for'',
     +   '' complex mode'')')
        call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
#   endif
c
      if (i_nonlin .ne. 0 .and. (ivmx .lt. 2 .or. ivmx .eq. 8
     .  .or. ivmx .eq. 9 .or. ivmx .eq. 11 .or. ivmx .eq. 12 .or.
     .  ivmx .eq. 13 .or. ivmx .eq. 14 .or. ivmx .gt. 15)) then
        write(iunit11,'('' i_nonlin cannot be used with'',
     .    '' this ivisc'')')
        call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
      if (i_tauijs .eq. 1) then
        if (ivmx .lt. 70 .and. ngrid .ne. 1) then
          write(iunit11,'('' ERROR: across-boundary BCs not correct'',
     .      '' for i_tauijs=1, ivmx<70;'',
     .      '' cannot use multiple zones currently'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
        end if
        if (ifullns .ne. 0) then
          write(iunit11,'('' ERROR: ifullns must be 0 when'',
     .      '' i_tauijs=1'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
        end if
        if (ivmx .eq. 11 .or. ivmx .eq. 12 .or. ivmx .eq. 13 .or.
     .      ivmx .eq. 14) then
          write(iunit11,'('' i_tauijs=1 cannot be used with'',
     .    '' ivisc=11,12,13,14'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
        end if
        if (i_nonlin .ne. 0) then
          write(iunit11,'('' i_tauijs=1 cannot be used with'',
     .    '' i_nonlin'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
        end if
      end if
c
      if (i_qcr2000 .eq. 1 .and. i_tauijs .ne. 1) then
        write(iunit11,'('' i_qcr2000 must be used with'',
     .    '' i_tauijs=1'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (i_qcr2013 .eq. 1 .and. i_tauijs .ne. 1) then
        write(iunit11,'('' i_qcr2013 must be used with'',
     .    '' i_tauijs=1'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (i_qcr2013v .eq. 1 .and. i_tauijs .ne. 1) then
        write(iunit11,'('' i_qcr2013v must be used with'',
     .    '' i_tauijs=1'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (i_qcr2000 .eq. 1 .and. i_qcr2013 .eq. 1) then
        write(iunit11,'('' cannot have both i_qcr2000=1'',
     .    '' and i_qcr2013=1'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (i_qcr2000 .eq. 1 .and. i_qcr2013v .eq. 1) then
        write(iunit11,'('' cannot have both i_qcr2000=1'',
     .    '' and i_qcr2013v=1'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (i_qcr2013 .eq. 1 .and. i_qcr2013v .eq. 1) then
        write(iunit11,'('' cannot have both i_qcr2013=1'',
     .    '' and i_qcr2013v=1'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (i_qcr2013 .eq. 1 .and. ivmx .ne. 5) then
        write(iunit11,'('' i_qcr2013 is only usable with SA model'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (i_qcr2013v .eq. 1 .and. ivmx .ne. 5) then
        write(iunit11,'('' i_qcr2013v is only usable with SA model'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
      if ((ivmx .eq. 8 .or. ivmx .eq. 9 .or. ivmx .eq. 13 .or.
     .     ivmx .eq. 14) .and. real(cmulim) .gt. .04) then
        write(iunit11,'('' error... cmulim set too high'')')
        call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if 
c
c     grid dimensions:
c     idim - number of points in I-direction
c
c     jdim - number of points in J-direction
c
c     kdim - number of points in K-direction
c
c***** repeat above for each grid, 1 through ngrid *****
c
      read(iunit5,10)
      nbl = 0
      write(iunit11,1631)
 1631 format(6x,4hidim,6x,4hjdim,6x,4hkdim)
      do 7001 igrid=1,ngrid
      nbl = nbl+1
      read(iunit5,*) idim,jdim,kdim
      write(iunit11,336) idim,jdim,kdim
  336 format(3i10)
      if (i2d.gt.0 .and. idim.ne.2) then
         write(iunit11,337)
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
  337    format(42h stopping...must have idim=2 for 2d cases!)
      end if
      if (i2d .eq. 0 .and. idim .eq. 2) then
         write(iunit11,338)
  338    format(49h WARNING: idim=2 for i2d=0 may not be a good idea)
      end if
      if (iaxi2plane.eq.1 .and. (i2d .ne. 0 .or. idim .ne. 2)) then
        write(iunit11,'('' WARNING: iaxi2plane set back to 0'')')
        write(iunit11,'(''    (i2d must=0 and idim must=2)'')')
        iaxi2plane=0
      end if
      if (iaxi2planeturb.eq.1 .and. (i2d .ne. 0 .or. idim .ne. 2)) then
        write(iunit11,'('' WARNING: iaxi2planeturb set back to 0'')')
        write(iunit11,'(''    (i2d must=0 and idim must=2)'')')
        iaxi2planeturb=0
      end if
c   check usage of full N-S
      if (ifullns .ne. 0) then
        write(iunit11,'('' Note that even with ifullns=1, diffusion'',
     .    '' terms in turb eqns are still thin-layer)'')')
        if (i2d .eq. 0 .and. idim .gt. 2 .and. (ivisc(1).eq.0 .or. 
     .      ivisc(2).eq.0 .or. ivisc(3).eq.0)) then
          write(iunit11,'('' WARNING: When using full N-S in 3-D, '',
     .      ''viscous terms (ivisc) usually '')')
          write(iunit11,'(''must be ON in all 3 directions'')')
        end if
        if (i2d .eq. 0 .and. idim .eq. 2 .and. (ivisc(2).eq.0 .or. 
     .    ivisc(3).eq.0)) then
          write(iunit11,'('' ERROR: When using full N-S in 3-D, '',
     .      ''with idim=2, viscous terms (ivisc) must be'')')
          write(iunit11,'('' ON in both J and K directions'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
        end if
        if (i2d .ne. 0 .and. (ivisc(2).eq.0 .or. ivisc(3).eq.0)) then
          write(iunit11,'('' ERROR: When using full N-S in 2-D, '',
     .      ''viscous terms (ivisc) must be'')')
          write(iunit11,'('' ON in both J and K directions'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
        end if
        if (idim .eq. 2 .and. ivisc(1).ne.0) then
          write(iunit11,'('' ERROR: Having ifullns .ne. 0 with'',
     .      '' idim=2 and ivisc(1) .ne. 0'')')
          write(iunit11,'('' can cause major problems; recommend'',
     .      '' setting ivisc(1)=0'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
        end if
      end if
      if (ides .ge. 1 .and. idim .eq. 2) then
         write(iunit11,339)
  339    format(52h WARNING: idim=2 for DES/DDES may not be a good idea)
      end if
      nblg(igrid) = nbl
      idimg(nbl)  = idim
      jdimg(nbl)  = jdim
      kdimg(nbl)  = kdim
      if (jdim .eq. 2 .or. kdim .eq. 2) then
        write(iunit11,'('' Error.  Neither JDIM nor KDIM can be 2.'')')
        write(iunit11,
     .  '('' (can cause problems when need 2 ghost BCs)'')')
        call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
      ncg = ncgg(igrid)
      if (ncg.gt.0) then
         do 6885 n=1,ncg
         nbl        = nbl+1
         idimg(nbl) = idimg(nbl-1)/2+1
         jdimg(nbl) = jdimg(nbl-1)/2+1
         kdimg(nbl) = kdimg(nbl-1)/2+1
         if (idim.eq.2) then
c        2-d meshes
            idimg(nbl) = 2
         end if
         if (jdimg(nbl) .eq. 2 .or. kdimg(nbl) .eq. 2) then
           write(iunit11,'('' Error.  Coarser levels of JDIM and KDIM'',
     .      '' cannot be 2.'')')
           write(iunit11,'('' (can cause problems when need 2 ghost'',
     .      '' BCs)'')')
           call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
         istop=0
         if (float(idimg(nbl-1)/2) .eq. float(idimg(nbl-1))/2. .and.
     .    idim .gt. 2) then
           write(iunit11,
     .     '('' Cannot create coarser level for idim past'',
     .     i6)') idimg(nbl-1)
           istop=1
         end if
         if (float(jdimg(nbl-1)/2) .eq. float(jdimg(nbl-1))/2.) then
           write(iunit11,
     .     '('' Cannot create coarser level for jdim past'',
     .     i6)') jdimg(nbl-1)
           istop=1
         end if
         if (float(kdimg(nbl-1)/2) .eq. float(kdimg(nbl-1))/2.) then
           write(iunit11,
     .     '('' Cannot create coarser level for kdim past'',
     .     i6)') kdimg(nbl-1)
           istop=1
         end if
         if (istop .eq. 1) then
           call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
 6885    continue
      end if
 7001 continue
c
c     Indices that define a region of laminar flow within each block 
c     (if want to simulate transition) for turbulent flow computations
c
c     ilamlo - lower i index (0 if no laminar region)
c     ilamhi - upper i index (0 if no laminar region)
c     jlamlo - lower j index (0 if no laminar region)
c     jlamhi - upper j index (0 if no laminar region)
c     klamlo - lower k index (0 if no laminar region)
c     klamhi - upper k index (0 if no laminar region)
c
c***** repeat above for each grid, 1 through ngrid *****
c
      read(iunit5,10)
      nbl = 0
      write(iunit11,1939)
 1939 format(4x,6hilamlo,4x,6hilamhi,4x,6hjlamlo,
     .       4x,6hjlamhi,4x,6hklamlo,4x,6hklamhi)
      do 7201 igrid=1,ngrid
      nbl = nbl+1
      read(iunit5,*) ilamlo,ilamhi,jlamlo,jlamhi,klamlo,klamhi
      if (ilamhi .gt. idimg(nbl)) ilamhi=idimg(nbl)
      if (jlamhi .gt. jdimg(nbl)) jlamhi=jdimg(nbl)
      if (klamhi .gt. kdimg(nbl)) klamhi=kdimg(nbl)
      if (ilamlo .le. 0) then
         ilamlo = 0
         ilamhi = 0
      end if
      if (ilamhi .le. 0) then
         ilamlo = 0
         ilamhi = 0
      end if
      if (jlamlo .le. 0) then
         jlamlo = 0
         jlamhi = 0
      end if
      if (jlamhi .le. 0) then
         jlamlo = 0
         jlamhi = 0
      end if
      if (klamlo .le. 0) then
         klamlo = 0
         klamhi = 0
      end if
      if (klamhi .le. 0) then
         klamlo = 0
         klamhi = 0
      end if
      write(iunit11,36) ilamlo,ilamhi,jlamlo,jlamhi,klamlo,klamhi
      ilamlog(nbl)  = ilamlo
      ilamhig(nbl)  = ilamhi
      jlamlog(nbl)  = jlamlo
      jlamhig(nbl)  = jlamhi
      klamlog(nbl)  = klamlo
      klamhig(nbl)  = klamhi
c
      ncg = ncgg(igrid)
      if (ncg.gt.0) then
         do 6285 n=1,ncg
         nbl        = nbl+1
         ilamlog(nbl) = ilamlog(nbl-1)/2+1
         ilamhig(nbl) = ilamhig(nbl-1)/2+1
         jlamlog(nbl) = jlamlog(nbl-1)/2+1
         jlamhig(nbl) = jlamhig(nbl-1)/2+1
         klamlog(nbl) = klamlog(nbl-1)/2+1
         klamhig(nbl) = klamhig(nbl-1)/2+1
         if (i2d .eq. 1) then
c        2-d meshes
            ilamlog(nbl) = ilamlog(nbl-1)
            ilamhig(nbl) = ilamhig(nbl-1)
         end if
         if (ilamlog(nbl-1) .eq. 0) then
            ilamlog(nbl)=0
            ilamhig(nbl)=0
         end if
         if (jlamlog(nbl-1) .eq. 0) then
            jlamlog(nbl)=0
            jlamhig(nbl)=0
         end if
         if (klamlog(nbl-1) .eq. 0) then
            klamlog(nbl)=0
            klamhig(nbl)=0
         end if
 6285    continue
      end if
 7201 continue
c
c     inewg - restart flag for grid (not needed if irest=0)
c           = 0  read flowfield data from restart file
c           = 1  initialize at freestream or by linear interpolation
c                 from coarser grids
c     igridc - grid to which this grid connects (input 0 for global
c              mesh(iem=0) and the grid number in which the embedded
c              mesh fits for embedded meshes(iem>0))
c     js,ks,is - starting indices in connecting grid for placement of 
c                embedded mesh (input 0 for global meshes)
c     je,ke,ie - ending indices in connecting grid for placement of 
c                embedded mesh (input 0 for global meshes)
c     NOTE:  The embedded meshes must be a regular refinement in all 
c            directions of the grid to which it connects.
c
c***** repeat above for each grid, 1 through ngrid *****
c
      read(iunit5,10)
      nbl = 0
      write(iunit11,1931)
 1931 format(5x,5hinewg,4x,6higridc,
     .       8x,2his,8x,2hjs,8x,2hks,8x,2hie,8x,2hje,8x,2hke)
      do 7002 igrid=1,ngrid
      nbl = nbl+1
      read(iunit5,*) inewg,igridc,is,js,ks,ie,je,ke
      write(iunit11,36) inewg,igridc,is,js,ks,ie,je,ke
      if (igridc .gt. 0) then
         nblcg(nbl) = nblg(igridc)
      else
         nblcg(nbl) = nbl
      end if
      inewgg(igrid) = inewg
      jsg(nbl)      = js
      ksg(nbl)      = ks
      isg(nbl)      = is
      jeg(nbl)      = je
      keg(nbl)      = ke
      ieg(nbl)      = ie
c
      ncg = ncgg(igrid)
      if (ncg.gt.0) then
         do 6887 n=1,ncg
         nbl        = nbl+1
         nblcg(nbl) = nblcg(nbl-1)+1
         jsg(nbl)   = jsg(nbl-1)/2+1
         ksg(nbl)   = ksg(nbl-1)/2+1
         isg(nbl)   = isg(nbl-1)/2+1
         jeg(nbl)   = jeg(nbl-1)/2+1
         keg(nbl)   = keg(nbl-1)/2+1
         ieg(nbl)   = ieg(nbl-1)/2+1
 6887    continue
      end if
 7002 continue
c
c
c     idiag(m) - matrix inversion flag
c              = 0  5x5 block tridiagonal inversion
c              = 1  scalar tridiagonal inversions (recommended)
c     iflim(m) - flux limiter flag                     m=1 : I-direction
c              = 0  unlimited                           =2 : J-direction
c              = 1  smooth limiter                      =3 : K-direction
c              = 2  min-mod scheme (recommended)
c              = 3  smooth limiter tuned to Kappa = 1/3 with a cut-off
c                   to eliminate limiting in regions of small gradients
c              = 4  a corrected version of 3 that is invariant with
c                   block splittng
c
c***** repeat above for each grid, 1 through ngrid *****
c
      read(iunit5,10)
      nbl = 0
      write(iunit11,9004)
 9004 format(2x,8hidiag(i),2x,8hidiag(j),2x,8hidiag(k),
     .       2x,8hiflim(i),2x,8hiflim(j),2x,8hiflim(k))
      do 9002 igrid=1,ngrid
      nbl = nbl+1
      read(iunit5,*) idiag,iflim
      write(iunit11,36) idiag,iflim
c
c  limit on ngc when iflim=4 (set via dimension of ncell):
c
      if ((iflim(1).eq.4 .or. iflim(2).eq.4 .or.
     .   iflim(3).eq.4) .and. ncg.ge.10) then
         write(iunit11,'('' currently ncg limited < 10 for iflim=4'')')
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
      idiagg(nbl,1) = idiag(1)
      idiagg(nbl,2) = idiag(2)
      idiagg(nbl,3) = idiag(3)
      iflimg(nbl,1) = iflim(1)
      iflimg(nbl,2) = iflim(2)
      iflimg(nbl,3) = iflim(3)
c
      ncg = ncgg(igrid)
      if (ncg.gt.0) then
         do 6889 n=1,ncg
         nbl           = nbl+1
         idiagg(nbl,1) = idiag(1)
         idiagg(nbl,2) = idiag(2)
         idiagg(nbl,3) = idiag(3)
         iflimg(nbl,1) = iflim(1)
         iflimg(nbl,2) = iflim(2)
         iflimg(nbl,3) = iflim(3)
 6889    continue
      end if
 9002 continue
c
c
c     ifds(m)  - spatial differencing parameter for Euler fluxes
c              = 0  flux-vector splitting
c              = 1  flux-difference splitting (Roe's scheme) 
c                   (recommended)
c     rkap0(m) - spatial differencing parameter for Euler fluxes
c              =  -1  fully upwind
c              =   0  Frommes's scheme
c              =   1  central
c              = 1/3  upwind-biased third order
c                     (recommended)
c
c***** repeat above for each grid, 1 through ngrid *****
c
      read(iunit5,10)
      nbl = 0
      write(iunit11,1951)
 1951 format(3x,7hifds(i),3x,7hifds(j),3x,7hifds(k),
     .       2x,8hrkap0(i),2x,8hrkap0(j),2x,8hrkap0(k))
      do 7012 igrid=1,ngrid
      nbl = nbl+1
      read(iunit5,*) ifds(1),ifds(2),ifds(3),realval(1),realval(2),
     .               realval(3)
      rkap0(1) = realval(1)
      rkap0(2) = realval(2)
      rkap0(3) = realval(3)
      write(iunit11,1950) ifds(1),ifds(2),ifds(3),real(rkap0(1)),
     .                    real(rkap0(2)),real(rkap0(3))
 1950 format(3i10,3f10.4)
      if (iaxi2plane.eq.1 .and. (ifds(1).eq.1.or.ifds(2).eq.1.or.
     .    ifds(3).eq.1)) then
        write(iunit11,'('' WARNING: if problems occur when running'',
     .    '' two-plane axisymmetric case with'')')
        write(iunit11,'(''   singular or near-singular axis,'',
     .    '' the use of IFDS=0 may help'')')
      end if
      ifdsg(nbl,1)  = ifds(1)
      ifdsg(nbl,2)  = ifds(2)
      ifdsg(nbl,3)  = ifds(3)
      rkap0g(nbl,1) = rkap0(1)
      rkap0g(nbl,2) = rkap0(2)
      rkap0g(nbl,3) = rkap0(3)
c
      ncg = ncgg(igrid)
      if (ncg.gt.0) then
         do 6987 n=1,ncg
         nbl           = nbl+1
         ifdsg(nbl,1)  = ifds(1)
         ifdsg(nbl,2)  = ifds(2)
         ifdsg(nbl,3)  = ifds(3)
         rkap0g(nbl,1) = rkap0(1)
         rkap0g(nbl,2) = rkap0(2)
         rkap0g(nbl,3) = rkap0(3)
 6987    continue
      end if
 7012 continue
c
c
c     grid      - grid number (must be in order)
c     nbci0     - number of segments on i=0 boundary
c     nbcidim   - number of segments on i=idim boundary
c     nbcj0     - number of segments on j=0 boundary
c     nbcjdim   - number of segments on j=jdim boundary
c     nbck0     - number of segments on k=0 boundary
c     nbckdim   - number of segments on k=kdim boundary
c     iovrlp    - grid-overlapping flag
c               = 1 the current grid is involved in grid-overlapping
c
c***** repeat above for each grid, 1 through ngrid *****
c
      read(iunit5,10)
      nbl = 0
      write(iunit11,9009)
 9009 format(6x,4hgrid,5x,5hnbci0,3x,7hnbcidim,5x,5hnbcj0,3x,7hnbcjdim,
     .                 5x,5hnbck0,3x,7hnbckdim,4x,6hiovrlp)
      msegment = 0
      lig(1) = 1
      lbg(1) = 1
      do 9007 igrid=1,ngrid
      nbl = nbl+1
      read(iunit5,*) ig,nbci0(nbl),nbcidim(nbl),nbcj0(nbl),nbcjdim(nbl),
     .          nbck0(nbl),nbckdim(nbl),iover
      if (ig.ne.igrid) then
         write(iunit11,*)' you must put these lines in order by grid!'
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      write(iunit11,36) ig,nbci0(nbl),nbcidim(nbl),nbcj0(nbl),
     .             nbcjdim(nbl),nbck0(nbl),nbckdim(nbl),iover
      if(nbci0(nbl).lt.1 .or. nbcidim(nbl).lt.1 .or.
     .   nbcj0(nbl).lt.1 .or. nbcjdim(nbl).lt.1 .or.
     .   nbck0(nbl).lt.1 .or. nbckdim(nbl).lt.1) then
        write(iunit11,'('' Error... nbci0 etc must be at least 1'')')
        call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
      iovrlp(nbl) = iover
c
c     determine size requirements for overset applications
c
      if (icall.eq.0 .and. iovrlp(nbl).eq.1) then
         if (nbl.gt.1) then
            lig(nbl) = lig(nbl-ncgg(igrid))
            lbg(nbl) = lbg(nbl-ncgg(igrid))
         end if
         call getibk0(jdimg(nbl),kdimg(nbl),idimg(nbl),nbl,itotb,
     .                itoti,maxbl,lig,lbg,ibpntsg,iipntsg,nou,bou,
     .                ibufdim,nbuf,ierrflg,myid)
      end if
c
      msegment = max(msegment,nbci0(nbl))
      msegment = max(msegment,nbcidim(nbl))
      msegment = max(msegment,nbcj0(nbl))
      msegment = max(msegment,nbcjdim(nbl))
      msegment = max(msegment,nbck0(nbl))
      msegment = max(msegment,nbckdim(nbl))
c
      ncg = ncgg(igrid)
      if (ncg.gt.0) then
         do 6879 n=1,ncg
         nbl           = nbl+1
         iovrlp(nbl)   = 0
c
c        temporarily set iovrlp(nbl) to negative to flag boundary
c        conditions for multigrid chimera
         if(iover.gt.0) iovrlp(nbl) = -1
c
         nbci0(nbl)    = nbci0(nbl-1)
         nbcidim(nbl)  = nbcidim(nbl-1)
         nbcj0(nbl)    = nbcj0(nbl-1)
         nbcjdim(nbl)  = nbcjdim(nbl-1)
         nbck0(nbl)    = nbck0(nbl-1)
         nbckdim(nbl)  = nbckdim(nbl-1)
 6879    continue
      end if
 9007 continue
c
      if (icall.eq.0) then
         rewind(21)
      end if
c
      if (msegment.gt.maxseg) then
         write(iunit11,*)' increase maxseg parameter; should be: ',
     .   msegment
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
c
c     I0 Boundary:
c     grid       - grid number (must be in order)
c     ibcinfo(1) - boundary condition type for current segment
c     ibcinfo(2) - j starting location for current segment
c     ibcinfo(3) - j ending location for current segment
c     ibcinfo(4) - k starting location for current segment
c     ibcinfo(5) - k ending location for current segment
c     ibcinfo(6) - > 0...compute force/moment on this segment
c                  < 0...do not compute force/moment on this segment
c     ibcinfo(7) - > 0...number of additional data *values* required for
c                        the boundary condition on this segment
c                  < 0...number of additional data *arrays* required for
c                        the boundary condition on this segment
c                  = 0...no additional data required for this segment
c     bcvali(m)  - data values for this segment (ibcinfo(7) > 0 only)
c     bcfilei    - pointer to name of data file containing the data 
c                  arrays for this segment (ibcinfo(7) < 0 only)
c                  (a value of 1 points to 'null')
c
c***** repeat above for 1 through nbci0 for each grid *****
c
      read(iunit5,10)
      nbl = 0
      write(iunit11,1209)
 1209 format(3hi0:,3x,4hgrid,3x,7hsegment,4x,6hbctype,6x,4hjsta,
     .             6x,4hjend,6x,4hksta,6x,4hkend,5x,5hndata)
      do 1207 igrid=1,ngrid
      nbl = nbl+1
      ifor  = iforce(nbl)
      ifo   = ifor/100
      jfo   = (ifor - ifo*100)/10
      kfo   = (ifor - ifo*100 - jfo*10)
      do 1206 iseg=1,nbci0(nbl)
      ibcinfo(nbl,iseg,1,1) = -99
      ibcinfo(nbl,iseg,6,1) = 0
      if(ifo.eq.1 .or. ifo.eq.3)ibcinfo(nbl,iseg,6,1) = 1
      bcfilei(nbl,iseg,1) = 1
      do 1206 mm=1,12
      bcvali(nbl,iseg,mm,1) = -1.e15
 1206 continue
      do 1201 iseg=1,nbci0(nbl)
      read(iunit5,*) ig,nseg,ibctyp,ibcinfo(nbl,abs(nseg),2,1),
     .                  ibcinfo(nbl,abs(nseg),3,1),
     .                  ibcinfo(nbl,abs(nseg),4,1),
     .                  ibcinfo(nbl,abs(nseg),5,1),ndata
c
      if (ibctyp .eq. 1004) then
          write(iunit11,2201)
 2201     format(44h stopping...bc1004 no longer available...use,
     .           13h 2004 instead,/,25h see bc.f for usage notes)
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if ((abs(ibctyp) .eq. 2004 .or. abs(ibctyp) .eq. 2014 .or.
     .     abs(ibctyp) .eq. 2024 .or. abs(ibctyp) .eq. 2034 .or.
     .     abs(ibctyp) .eq. 2016) .and. iviscg(nbl,1) .eq. 0) then
       write(iunit11,'('' stopping...must have ivisc(i) > 0 if'')')
       write(iunit11,'(''   bc2004/14/24/34/16 is being used on'',
     .   '' an i-face'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (abs(ibctyp) .eq. 2014 .and. ivmx .le. 3) then
          write(iunit11,'('' stopping...2014 can only be used'')')
          write(iunit11,'(''   for turb models ivisc>3'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (abs(ibctyp) .eq. 2024 .and. ivmx .ne. 30) then
          write(iunit11,'('' stopping...2024 can only be used'')')
          write(iunit11,'(''   for turb models ivisc=30'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
      if (ibctyp.ge.2000 .and. ibctyp.lt.3000 .and. ndata.eq.0) then
          write(iunit11,2202)
 2202     format(36h stopping...2000 series bc's require,
     .    15h abs(ndata) > 0)
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
      if (ibctyp.eq.9999) then
        if (iexact_trunc .eq. 0 .and. iexact_disc .eq. 0) then
          write(iunit11,2203)
 2203     format(36h stopping...9999 bc requires keyword,
     .    28h iexact_trunc or iexact_disc)
        end if
      end if
c
c     check and correct for reverse input of indicies
      itemp1 = ibcinfo(nbl,abs(nseg),2,1)
      itemp2 = ibcinfo(nbl,abs(nseg),3,1)
      if (itemp1 .gt. itemp2) then
         ibcinfo(nbl,abs(nseg),2,1) = itemp2
         ibcinfo(nbl,abs(nseg),3,1) = itemp1
      end if
      itemp1 = ibcinfo(nbl,abs(nseg),4,1)
      itemp2 = ibcinfo(nbl,abs(nseg),5,1)
      if (itemp1 .gt. itemp2) then
         ibcinfo(nbl,abs(nseg),4,1) = itemp2
         ibcinfo(nbl,abs(nseg),5,1) = itemp1
      end if
c
      ifoseg = 1
      if(nseg.lt.0) ifoseg = 0
      mfoseg = 0
      if(abs(ibctyp).eq.2004 .or. abs(ibctyp).eq.1005 .or.
     .   abs(ibctyp).eq.1006 .or. abs(ibctyp).eq.2014 .or.
     .   abs(ibctyp).eq.2024 .or. abs(ibctyp).eq.2034 .or.
     .   abs(ibctyp).eq.2016) mfoseg = 1
      ifoseg = ifoseg*mfoseg
      nseg = abs(nseg)
      if (ndata.gt.0) then
         read(iunit5,*)
         read(iunit5,*) (realval(mm),mm=1,ndata)
         do mm=1,ndata
            bcvali(nbl,nseg,mm,1) = realval(mm)
         end do
      end if
      if (ndata.lt.0) then
         nfiles = nfiles+1
         if (nfiles.gt.mxbcfil) then 
            write(iunit11,*)'too many bc files specified...increase ',
     .      'parameter mxbcfil'
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
         bcfilei(nbl,nseg,1) = nfiles
         read(iunit5,*)
         read(iunit5,'(a60)') bcfiles(nfiles)
      end if
      if (ig.ne.igrid) then
         write(iunit11,*)' you must put these lines in order by grid!'
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (ibcinfo(nbl,nseg,1,1).ne.-99) then
         write(iunit11,*) 
     .   ' stopping...attempting to set data for segment ',
     .   nseg,' more than once'
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
      ibcinfo(nbl,nseg,1,1) = ibctyp
      if (ibcinfo(nbl,nseg,2,1).eq.0 .and.
     .    ibcinfo(nbl,nseg,3,1).eq.0) then
          ibcinfo(nbl,nseg,2,1) = 1
          ibcinfo(nbl,nseg,3,1) = jdimg(nbl)
      end if
      if (ibcinfo(nbl,nseg,4,1).eq.0 .and.
     .    ibcinfo(nbl,nseg,5,1).eq.0) then
          ibcinfo(nbl,nseg,4,1) = 1
          ibcinfo(nbl,nseg,5,1) = kdimg(nbl)
      end if
      if (ifoseg.eq.0 .or. ibctyp.eq.0) then
         ibcinfo(nbl,nseg,6,1) = 0
      end if
      ibcinfo(nbl,nseg,7,1) =  ndata
c
      write(iunit11,36) ig,nseg,ibcinfo(nbl,nseg,1,1),
     .  ibcinfo(nbl,nseg,2,1),ibcinfo(nbl,nseg,3,1),
     .  ibcinfo(nbl,nseg,4,1),ibcinfo(nbl,nseg,5,1),
     .  ibcinfo(nbl,nseg,7,1)
      if (ndata.gt.0) then
         call getdhdr(datahdr,ibctyp,ndata)
         write(iunit11,1210) (datahdr(mm),mm=1,ndata)
         write(iunit11,1212) (real(bcvali(nbl,nseg,mm,1)),mm=1,ndata)
      end if
      if (ndata.lt.0) then
         write(iunit11,1213)
         nfl = bcfilei(nbl,nseg,1)
         write(iunit11,'(''     '',a60)') bcfiles(nfl)
      end if
 1210 format(10a10)
 1212 format(4f10.5,8e10.3)
 1213 format(17h     bc data file)
c
      if (ibcinfo(nbl,nseg,2,1).lt.1 .or.
     .    ibcinfo(nbl,nseg,2,1).gt.jdimg(nbl)) then
          write(iunit11,*)' stopping...jsta is out of range'
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (ibcinfo(nbl,nseg,3,1).lt.1 .or.
     .    ibcinfo(nbl,nseg,3,1).gt.jdimg(nbl)) then
          write(iunit11,*)' stopping...jend is out of range'
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (ibcinfo(nbl,nseg,4,1).lt.1 .or.
     .    ibcinfo(nbl,nseg,4,1).gt.kdimg(nbl)) then
          write(iunit11,*)' stopping...ksta is out of range'
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (ibcinfo(nbl,nseg,5,1).lt.1 .or.
     .    ibcinfo(nbl,nseg,5,1).gt.kdimg(nbl)) then
          write(iunit11,*)' stopping...kend is out of range'
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
 1201 continue
c
c     set boundary condition type array to 21 to allow
c     conservative flux accumulations
c 
      if (iemg(igrid).gt.0 .and. isg(nbl).ne.1) then
c        i=1 face of current block is embedded
         if (nbci0(nbl).gt.1) then
            write(iunit11,*)' error: embedded mesh boundary ',
     .      'at i=1 must extend over entire block face'
            write(iunit11,*)'        segmentation not allowed'
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         else
            ibcinfo(nbl,1,1,1) = 21
            ibcinfo(nbl,1,2,1) = 1
            ibcinfo(nbl,1,3,1) = jdimg(nbl)
            ibcinfo(nbl,1,4,1) = 1
            ibcinfo(nbl,1,5,1) = kdimg(nbl)
         end if
      end if
c
      ncg = ncgg(igrid)
      if (ncg.gt.0) then
         do 1279 n=1,ncg
         nbl           = nbl+1
         do 1279 nseg=1,nbci0(nbl)
         ibcinfo(nbl,nseg,1,1) = ibcinfo(nbl-1,nseg,1,1)
c
c        set extrapolation boundary condition on coarser blocks for
c        multigrid chimera. NOTE: at this stage, ALL interface bc's
c        (bctype = 0) are set to extrapolation on coarser meshes. Later
c        in subroutine bc, if any of these are 1-1 or patched interfaces,
c        the correct interface conditions WILL be applied. For 1-1/patched
c        interfaces, extrapolation will first be applied, THEN the 
c        correct interface condition will be used to rewrite the bc's. 
c        This will be reflected in the main output file under the header
c        "boundary conditions for block...". The net result is that
c        only chimera interfaces are truly set with extrapolation bc's
c
         if (iovrlp(nbl).lt.0) then
            if (ibcinfo(nbl,nseg,1,1).eq.0) ibcinfo(nbl,nseg,1,1) = 1002
         end if
c
         ibcinfo(nbl,nseg,2,1) = ibcinfo(nbl-1,nseg,2,1)/2+1
         ibcinfo(nbl,nseg,3,1) = ibcinfo(nbl-1,nseg,3,1)/2+1
         ibcinfo(nbl,nseg,4,1) = ibcinfo(nbl-1,nseg,4,1)/2+1
         ibcinfo(nbl,nseg,5,1) = ibcinfo(nbl-1,nseg,5,1)/2+1
         ibcinfo(nbl,nseg,6,1) = ibcinfo(nbl-1,nseg,6,1)
         ibcinfo(nbl,nseg,7,1) = ibcinfo(nbl-1,nseg,7,1)
         bcfilei(nbl,nseg,1)   = bcfilei(nbl-1,nseg,1)
         do 1279 l=1,12
         bcvali(nbl,nseg,l,1) = bcvali(nbl-1,nseg,l,1)
 1279    continue
      end if
c
      nblt=(igrid-1)*(ncgg(igrid)+1)+1
      isum=0
      do nseg=1,nbci0(nblt)
        isum=isum+(ibcinfo(nblt,nseg,3,1)-ibcinfo(nblt,nseg,2,1))*
     .            (ibcinfo(nblt,nseg,5,1)-ibcinfo(nblt,nseg,4,1))
      end do
      if(isum .lt. (jdimg(nblt)-1)*(kdimg(nblt)-1)) then
        write(iunit11,'('' Error.  I0 BCs do not span the space'',
     .  '' for grid '',i5)') igrid
        write(iunit11,
     .  '(''   Must exactly cover 1-'',i4,'' in j and 1-'',
     .  i4,'' in k'')') jdimg(nblt),kdimg(nblt)
        call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if(isum .gt. (jdimg(nblt)-1)*(kdimg(nblt)-1)) then
        write(iunit11,'('' Error.  I0 BCs overspan the space'',
     .  '' for grid '',i5)') igrid
        write(iunit11,
     .  '(''   Must exactly cover 1-'',i4,'' in j and 1-'',
     .  i4,'' in k'')') jdimg(nblt),kdimg(nblt)
        call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
 1207 continue
c
c
c     IDIM Boundary:
c     grid       - grid number (must be in order)
c     ibcinfo(1) - boundary condition type for current segment
c     ibcinfo(2) - j starting location for current segment
c     ibcinfo(3) - j ending location for current segment
c     ibcinfo(4) - k starting location for current segment
c     ibcinfo(5) - k ending location for current segment
c     ibcinfo(6) - > 0...compute force/moment on this segment
c                  < 0...do not compute force/moment on this segment
c     ibcinfo(7) - > 0...number of additional data *values* required for
c                        the boundary condition on this segment
c                  < 0...number of additional data *arrays* required for
c                        the boundary condition on this segment
c                  = 0...no additional data required for this segment
c     bcvali(m)  - data values for this segment (ibcinfo(7) > 0 only)
c     bcfilei    - pointer to name of data file containing the data
c                  arrays for this segment (ibcinfo(7) < 0 only)
c                  (a value of 1 points to 'null')
c
c***** repeat above for 1 through nbcidim for each grid *****
c
      read(iunit5,10)
      nbl = 0
      write(iunit11,1309)
 1309 format(5hidim:,1x,4hgrid,3x,7hsegment,4x,6hbctype,6x,4hjsta,
     .               6x,4hjend,6x,4hksta,6x,4hkend,5x,5hndata)
      do 1307 igrid=1,ngrid
      nbl = nbl+1
      ifor  = iforce(nbl)
      ifo   = ifor/100
      jfo   = (ifor - ifo*100)/10
      kfo   = (ifor - ifo*100 - jfo*10)
      do 1306 iseg=1,nbcidim(nbl)
      ibcinfo(nbl,iseg,1,2) = -99
      ibcinfo(nbl,iseg,6,2) = 0
      if(ifo.eq.2 .or. ifo.eq.3)ibcinfo(nbl,iseg,6,2) = 1
      bcfilei(nbl,iseg,2) = 1
      do 1306 mm=1,12
      bcvali(nbl,iseg,mm,2) = -1.e15
 1306 continue
      do 1301 iseg=1,nbcidim(nbl)
      read(iunit5,*) ig,nseg,ibctyp,ibcinfo(nbl,abs(nseg),2,2),
     .                  ibcinfo(nbl,abs(nseg),3,2),
     .                  ibcinfo(nbl,abs(nseg),4,2),
     .                  ibcinfo(nbl,abs(nseg),5,2),ndata
c
      if (ibctyp .eq. 1004) then
          write(iunit11,2201)
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if ((abs(ibctyp) .eq. 2004 .or. abs(ibctyp) .eq. 2014 .or.
     .     abs(ibctyp) .eq. 2024 .or. abs(ibctyp) .eq. 2034 .or.
     .     abs(ibctyp) .eq. 2016) .and. iviscg(nbl,1) .eq. 0) then
       write(iunit11,'('' stopping...must have ivisc(i) > 0 if'')')
       write(iunit11,'(''   bc2004/14/24/34/16 is being used on'',
     .   '' an i-face'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (abs(ibctyp) .eq. 2014 .and. ivmx .le. 3) then
          write(iunit11,'('' stopping...2014 can only be used'')')
          write(iunit11,'(''   for turb models ivisc>3'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (abs(ibctyp) .eq. 2024 .and. ivmx .ne. 30) then
          write(iunit11,'('' stopping...2024 can only be used'')')
          write(iunit11,'(''   for turb models ivisc=30'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
      if (ibctyp.ge.2000 .and. ibctyp.lt.3000 .and. ndata.eq.0) then
          write(iunit11,2202)
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
      if (ibctyp.eq.9999) then
        if (iexact_trunc .eq. 0 .and. iexact_disc .eq. 0) then
          write(iunit11,2203)
        end if
      end if
c
c     check and correct for reverse input of indicies
      itemp1 = ibcinfo(nbl,abs(nseg),2,2)
      itemp2 = ibcinfo(nbl,abs(nseg),3,2)
      if (itemp1 .gt. itemp2) then
         ibcinfo(nbl,abs(nseg),2,2) = itemp2
         ibcinfo(nbl,abs(nseg),3,2) = itemp1
      end if
      itemp1 = ibcinfo(nbl,abs(nseg),4,2)
      itemp2 = ibcinfo(nbl,abs(nseg),5,2)
      if (itemp1 .gt. itemp2) then
         ibcinfo(nbl,abs(nseg),4,2) = itemp2
         ibcinfo(nbl,abs(nseg),5,2) = itemp1
      end if
c
      ifoseg = 1
      if(nseg.lt.0) ifoseg = 0
      mfoseg = 0
      if(abs(ibctyp).eq.2004 .or. abs(ibctyp).eq.1005 .or.
     .   abs(ibctyp).eq.1006 .or. abs(ibctyp).eq.2014 .or.
     .   abs(ibctyp).eq.2024 .or. abs(ibctyp).eq.2034 .or.
     .   abs(ibctyp).eq.2016) mfoseg = 1
      ifoseg = ifoseg*mfoseg
      nseg = abs(nseg)
      if (ndata.gt.0) then
         read(iunit5,*)
         read(iunit5,*) (realval(mm),mm=1,ndata)
         do mm=1,ndata
            bcvali(nbl,nseg,mm,2) = realval(mm)
         end do
      end if
      if (ndata.lt.0) then
         nfiles = nfiles+1
         if (nfiles.gt.mxbcfil) then
            write(iunit11,*)'too many bc files specified...increase ',
     .      'parameter mxbcfil'
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
         bcfilei(nbl,nseg,2) = nfiles
         read(iunit5,*)
         read(iunit5,'(a60)') bcfiles(nfiles)
      end if
      if (ig.ne.igrid) then
         write(iunit11,*)' you must put these lines in order by grid!'
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (ibcinfo(nbl,nseg,1,2).ne.-99) then
         write(iunit11,*) 
     .   ' stopping...attempting to set data for segment ',
     .   nseg,' more than once'
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
      ibcinfo(nbl,nseg,1,2) = ibctyp
      if (ibcinfo(nbl,nseg,2,2).eq.0 .and.
     .    ibcinfo(nbl,nseg,3,2).eq.0) then
          ibcinfo(nbl,nseg,2,2) = 1
          ibcinfo(nbl,nseg,3,2) = jdimg(nbl)
      end if
      if (ibcinfo(nbl,nseg,4,2).eq.0 .and.
     .    ibcinfo(nbl,nseg,5,2).eq.0) then
          ibcinfo(nbl,nseg,4,2) = 1
          ibcinfo(nbl,nseg,5,2) = kdimg(nbl)
      end if
      if (ifoseg.eq.0 .or. ibctyp.eq.0) then
         ibcinfo(nbl,nseg,6,2) = 0
      end if
      ibcinfo(nbl,nseg,7,2) = ndata
c
      write(iunit11,36) ig,nseg,ibcinfo(nbl,nseg,1,2),
     .   ibcinfo(nbl,nseg,2,2),ibcinfo(nbl,nseg,3,2),
     .   ibcinfo(nbl,nseg,4,2),ibcinfo(nbl,nseg,5,2),
     .   ibcinfo(nbl,nseg,7,2)
      if (ndata.gt.0) then
         call getdhdr(datahdr,ibctyp,ndata)
         write(iunit11,1210) (datahdr(mm),mm=1,ndata)
         write(iunit11,1212) (real(bcvali(nbl,nseg,mm,2)),mm=1,ndata)
      end if
      if (ndata.lt.0) then
         nfl = bcfilei(nbl,nseg,2)
         write(iunit11,1213)
         write(iunit11,'(''     '',a60)') bcfiles(nfl)
      end if
c
      if (ibcinfo(nbl,nseg,2,2).lt.1 .or.
     .    ibcinfo(nbl,nseg,2,2).gt.jdimg(nbl)) then
          write(iunit11,*)' stopping...jsta is out of range'
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (ibcinfo(nbl,nseg,3,2).lt.1 .or.
     .    ibcinfo(nbl,nseg,3,2).gt.jdimg(nbl)) then
          write(iunit11,*)' stopping...jend is out of range'
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (ibcinfo(nbl,nseg,4,2).lt.1 .or.
     .    ibcinfo(nbl,nseg,4,2).gt.kdimg(nbl)) then
          write(iunit11,*)' stopping...ksta is out of range'
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (ibcinfo(nbl,nseg,5,2).lt.1 .or.
     .    ibcinfo(nbl,nseg,5,2).gt.kdimg(nbl)) then
          write(iunit11,*)' stopping...kend is out of range'
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
 1301 continue
c
c     set boundary condition type array to 21 to allow
c     conservative flux accumulations
c 
      if (iemg(igrid).gt.0 .and. ieg(nbl).ne.idimg(nblcg(nbl))) then
c        i=idim face of current block is embedded
         if (nbcidim(nbl).gt.1) then
            write(iunit11,*)' error: embedded mesh boundary ',
     .      'at i=idim must extend over entire block face'
            write(iunit11,*)'        segmentation not allowed'
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         else
            ibcinfo(nbl,1,1,2) = 21
            ibcinfo(nbl,1,2,2) = 1
            ibcinfo(nbl,1,3,2) = jdimg(nbl)
            ibcinfo(nbl,1,4,2) = 1
            ibcinfo(nbl,1,5,2) = kdimg(nbl)
         end if
      end if
c
      ncg = ncgg(igrid)
      if (ncg.gt.0) then
         do 1379 n=1,ncg
         nbl           = nbl+1
         do 1379 nseg=1,nbcidim(nbl)
         ibcinfo(nbl,nseg,1,2) = ibcinfo(nbl-1,nseg,1,2)
c
c        set extrapolation boundary condition on coarser blocks for
c        multigrid chimera. SEE NOTE IN THE I=I0 BOUNDARY SECTION
         if (iovrlp(nbl).lt.0) then
            if (ibcinfo(nbl,nseg,1,2).eq.0) ibcinfo(nbl,nseg,1,2) = 1002
         end if
c
         ibcinfo(nbl,nseg,2,2) = ibcinfo(nbl-1,nseg,2,2)/2+1
         ibcinfo(nbl,nseg,3,2) = ibcinfo(nbl-1,nseg,3,2)/2+1
         ibcinfo(nbl,nseg,4,2) = ibcinfo(nbl-1,nseg,4,2)/2+1
         ibcinfo(nbl,nseg,5,2) = ibcinfo(nbl-1,nseg,5,2)/2+1
         ibcinfo(nbl,nseg,6,2) = ibcinfo(nbl-1,nseg,6,2)
         ibcinfo(nbl,nseg,7,2) = ibcinfo(nbl-1,nseg,7,2)
         bcfilei(nbl,nseg,2)   = bcfilei(nbl-1,nseg,2)
         do 1379 l=1,12
         bcvali(nbl,nseg,l,2) = bcvali(nbl-1,nseg,l,2)
 1379    continue
      end if
c
      nblt=(igrid-1)*(ncgg(igrid)+1)+1
      isum=0
      do nseg=1,nbcidim(nblt)
        isum=isum+(ibcinfo(nblt,nseg,3,2)-ibcinfo(nblt,nseg,2,2))*
     .            (ibcinfo(nblt,nseg,5,2)-ibcinfo(nblt,nseg,4,2))
      enddo
      if(isum .lt. (jdimg(nblt)-1)*(kdimg(nblt)-1)) then
        write(iunit11,'('' Error.  IDIM BCs do not span the space'',
     .  '' for grid '',i5)') igrid
        write(iunit11,
     .  '(''   Must exactly cover 1-'',i4,'' in j and 1-'',
     .  i4,'' in k'')') jdimg(nblt),kdimg(nblt)
        call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if(isum .gt. (jdimg(nblt)-1)*(kdimg(nblt)-1)) then
        write(iunit11,'('' Error.  IDIM BCs overspan the space'',
     .  '' for grid '',i5)') igrid
        write(iunit11,
     .  '(''   Must exactly cover 1-'',i4,'' in j and 1-'',
     .  i4,'' in k'')') jdimg(nblt),kdimg(nblt)
        call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
 1307 continue
c
c
c     J0 Boundary:
c     grid       - grid number (must be in order)
c     jbcinfo(1) - boundary condition type for current segment
c     jbcinfo(2) - i starting location for current segment
c     jbcinfo(3) - i ending location for current segment
c     jbcinfo(4) - k starting location for current segment
c     jbcinfo(5) - k ending location for current segment
c     jbcinfo(6) - > 0...compute force/moment on this segment
c                  < 0...do not compute force/moment on this segment
c     jbcinfo(7) - > 0...number of additional data *values* required for
c                        the boundary condition on this segment
c                  < 0...number of additional data *arrays* required for
c                        the boundary condition on this segment
c                  = 0...no additional data required for this segment
c     bcvalj(m)  - data values for this segment (jbcinfo(7) > 0 only)
c     bcfilej    - pointer to name of data file containing the data
c                  arrays for this segment (jbcinfo(7) < 0 only)
c                  (a value of 1 points to 'null')
c
c***** repeat above for 1 through nbcj0 for each grid *****
c
      read(iunit5,10)
      nbl = 0
      write(iunit11,1409)
 1409 format(3hj0:,3x,4hgrid,3x,7hsegment,4x,6hbctype,6x,4hista,
     .             6x,4hiend,6x,4hksta,6x,4hkend,5x,5hndata)
      do 1407 igrid=1,ngrid
      nbl = nbl+1
      ifor  = iforce(nbl)
      ifo   = ifor/100
      jfo   = (ifor - ifo*100)/10
      kfo   = (ifor - ifo*100 - jfo*10)
      do 1406 iseg=1,nbcj0(nbl)
      jbcinfo(nbl,iseg,1,1) = -99
      jbcinfo(nbl,iseg,6,1) = 0
      if(jfo.eq.1 .or. jfo.eq.3)jbcinfo(nbl,iseg,6,1) = 1
      bcfilej(nbl,iseg,1) = 1
      do 1406 mm=1,12
      bcvalj(nbl,iseg,mm,1) = -1.e15
 1406 continue
      do 1401 iseg=1,nbcj0(nbl)
      read(iunit5,*) ig,nseg,jbctyp,jbcinfo(nbl,abs(nseg),2,1),
     .                  jbcinfo(nbl,abs(nseg),3,1),
     .                  jbcinfo(nbl,abs(nseg),4,1),
     .                  jbcinfo(nbl,abs(nseg),5,1),ndata
c
      if (jbctyp .eq. 1004) then
          write(iunit11,2201)
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if ((abs(jbctyp) .eq. 2004 .or. abs(jbctyp) .eq. 2014 .or.
     .     abs(jbctyp) .eq. 2024 .or. abs(jbctyp) .eq. 2034 .or.
     .     abs(jbctyp) .eq. 2016) .and. iviscg(nbl,2) .eq. 0) then
       write(iunit11,'('' stopping...must have ivisc(j) > 0 if'')')
       write(iunit11,'(''   bc2004/14/24/34/16 is being used on'',
     .   '' a j-face'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (abs(jbctyp) .eq. 2014 .and. ivmx .le. 3) then
          write(iunit11,'('' stopping...2014 can only be used'')')
          write(iunit11,'(''   for turb models ivisc>3'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (abs(jbctyp) .eq. 2024 .and. ivmx .ne. 30) then
          write(iunit11,'('' stopping...2024 can only be used'')')
          write(iunit11,'(''   for turb models ivisc=30'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
      if (jbctyp.ge.2000 .and. jbctyp.lt.3000 .and. ndata.eq.0) then
          write(iunit11,2202)
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
      if (jbctyp.eq.9999) then
        if (iexact_trunc .eq. 0 .and. iexact_disc .eq. 0) then
          write(iunit11,2203)
        end if
      end if
c
c     check and correct for reverse input of indicies
      jtemp1 = jbcinfo(nbl,abs(nseg),2,1)
      jtemp2 = jbcinfo(nbl,abs(nseg),3,1)
      if (jtemp1 .gt. jtemp2) then
         jbcinfo(nbl,abs(nseg),2,1) = jtemp2
         jbcinfo(nbl,abs(nseg),3,1) = jtemp1
      end if
      jtemp1 = jbcinfo(nbl,abs(nseg),4,1)
      jtemp2 = jbcinfo(nbl,abs(nseg),5,1)
      if (jtemp1 .gt. jtemp2) then
         jbcinfo(nbl,abs(nseg),4,1) = jtemp2
         jbcinfo(nbl,abs(nseg),5,1) = jtemp1
      end if
c
      jfoseg = 1
      if(nseg.lt.0) jfoseg = 0
      mfoseg = 0
      if(abs(jbctyp).eq.2004 .or. abs(jbctyp).eq.1005 .or.
     .   abs(jbctyp).eq.1006 .or. abs(jbctyp).eq.2014 .or.
     .   abs(jbctyp).eq.2024 .or. abs(jbctyp).eq.2034 .or.
     .   abs(jbctyp).eq.2016) mfoseg = 1
      jfoseg = jfoseg*mfoseg
      nseg = abs(nseg)
      if (ndata.gt.0) then
         read(iunit5,*)
         read(iunit5,*) (realval(mm),mm=1,ndata)
         do mm=1,ndata
            bcvalj(nbl,nseg,mm,1) = realval(mm)
         end do
      end if
      if (ndata.lt.0) then
         nfiles = nfiles+1
         if (nfiles.gt.mxbcfil) then
            write(iunit11,*)'too many bc files specified...increase ',
     .      'parameter mxbcfil'
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
         bcfilej(nbl,nseg,1) = nfiles
         read(iunit5,*)
         read(iunit5,'(a60)') bcfiles(nfiles)
      end if
      if (ig.ne.igrid) then
         write(iunit11,*)' you must put these lines in order by grid!'
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (jbcinfo(nbl,nseg,1,1).ne.-99) then
         write(iunit11,*) 
     .   ' stopping...attempting to set data for segment ',
     .   nseg,' more than once'
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
      jbcinfo(nbl,nseg,1,1) = jbctyp
      if (jbcinfo(nbl,nseg,2,1).eq.0 .and.
     .    jbcinfo(nbl,nseg,3,1).eq.0) then
          jbcinfo(nbl,nseg,2,1) = 1
          jbcinfo(nbl,nseg,3,1) = idimg(nbl)
      end if
      if (jbcinfo(nbl,nseg,4,1).eq.0 .and.
     .    jbcinfo(nbl,nseg,5,1).eq.0) then
          jbcinfo(nbl,nseg,4,1) = 1
          jbcinfo(nbl,nseg,5,1) = kdimg(nbl)
      end if
      if (jfoseg.eq.0 .or. jbctyp.eq.0) then
         jbcinfo(nbl,nseg,6,1) = 0
      end if
      jbcinfo(nbl,nseg,7,1) = ndata
c
      write(iunit11,36) ig,nseg,jbcinfo(nbl,nseg,1,1),
     . jbcinfo(nbl,nseg,2,1),jbcinfo(nbl,nseg,3,1),
     . jbcinfo(nbl,nseg,4,1),jbcinfo(nbl,nseg,5,1),
     . jbcinfo(nbl,nseg,7,1)
      if (ndata.gt.0) then
         call getdhdr(datahdr,jbctyp,ndata)
         write(iunit11,1210) (datahdr(mm),mm=1,ndata)
         write(iunit11,1212) (real(bcvalj(nbl,nseg,mm,1)),mm=1,ndata)
      end if
      if (ndata.lt.0) then
         nfl = bcfilej(nbl,nseg,1)
         write(iunit11,1213)
         write(iunit11,'(''     '',a60)') bcfiles(nfl)
      end if
c
      if (jbcinfo(nbl,nseg,2,1).lt.1 .or.
     .    jbcinfo(nbl,nseg,2,1).gt.idimg(nbl)) then
          write(iunit11,*)' stopping...ista is out of range'
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (jbcinfo(nbl,nseg,3,1).lt.1 .or.
     .    jbcinfo(nbl,nseg,3,1).gt.idimg(nbl)) then
          write(iunit11,*)' stopping...iend is out of range'
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (jbcinfo(nbl,nseg,4,1).lt.1 .or.
     .    jbcinfo(nbl,nseg,4,1).gt.kdimg(nbl)) then
          write(iunit11,*)' stopping...ksta is out of range'
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (jbcinfo(nbl,nseg,5,1).lt.1 .or.
     .    jbcinfo(nbl,nseg,5,1).gt.kdimg(nbl)) then
          write(iunit11,*)' stopping...kend is out of range'
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
 1401 continue
c
c     set boundary condition type array to 21 to allow
c     conservative flux accumulations
c
      if (iemg(igrid).gt.0 .and. jsg(nbl).ne.1) then
c        j=1 face of current block is embedded
         if (nbcj0(nbl).gt.1) then
            write(iunit11,*)' error: embedded mesh boundary ',
     .      'at j=1 must extend over entire block face'
            write(iunit11,*)'        segmentation not allowed'
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         else
            jbcinfo(nbl,1,1,1) = 21
            jbcinfo(nbl,1,2,1) = 1
            jbcinfo(nbl,1,3,1) = idimg(nbl)
            jbcinfo(nbl,1,4,1) = 1
            jbcinfo(nbl,1,5,1) = kdimg(nbl)
         end if
      end if
c
      ncg = ncgg(igrid)
      if (ncg.gt.0) then
         do 1479 n=1,ncg
         nbl           = nbl+1
         do 1479 nseg=1,nbcj0(nbl)
         jbcinfo(nbl,nseg,1,1) = jbcinfo(nbl-1,nseg,1,1)
c
c        set extrapolation boundary condition on coarser blocks for
c        multigrid chimera. SEE NOTE IN THE I=I0 BOUNDARY SECTION
         if (iovrlp(nbl).lt.0) then
            if (jbcinfo(nbl,nseg,1,1).eq.0) jbcinfo(nbl,nseg,1,1) = 1002
         end if
c
         jbcinfo(nbl,nseg,2,1) = jbcinfo(nbl-1,nseg,2,1)/2+1
         jbcinfo(nbl,nseg,3,1) = jbcinfo(nbl-1,nseg,3,1)/2+1
         jbcinfo(nbl,nseg,4,1) = jbcinfo(nbl-1,nseg,4,1)/2+1
         jbcinfo(nbl,nseg,5,1) = jbcinfo(nbl-1,nseg,5,1)/2+1
         jbcinfo(nbl,nseg,6,1) = jbcinfo(nbl-1,nseg,6,1)
         jbcinfo(nbl,nseg,7,1) = jbcinfo(nbl-1,nseg,7,1)
         bcfilej(nbl,nseg,1)   = bcfilej(nbl-1,nseg,1)
         do 1479 l=1,12
         bcvalj(nbl,nseg,l,1) = bcvalj(nbl-1,nseg,l,1)
 1479    continue
      end if
c
      nblt=(igrid-1)*(ncgg(igrid)+1)+1
      isum=0
      do nseg=1,nbcj0(nblt)
        isum=isum+(jbcinfo(nblt,nseg,3,1)-jbcinfo(nblt,nseg,2,1))*
     .            (jbcinfo(nblt,nseg,5,1)-jbcinfo(nblt,nseg,4,1))
      enddo
      if(isum .lt. (idimg(nblt)-1)*(kdimg(nblt)-1)) then
        write(iunit11,'('' Error.  J0 BCs do not span the space'',
     .  '' for grid '',i5)') igrid
        write(iunit11,
     .  '(''   Must exactly cover 1-'',i4,'' in i and 1-'',
     .  i4,'' in k'')') idimg(nblt),kdimg(nblt)
        call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if(isum .gt. (idimg(nblt)-1)*(kdimg(nblt)-1)) then
        write(iunit11,'('' Error.  J0 BCs overspan the space'',
     .  '' for grid '',i5)') igrid
        write(iunit11,
     .  '(''   Must exactly cover 1-'',i4,'' in i and 1-'',
     .  i4,'' in k'')') idimg(nblt),kdimg(nblt)
        call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
 1407 continue
c
c     JDIM Boundary:
c     grid       - grid number (must be in order)
c     jbcinfo(1) - boundary condition type for current segment
c     jbcinfo(2) - i starting location for current segment
c     jbcinfo(3) - i ending location for current segment
c     jbcinfo(4) - k starting location for current segment
c     jbcinfo(5) - k ending location for current segment
c     jbcinfo(6) - > 0...compute force/moment on this segment
c                  < 0...do not compute force/moment on this segment
c     jbcinfo(7) - > 0...number of additional data *values* required for
c                        the boundary condition on this segment
c                  < 0...number of additional data *arrays* required for
c                        the boundary condition on this segment
c                  = 0...no additional data required for this segment
c     bcvalj(m)  - data values for this segment (jbcinfo(7) > 0 only)
c     bcfilej    - pointer to name of data file containing the data
c                  arrays for this segment (jbcinfo(7) < 0 only)
c                  (a value of 1 points to 'null')
c
c***** repeat above for 1 through nbcjdim for each grid *****
c
      read(iunit5,10)
      nbl = 0
      write(iunit11,1509)
 1509 format(5hjdim:,1x,4hgrid,3x,7hsegment,4x,6hbctype,6x,4hista,
     .               6x,4hiend,6x,4hksta,6x,4hkend,5x,5hndata)
      do 1507 igrid=1,ngrid
      nbl = nbl+1
      ifor  = iforce(nbl)
      ifo   = ifor/100
      jfo   = (ifor - ifo*100)/10
      kfo   = (ifor - ifo*100 - jfo*10)
      do 1506 iseg=1,nbcjdim(nbl)
      jbcinfo(nbl,iseg,1,2) = -99
      jbcinfo(nbl,iseg,6,2) = 0
      if(jfo.eq.2 .or. jfo.eq.3)jbcinfo(nbl,iseg,6,2) = 1
      bcfilej(nbl,iseg,2) = 1
      do 1506 mm=1,12
      bcvalj(nbl,iseg,mm,2) = -1.e15
 1506 continue
      do 1501 iseg=1,nbcjdim(nbl)
      read(iunit5,*) ig,nseg,jbctyp,jbcinfo(nbl,abs(nseg),2,2),
     .                  jbcinfo(nbl,abs(nseg),3,2),
     .                  jbcinfo(nbl,abs(nseg),4,2),
     .                  jbcinfo(nbl,abs(nseg),5,2),ndata
c
      if (jbctyp .eq. 1004) then
          write(iunit11,2201)
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if ((abs(jbctyp) .eq. 2004 .or. abs(jbctyp) .eq. 2014 .or.
     .     abs(jbctyp) .eq. 2024 .or. abs(jbctyp) .eq. 2034 .or.
     .     abs(jbctyp) .eq. 2016) .and. iviscg(nbl,2) .eq. 0) then
       write(iunit11,'('' stopping...must have ivisc(j) > 0 if'')')
       write(iunit11,'(''   bc2004/14/24/34/16 is being used on'',
     .   '' a j-face'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (abs(jbctyp) .eq. 2014 .and. ivmx .le. 3) then
          write(iunit11,'('' stopping...2014 can only be used'')')
          write(iunit11,'(''   for turb models ivisc>3'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (abs(jbctyp) .eq. 2024 .and. ivmx .ne. 30) then
          write(iunit11,'('' stopping...2024 can only be used'')')
          write(iunit11,'(''   for turb models ivisc=30'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
      if (jbctyp.ge.2000 .and. jbctyp.lt.3000 .and. ndata.eq.0) then
          write(iunit11,2202)
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
      if (jbctyp.eq.9999) then
        if (iexact_trunc .eq. 0 .and. iexact_disc .eq. 0) then
          write(iunit11,2203)
        end if
      end if
c
c     check and correct for reverse input of indicies
      jtemp1 = jbcinfo(nbl,abs(nseg),2,2)
      jtemp2 = jbcinfo(nbl,abs(nseg),3,2)
      if (jtemp1 .gt. jtemp2) then
         jbcinfo(nbl,abs(nseg),2,2) = jtemp2
         jbcinfo(nbl,abs(nseg),3,2) = jtemp1
      end if
      jtemp1 = jbcinfo(nbl,abs(nseg),4,2)
      jtemp2 = jbcinfo(nbl,abs(nseg),5,2)
      if (jtemp1 .gt. jtemp2) then
         jbcinfo(nbl,abs(nseg),4,2) = jtemp2
         jbcinfo(nbl,abs(nseg),5,2) = jtemp1
      end if
c
      jfoseg = 1
      if(nseg.lt.0) jfoseg = 0
      mfoseg = 0
      if(abs(jbctyp).eq.2004 .or. abs(jbctyp).eq.1005 .or.
     .   abs(jbctyp).eq.1006 .or. abs(jbctyp).eq.2014 .or.
     .   abs(jbctyp).eq.2024 .or. abs(jbctyp).eq.2034 .or.
     .   abs(jbctyp).eq.2016) mfoseg = 1
      jfoseg = jfoseg*mfoseg
      nseg = abs(nseg)
      if (ndata.gt.0) then
         read(iunit5,*)
         read(iunit5,*) (realval(mm),mm=1,ndata)
         do mm=1,ndata
            bcvalj(nbl,nseg,mm,2) = realval(mm)
         end do
      end if
      if (ndata.lt.0) then
         nfiles = nfiles+1
         if (nfiles.gt.mxbcfil) then
            write(iunit11,*)'too many bc files specified...increase ',
     .      'parameter mxbcfil'
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
         bcfilej(nbl,nseg,2) = nfiles
         read(iunit5,*)
         read(iunit5,'(a60)') bcfiles(nfiles)
      end if
      if (ig.ne.igrid) then
         write(iunit11,*)' you must put these lines in order by grid!'
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (jbcinfo(nbl,nseg,1,2).ne.-99) then
         write(iunit11,*) 
     .   ' stopping...attempting to set data for segment ',
     .   nseg,' more than once'
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
      jbcinfo(nbl,nseg,1,2) = jbctyp
      if (jbcinfo(nbl,nseg,2,2).eq.0 .and.
     .    jbcinfo(nbl,nseg,3,2).eq.0) then
          jbcinfo(nbl,nseg,2,2) = 1
          jbcinfo(nbl,nseg,3,2) = idimg(nbl)
      end if
      if (jbcinfo(nbl,nseg,4,2).eq.0 .and.
     .    jbcinfo(nbl,nseg,5,2).eq.0) then
          jbcinfo(nbl,nseg,4,2) = 1
          jbcinfo(nbl,nseg,5,2) = kdimg(nbl)
      end if
      if (jfoseg.eq.0 .or. jbctyp.eq.0) then
         jbcinfo(nbl,nseg,6,2) = 0
      end if
      jbcinfo(nbl,nseg,7,2) = ndata
c
      write(iunit11,36) ig,nseg,jbcinfo(nbl,nseg,1,2),
     . jbcinfo(nbl,nseg,2,2),jbcinfo(nbl,nseg,3,2),
     . jbcinfo(nbl,nseg,4,2),jbcinfo(nbl,nseg,5,2),
     . jbcinfo(nbl,nseg,7,2)
      if (ndata.gt.0) then
         call getdhdr(datahdr,jbctyp,ndata)
         write(iunit11,1210) (datahdr(mm),mm=1,ndata)
         write(iunit11,1212) (real(bcvalj(nbl,nseg,mm,2)),mm=1,ndata)
      end if
      if (ndata.lt.0) then
         nfl = bcfilej(nbl,nseg,2)
         write(iunit11,1213)
         write(iunit11,'(''     '',a60)') bcfiles(nfl)
      end if
c
      if (jbcinfo(nbl,nseg,2,2).lt.1 .or.
     .    jbcinfo(nbl,nseg,2,2).gt.idimg(nbl)) then
          write(iunit11,*)' stopping...ista is out of range'
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (jbcinfo(nbl,nseg,3,2).lt.1 .or.
     .    jbcinfo(nbl,nseg,3,2).gt.idimg(nbl)) then
          write(iunit11,*)' stopping...iend is out of range'
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (jbcinfo(nbl,nseg,4,2).lt.1 .or.
     .    jbcinfo(nbl,nseg,4,2).gt.kdimg(nbl)) then
          write(iunit11,*)' stopping...ksta is out of range'
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (jbcinfo(nbl,nseg,5,2).lt.1 .or.
     .    jbcinfo(nbl,nseg,5,2).gt.kdimg(nbl)) then
          write(iunit11,*)' stopping...kend is out of range'
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
 1501 continue
c
c     set boundary condition type array to 21 to allow
c     conservative flux accumulations
c
      if (iemg(igrid).gt.0 .and. jeg(nbl).ne.jdimg(nblcg(nbl))) then
c        j=jdim face of current block is embedded
         if (nbcjdim(nbl).gt.1) then
            write(iunit11,*)' error: embedded mesh boundary ',
     .      'at j=jdim must extend over entire block face'
            write(iunit11,*)'        segmentation not allowed'
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         else
            jbcinfo(nbl,1,1,2) = 21
            jbcinfo(nbl,1,2,2) = 1
            jbcinfo(nbl,1,3,2) = idimg(nbl)
            jbcinfo(nbl,1,4,2) = 1
            jbcinfo(nbl,1,5,2) = kdimg(nbl)
         end if
      end if
c
      ncg = ncgg(igrid)
      if (ncg.gt.0) then
         do 1579 n=1,ncg
         nbl           = nbl+1
         do 1579 nseg=1,nbcjdim(nbl)
         jbcinfo(nbl,nseg,1,2) = jbcinfo(nbl-1,nseg,1,2)
c
c        set extrapolation boundary condition on coarser blocks for
c        multigrid chimera. SEE NOTE IN THE I=I0 BOUNDARY SECTION
         if (iovrlp(nbl).lt.0) then
            if (jbcinfo(nbl,nseg,1,2).eq.0) jbcinfo(nbl,nseg,1,2) = 1002
         end if
c
         jbcinfo(nbl,nseg,2,2) = jbcinfo(nbl-1,nseg,2,2)/2+1
         jbcinfo(nbl,nseg,3,2) = jbcinfo(nbl-1,nseg,3,2)/2+1
         jbcinfo(nbl,nseg,4,2) = jbcinfo(nbl-1,nseg,4,2)/2+1
         jbcinfo(nbl,nseg,5,2) = jbcinfo(nbl-1,nseg,5,2)/2+1
         jbcinfo(nbl,nseg,6,2) = jbcinfo(nbl-1,nseg,6,2)
         jbcinfo(nbl,nseg,7,2) = jbcinfo(nbl-1,nseg,7,2)
         bcfilej(nbl,nseg,2)   = bcfilej(nbl-1,nseg,2)
         do 1579 l=1,12
         bcvalj(nbl,nseg,l,2) = bcvalj(nbl-1,nseg,l,2)
 1579    continue
      end if
c
      nblt=(igrid-1)*(ncgg(igrid)+1)+1
      isum=0
      do nseg=1,nbcjdim(nblt)
        isum=isum+(jbcinfo(nblt,nseg,3,2)-jbcinfo(nblt,nseg,2,2))*
     .            (jbcinfo(nblt,nseg,5,2)-jbcinfo(nblt,nseg,4,2))
      enddo
      if(isum .lt. (idimg(nblt)-1)*(kdimg(nblt)-1)) then
        write(iunit11,'('' Error.  JDIM BCs do not span the space'',
     .  '' for grid '',i5)') igrid
        write(iunit11,
     .  '(''   Must exactly cover 1-'',i4,'' in i and 1-'',
     .  i4,'' in k'')') idimg(nblt),kdimg(nblt)
        call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if(isum .gt. (idimg(nblt)-1)*(kdimg(nblt)-1)) then
        write(iunit11,'('' Error.  JDIM BCs overspan the space'',
     .  '' for grid '',i5)') igrid
        write(iunit11,
     .  '(''   Must exactly cover 1-'',i4,'' in i and 1-'',
     .  i4,'' in k'')') idimg(nblt),kdimg(nblt)
        call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
 1507 continue
c
c     K0 Boundary:
c     grid       - grid number (must be in order)
c     kbcinfo(1) - boundary condition type for current segment
c     kbcinfo(2) - i starting location for current segment
c     kbcinfo(3) - i ending location for current segment
c     kbcinfo(4) - j starting location for current segment
c     kbcinfo(5) - j ending location for current segment
c     kbcinfo(6) - > 0...compute force/moment on this segment
c                  < 0...do not compute force/moment on this segment
c     kbcinfo(7) - > 0...number of additional data *values* required for
c                        the boundary condition on this segment
c                  < 0...number of additional data *arrays* required for
c                        the boundary condition on this segment
c                  = 0...no additional data required for this segment
c     bcvalk(m)  - data values for this segment (kbcinfo(7) > 0 only)
c     bcfilek    - pointer to name of data file containing the data
c                  arrays for this segment (kbcinfo(7) < 0 only)
c                  (a value of 1 points to 'null')
c
c***** repeat above for 1 through nbck0 for each grid *****
c
      read(iunit5,10)
      nbl = 0
      write(iunit11,1609)
 1609 format(3hk0:,3x,4hgrid,3x,7hsegment,4x,6hbctype,6x,4hista,
     .             6x,4hiend,6x,4hjsta,6x,4hjend,5x,5hndata)
      do 1607 igrid=1,ngrid
      nbl = nbl+1
      ifor  = iforce(nbl)
      ifo   = ifor/100
      jfo   = (ifor - ifo*100)/10
      kfo   = (ifor - ifo*100 - jfo*10)
      do 1606 iseg=1,nbck0(nbl)
      kbcinfo(nbl,iseg,1,1) = -99
      kbcinfo(nbl,iseg,6,1) = 0
      if(kfo.eq.1 .or. kfo.eq.3)kbcinfo(nbl,iseg,6,1) = 1
      bcfilek(nbl,iseg,1) = 1
      do 1606 mm=1,12
      bcvalk(nbl,iseg,mm,1) = -1.e15
 1606 continue
      do 1601 iseg=1,nbck0(nbl)
      read(iunit5,*) ig,nseg,kbctyp,kbcinfo(nbl,abs(nseg),2,1),
     .                  kbcinfo(nbl,abs(nseg),3,1),
     .                  kbcinfo(nbl,abs(nseg),4,1),
     .                  kbcinfo(nbl,abs(nseg),5,1),ndata
      if (kbctyp .eq. 1004) then
          write(iunit11,2201)
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if ((abs(kbctyp) .eq. 2004 .or. abs(kbctyp) .eq. 2014 .or.
     .     abs(kbctyp) .eq. 2024 .or. abs(kbctyp) .eq. 2034 .or.
     .     abs(kbctyp) .eq. 2016) .and. iviscg(nbl,3) .eq. 0) then
       write(iunit11,'('' stopping...must have ivisc(k) > 0 if'')')
       write(iunit11,'(''   bc2004/14/24/34/16 is being used on'',
     .   '' a k-face'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (abs(kbctyp) .eq. 2014 .and. ivmx .le. 3) then
          write(iunit11,'('' stopping...2014 can only be used'')')
          write(iunit11,'(''   for turb models ivisc>3'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (abs(kbctyp) .eq. 2024 .and. ivmx .ne. 30) then
          write(iunit11,'('' stopping...2024 can only be used'')')
          write(iunit11,'(''   for turb models ivisc=30'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
      if (kbctyp.ge.2000 .and. kbctyp.lt.3000 .and. ndata.eq.0) then
          write(iunit11,2202)
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
      if (kbctyp.eq.9999) then
        if (iexact_trunc .eq. 0 .and. iexact_disc .eq. 0) then
          write(iunit11,2203)
        end if
      end if
c
c     check and correct for reverse input of indicies
      ktemp1 = kbcinfo(nbl,abs(nseg),2,1)
      ktemp2 = kbcinfo(nbl,abs(nseg),3,1)
      if (ktemp1 .gt. ktemp2) then
         kbcinfo(nbl,abs(nseg),2,1) = ktemp2
         kbcinfo(nbl,abs(nseg),3,1) = ktemp1
      end if
      ktemp1 = kbcinfo(nbl,abs(nseg),4,1)
      ktemp2 = kbcinfo(nbl,abs(nseg),5,1)
      if (ktemp1 .gt. ktemp2) then
         kbcinfo(nbl,abs(nseg),4,1) = ktemp2
         kbcinfo(nbl,abs(nseg),5,1) = ktemp1
      end if
c
      kfoseg = 1
      if(nseg.lt.0) kfoseg = 0
      mfoseg = 0
      if(abs(kbctyp).eq.2004 .or. abs(kbctyp).eq.1005 .or.
     .   abs(kbctyp).eq.1006 .or. abs(kbctyp).eq.2014 .or.
     .   abs(kbctyp).eq.2024 .or. abs(kbctyp).eq.2034 .or.
     .   abs(kbctyp).eq.2016) mfoseg = 1
      kfoseg = kfoseg*mfoseg
      nseg = abs(nseg)
      if (ndata.gt.0) then
         read(iunit5,*)
         read(iunit5,*) (realval(mm),mm=1,ndata)
         do mm=1,ndata
            bcvalk(nbl,nseg,mm,1) = realval(mm)
         end do
      end if
      if (ndata.lt.0) then
         nfiles = nfiles+1
         if (nfiles.gt.mxbcfil) then
            write(iunit11,*)'too many bc files specified...increase ',
     .      'parameter mxbcfil'
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
         bcfilek(nbl,nseg,1) = nfiles
         read(iunit5,*)
         read(iunit5,'(a60)') bcfiles(nfiles)
      end if
      if (ig.ne.igrid) then
         write(iunit11,*)' you must put these lines in order by grid!'
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (kbcinfo(nbl,nseg,1,1).ne.-99) then
         write(iunit11,*) 
     .   ' stopping...attempting to set data for segment ',
     .   nseg,' more than once'
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
      kbcinfo(nbl,nseg,1,1) = kbctyp
      if (kbcinfo(nbl,nseg,2,1).eq.0 .and.
     .    kbcinfo(nbl,nseg,3,1).eq.0) then
          kbcinfo(nbl,nseg,2,1) = 1
          kbcinfo(nbl,nseg,3,1) = idimg(nbl)
      end if
      if (kbcinfo(nbl,nseg,4,1).eq.0 .and.
     .    kbcinfo(nbl,nseg,5,1).eq.0) then
          kbcinfo(nbl,nseg,4,1) = 1
          kbcinfo(nbl,nseg,5,1) = jdimg(nbl)
      end if
      if (kfoseg.eq.0 .or. kbctyp.eq.0) then
          kbcinfo(nbl,nseg,6,1) = 0
      end if
      kbcinfo(nbl,nseg,7,1) = ndata
c
      write(iunit11,36) ig,nseg,kbcinfo(nbl,nseg,1,1),
     .  kbcinfo(nbl,nseg,2,1),kbcinfo(nbl,nseg,3,1),
     .  kbcinfo(nbl,nseg,4,1),kbcinfo(nbl,nseg,5,1),
     .  kbcinfo(nbl,nseg,7,1)
      if (ndata.gt.0) then
         call getdhdr(datahdr,kbctyp,ndata)
         write(iunit11,1210) (datahdr(mm),mm=1,ndata)
         write(iunit11,1212) (real(bcvalk(nbl,nseg,mm,1)),mm=1,ndata)
      end if
      if (ndata.lt.0) then
         nfl = bcfilek(nbl,nseg,1)
         write(iunit11,1213)
         write(iunit11,'(''     '',a60)') bcfiles(nfl)
      end if
c
      if (kbcinfo(nbl,nseg,2,1).lt.1 .or.
     .    kbcinfo(nbl,nseg,2,1).gt.idimg(nbl)) then
          write(iunit11,*)' stopping...ista is out of range'
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (kbcinfo(nbl,nseg,3,1).lt.1 .or.
     .    kbcinfo(nbl,nseg,3,1).gt.idimg(nbl)) then
          write(iunit11,*)' stopping...iend is out of range'
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (kbcinfo(nbl,nseg,4,1).lt.1 .or.
     .    kbcinfo(nbl,nseg,4,1).gt.jdimg(nbl)) then
          write(iunit11,*)' stopping...jsta is out of range'
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (kbcinfo(nbl,nseg,5,1).lt.1 .or.
     .    kbcinfo(nbl,nseg,5,1).gt.jdimg(nbl)) then
          write(iunit11,*)' stopping...jend is out of range'
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
 1601 continue
c
c     set boundary condition type array to 21 to allow
c     conservative flux accumulations
c
      if (iemg(igrid).gt.0 .and. ksg(nbl).ne.1) then
c        k=1 face of current block is embedded
         if (nbck0(nbl).gt.1) then
            write(iunit11,*)' error: embedded mesh boundary ',
     .      'at k=1 must extend over entire block face'
            write(iunit11,*)'        segmentation not allowed'
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         else
            kbcinfo(nbl,1,1,1) = 21
            kbcinfo(nbl,1,2,1) = 1
            kbcinfo(nbl,1,3,1) = idimg(nbl)
            kbcinfo(nbl,1,4,1) = 1
            kbcinfo(nbl,1,5,1) = jdimg(nbl)
         end if
      end if
c
      ncg = ncgg(igrid)
      if (ncg.gt.0) then
         do 1679 n=1,ncg
         nbl           = nbl+1
         do 1679 nseg=1,nbck0(nbl)
         kbcinfo(nbl,nseg,1,1) = kbcinfo(nbl-1,nseg,1,1)
c
c        set extrapolation boundary condition on coarser blocks for
c        multigrid chimera. SEE NOTE IN THE I=I0 BOUNDARY SECTION
         if (iovrlp(nbl).lt.0) then
            if (kbcinfo(nbl,nseg,1,1).eq.0) kbcinfo(nbl,nseg,1,1) = 1002
         end if
c
         kbcinfo(nbl,nseg,2,1) = kbcinfo(nbl-1,nseg,2,1)/2+1
         kbcinfo(nbl,nseg,3,1) = kbcinfo(nbl-1,nseg,3,1)/2+1
         kbcinfo(nbl,nseg,4,1) = kbcinfo(nbl-1,nseg,4,1)/2+1
         kbcinfo(nbl,nseg,5,1) = kbcinfo(nbl-1,nseg,5,1)/2+1
         kbcinfo(nbl,nseg,6,1) = kbcinfo(nbl-1,nseg,6,1)
         kbcinfo(nbl,nseg,7,1) = kbcinfo(nbl-1,nseg,7,1)
         bcfilek(nbl,nseg,1)   = bcfilek(nbl-1,nseg,1)
         do 1679 l=1,12
         bcvalk(nbl,nseg,l,1) = bcvalk(nbl-1,nseg,l,1)
 1679    continue
      end if
c
      nblt=(igrid-1)*(ncgg(igrid)+1)+1
      isum=0
      do nseg=1,nbck0(nblt)
        isum=isum+(kbcinfo(nblt,nseg,3,1)-kbcinfo(nblt,nseg,2,1))*
     .            (kbcinfo(nblt,nseg,5,1)-kbcinfo(nblt,nseg,4,1))
      enddo
      if(isum .lt. (idimg(nblt)-1)*(jdimg(nblt)-1)) then
        write(iunit11,'('' Error.  K0 BCs do not span the space'',
     .  '' for grid '',i5)') igrid
        write(iunit11,
     .  '(''   Must exactly cover 1-'',i4,'' in i and 1-'',
     .  i4,'' in j'')') idimg(nblt),jdimg(nblt)
        call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if(isum .gt. (idimg(nblt)-1)*(jdimg(nblt)-1)) then
        write(iunit11,'('' Error.  K0 BCs overspan the space'',
     .  '' for grid '',i5)') igrid
        write(iunit11,
     .  '(''   Must exactly cover 1-'',i4,'' in i and 1-'',
     .  i4,'' in j'')') idimg(nblt),jdimg(nblt)
        call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
 1607 continue
c
c     KDIM Boundary:
c     grid       - grid number (must be in order)
c     kbcinfo(1) - boundary condition type for current segment
c     kbcinfo(2) - i starting location for current segment
c     kbcinfo(3) - i ending location for current segment
c     kbcinfo(4) - j starting location for current segment
c     kbcinfo(5) - j ending location for current segment
c     kbcinfo(6) - > 0...compute force/moment on this segment
c                  < 0...do not compute force/moment on this segment
c     kbcinfo(7) - > 0...number of additional data *values* required for
c                        the boundary condition on this segment
c                  < 0...number of additional data *arrays* required for
c                        the boundary condition on this segment
c                  = 0...no additional data required for this segment
c     bcvalk(m)  - data values for this segment (kbcinfo(7) > 0 only)
c     bcfilek    - pointer to name of data file containing the data
c                  arrays for this segment (kbcinfo(7) < 0 only)
c                  (a value of 1 points to 'null')
c
c***** repeat above for 1 through nbckdim for each grid *****
c
      read(iunit5,10)
      nbl = 0
      write(iunit11,1709)
 1709 format(5hkdim:,1x,4hgrid,3x,7hsegment,4x,6hbctype,6x,4hista,
     .               6x,4hiend,6x,4hjsta,6x,4hjend,5x,5hndata)
      do 1707 igrid=1,ngrid
      nbl = nbl+1
      ifor  = iforce(nbl)
      ifo   = ifor/100
      jfo   = (ifor - ifo*100)/10
      kfo   = (ifor - ifo*100 - jfo*10)
      do 1706 iseg=1,nbckdim(nbl)
      kbcinfo(nbl,iseg,1,2) = -99
      kbcinfo(nbl,iseg,6,2) = 0
      if(kfo.eq.2 .or. kfo.eq.3)kbcinfo(nbl,iseg,6,2) = 1
      bcfilek(nbl,iseg,2) = 1
      do 1706 mm=1,12
      bcvalk(nbl,iseg,mm,2) = -1.e15
 1706 continue
      do 1710 iseg=1,nbckdim(nbl)
      read(iunit5,*) ig,nseg,kbctyp,kbcinfo(nbl,abs(nseg),2,2),
     .                  kbcinfo(nbl,abs(nseg),3,2),
     .                  kbcinfo(nbl,abs(nseg),4,2),
     .                  kbcinfo(nbl,abs(nseg),5,2),ndata
      if (kbctyp .eq. 1004) then
          write(iunit11,2201)
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if ((abs(kbctyp) .eq. 2004 .or. abs(kbctyp) .eq. 2014 .or.
     .     abs(kbctyp) .eq. 2024 .or. abs(kbctyp) .eq. 2034 .or.
     .     abs(kbctyp) .eq. 2016) .and. iviscg(nbl,3) .eq. 0) then
       write(iunit11,'('' stopping...must have ivisc(k) > 0 if'')')
       write(iunit11,'(''   bc2004/14/24/34/16 is being used on'',
     .   '' a k-face'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (abs(kbctyp) .eq. 2014 .and. ivmx .le. 3) then
          write(iunit11,'('' stopping...2014 can only be used'')')
          write(iunit11,'(''   for turb models ivisc>3'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (abs(kbctyp) .eq. 2024 .and. ivmx .ne. 30) then
          write(iunit11,'('' stopping...2024 can only be used'')')
          write(iunit11,'(''   for turb models ivisc=30'')')
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
      if (kbctyp.ge.2000 .and. kbctyp.lt.3000 .and. ndata.eq.0) then
          write(iunit11,2202)
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
      if (kbctyp.eq.9999) then
        if (iexact_trunc .eq. 0 .and. iexact_disc .eq. 0) then
          write(iunit11,2203)
        end if
      end if
c
c     check and correct for reverse input of indicies
      ktemp1 = kbcinfo(nbl,abs(nseg),2,2)
      ktemp2 = kbcinfo(nbl,abs(nseg),3,2)
      if (ktemp1 .gt. ktemp2) then
         kbcinfo(nbl,abs(nseg),2,2) = ktemp2
         kbcinfo(nbl,abs(nseg),3,2) = ktemp1
      end if
      ktemp1 = kbcinfo(nbl,abs(nseg),4,2)
      ktemp2 = kbcinfo(nbl,abs(nseg),5,2)
      if (ktemp1 .gt. ktemp2) then
         kbcinfo(nbl,abs(nseg),4,2) = ktemp2
         kbcinfo(nbl,abs(nseg),5,2) = ktemp1
      end if
c
      kfoseg = 1
      if(nseg.lt.0) kfoseg = 0
      mfoseg = 0
      if(abs(kbctyp).eq.2004 .or. abs(kbctyp).eq.1005 .or.
     .   abs(kbctyp).eq.1006 .or. abs(kbctyp).eq.2014 .or.
     .   abs(kbctyp).eq.2024 .or. abs(kbctyp).eq.2034 .or.
     .   abs(kbctyp).eq.2016) mfoseg = 1
      kfoseg = kfoseg*mfoseg
      nseg = abs(nseg)
      if (ndata.gt.0) then
         read(iunit5,*)
         read(iunit5,*) (realval(mm),mm=1,ndata)
         do mm=1,ndata
            bcvalk(nbl,nseg,mm,2) = realval(mm)
         end do
      end if
      if (ndata.lt.0) then
         nfiles = nfiles+1
         if (nfiles.gt.mxbcfil) then
            write(iunit11,*)'too many bc files specified...increase ',
     .      'parameter mxbcfil'
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
         bcfilek(nbl,nseg,2) = nfiles
         read(iunit5,*)
         read(iunit5,'(a60)') bcfiles(nfiles)
      end if
      if (ig.ne.igrid) then
         write(iunit11,*)' you must put these lines in order by grid!'
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (kbcinfo(nbl,nseg,1,2).ne.-99) then
         write(iunit11,*) 
     .   ' stopping...attempting to set data for segment ',
     .   nseg,' more than once'
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
      kbcinfo(nbl,nseg,1,2) = kbctyp
      if (kbcinfo(nbl,nseg,2,2).eq.0 .and.
     .    kbcinfo(nbl,nseg,3,2).eq.0) then
          kbcinfo(nbl,nseg,2,2) = 1
          kbcinfo(nbl,nseg,3,2) = idimg(nbl)
      end if
      if (kbcinfo(nbl,nseg,4,2).eq.0 .and.
     .    kbcinfo(nbl,nseg,5,2).eq.0) then
          kbcinfo(nbl,nseg,4,2) = 1
          kbcinfo(nbl,nseg,5,2) = jdimg(nbl)
      end if
      if (kfoseg.eq.0 .or. kbctyp.eq.0) then
          kbcinfo(nbl,nseg,6,2) = 0
      end if
      kbcinfo(nbl,nseg,7,2) = ndata
c
      write(iunit11,36) ig,nseg,kbcinfo(nbl,nseg,1,2),
     . kbcinfo(nbl,nseg,2,2),kbcinfo(nbl,nseg,3,2),
     . kbcinfo(nbl,nseg,4,2),kbcinfo(nbl,nseg,5,2),
     . kbcinfo(nbl,nseg,7,2)
      if (ndata.gt.0) then
         call getdhdr(datahdr,kbctyp,ndata)
         write(iunit11,1210) (datahdr(mm),mm=1,ndata)
         write(iunit11,1212) (real(bcvalk(nbl,nseg,mm,2)),mm=1,ndata)
      end if
      if (ndata.lt.0) then
         nfl = bcfilek(nbl,nseg,2)
         write(iunit11,1213)
         write(iunit11,'(''     '',a60)') bcfiles(nfl)
      end if
c
      if (kbcinfo(nbl,nseg,2,2).lt.1 .or.
     .    kbcinfo(nbl,nseg,2,2).gt.idimg(nbl)) then
          write(iunit11,*)' stopping...ista is out of range'
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (kbcinfo(nbl,nseg,3,2).lt.1 .or.
     .    kbcinfo(nbl,nseg,3,2).gt.idimg(nbl)) then
          write(iunit11,*)' stopping...iend is out of range'
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (kbcinfo(nbl,nseg,4,2).lt.1 .or.
     .    kbcinfo(nbl,nseg,4,2).gt.jdimg(nbl)) then
          write(iunit11,*)' stopping...jsta is out of range'
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (kbcinfo(nbl,nseg,5,2).lt.1 .or.
     .    kbcinfo(nbl,nseg,5,2).gt.jdimg(nbl)) then
          write(iunit11,*)' stopping...jend is out of range'
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
 1710 continue
c
c     set boundary condition type array to 21 to allow
c     conservative flux accumulations
c
      if (iemg(igrid).gt.0 .and. keg(nbl).ne.kdimg(nblcg(nbl))) then
c        k=kdim face of current block is embedded
         if (nbckdim(nbl).gt.1) then
            write(iunit11,*)' error: embedded mesh boundary ',
     .      'at k=kdim must extend over entire block face'
            write(iunit11,*)'        segmentation not allowed'
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         else
            kbcinfo(nbl,1,1,2) = 21
            kbcinfo(nbl,1,2,2) = 1
            kbcinfo(nbl,1,3,2) = idimg(nbl)
            kbcinfo(nbl,1,4,2) = 1
            kbcinfo(nbl,1,5,2) = jdimg(nbl)
         end if
      end if
c
      ncg = ncgg(igrid)
      if (ncg.gt.0) then
         do 1780 n=1,ncg
         nbl           = nbl+1
         do 1779 nseg=1,nbckdim(nbl)
         kbcinfo(nbl,nseg,1,2) = kbcinfo(nbl-1,nseg,1,2)
c
c        set extrapolation boundary condition on coarser blocks for
c        multigrid chimera. SEE NOTE IN THE I=I0 BOUNDARY SECTION
         if (iovrlp(nbl).lt.0) then
            if (kbcinfo(nbl,nseg,1,2).eq.0) kbcinfo(nbl,nseg,1,2) = 1002
         end if
c
         kbcinfo(nbl,nseg,2,2) = kbcinfo(nbl-1,nseg,2,2)/2+1
         kbcinfo(nbl,nseg,3,2) = kbcinfo(nbl-1,nseg,3,2)/2+1
         kbcinfo(nbl,nseg,4,2) = kbcinfo(nbl-1,nseg,4,2)/2+1
         kbcinfo(nbl,nseg,5,2) = kbcinfo(nbl-1,nseg,5,2)/2+1
         kbcinfo(nbl,nseg,6,2) = kbcinfo(nbl-1,nseg,6,2)
         kbcinfo(nbl,nseg,7,2) = kbcinfo(nbl-1,nseg,7,2)
         bcfilek(nbl,nseg,2)   = bcfilek(nbl-1,nseg,2)
         do 1779 l=1,12
         bcvalk(nbl,nseg,l,2) = bcvalk(nbl-1,nseg,l,2)
 1779    continue
c
c        reset iovrlp flag for coarser meshes
         if (iovrlp(nbl).lt.0) iovrlp(nbl) = 0
c
c 
 1780    continue
      end if
      nblt=(igrid-1)*(ncgg(igrid)+1)+1
      isum=0
      do nseg=1,nbckdim(nblt)
        isum=isum+(kbcinfo(nblt,nseg,3,2)-kbcinfo(nblt,nseg,2,2))*
     .            (kbcinfo(nblt,nseg,5,2)-kbcinfo(nblt,nseg,4,2))
      enddo
      if(isum .lt. (idimg(nblt)-1)*(jdimg(nblt)-1)) then
        write(iunit11,'('' Error.  KDIM BCs do not span the space'',
     .  '' for grid '',i5)') igrid
        write(iunit11,
     .  '(''   Must exactly cover 1-'',i4,'' in i and 1-'',
     .  i4,'' in j'')') idimg(nblt),jdimg(nblt)
        call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if(isum .gt. (idimg(nblt)-1)*(jdimg(nblt)-1)) then
        write(iunit11,'('' Error.  KDIM BCs overspan the space'',
     .  '' for grid '',i5)') igrid
        write(iunit11,
     .  '(''   Must exactly cover 1-'',i4,'' in i and 1-'',
     .  i4,'' in j'')') idimg(nblt),jdimg(nblt)
        call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
 1707 continue
c
c     if nplot3d < 0, set plot3d output to contain all solid surfaces in
c     the grid if 3D, or all zones in 2D (do same for printed output if 
c     nprint < 0)
c
      if (nplot3d.lt.0 .or. nprint.lt.0) then 
c        just set iptype=0 as a placeholder - overwritten later
         iptype  = 0
         nplt3dtmp = 0
         do 4800 ng=1,ngrid
c
         nbl = nblg(ng)
         idim = idimg(nbl)
         jdim = jdimg(nbl)
         kdim = kdimg(nbl)
c
         if (i2d.eq.0) then
            do 4801 nseg=1,nbci0(nbl)
            if (abs(ibcinfo(nbl,nseg,1,1)).eq.2004 .or. 
     .          abs(ibcinfo(nbl,nseg,1,1)).eq.1005 .or.
     .          abs(ibcinfo(nbl,nseg,1,1)).eq.1006 .or.
     .          abs(ibcinfo(nbl,nseg,1,1)).eq.2014 .or.
     .          abs(ibcinfo(nbl,nseg,1,1)).eq.2024 .or.
     .          abs(ibcinfo(nbl,nseg,1,1)).eq.2034 .or.
     .          abs(ibcinfo(nbl,nseg,1,1)).eq.2016) then
                ii = 1
                nplt3dtmp = nplt3dtmp + 1
                ipl3dtmp(1,nplt3dtmp)   = ng
                ipl3dtmp(2,nplt3dtmp)   = iptype       
                ipl3dtmp(3,nplt3dtmp)   = ii
                ipl3dtmp(4,nplt3dtmp)   = ii
                ipl3dtmp(5,nplt3dtmp)   = 1
                ipl3dtmp(6,nplt3dtmp)   = ibcinfo(nbl,nseg,2,1)
                ipl3dtmp(7,nplt3dtmp)   = ibcinfo(nbl,nseg,3,1)
                ipl3dtmp(8,nplt3dtmp)   = 1
                ipl3dtmp(9,nplt3dtmp)   = ibcinfo(nbl,nseg,4,1)
                ipl3dtmp(10,nplt3dtmp)  = ibcinfo(nbl,nseg,5,1)
                ipl3dtmp(11,nplt3dtmp)  = 1
            end if
 4801       continue
            do 4802 nseg=1,nbcidim(nbl)
            if (abs(ibcinfo(nbl,nseg,1,2)).eq.2004 .or.
     .          abs(ibcinfo(nbl,nseg,1,2)).eq.1005 .or.
     .          abs(ibcinfo(nbl,nseg,1,2)).eq.1006 .or.
     .          abs(ibcinfo(nbl,nseg,1,2)).eq.2014 .or.
     .          abs(ibcinfo(nbl,nseg,1,2)).eq.2024 .or.
     .          abs(ibcinfo(nbl,nseg,1,2)).eq.2034 .or.
     .          abs(ibcinfo(nbl,nseg,1,2)).eq.2016) then
                ii = idim
                nplt3dtmp = nplt3dtmp + 1
                ipl3dtmp(1,nplt3dtmp)   = ng
                ipl3dtmp(2,nplt3dtmp)   = iptype
                ipl3dtmp(3,nplt3dtmp)   = ii
                ipl3dtmp(4,nplt3dtmp)   = ii
                ipl3dtmp(5,nplt3dtmp)   = 1
                ipl3dtmp(6,nplt3dtmp)   = ibcinfo(nbl,nseg,2,2)
                ipl3dtmp(7,nplt3dtmp)   = ibcinfo(nbl,nseg,3,2)
                ipl3dtmp(8,nplt3dtmp)   = 1
                ipl3dtmp(9,nplt3dtmp)   = ibcinfo(nbl,nseg,4,2)
                ipl3dtmp(10,nplt3dtmp)  = ibcinfo(nbl,nseg,5,2)
                ipl3dtmp(11,nplt3dtmp)  = 1
            end if
 4802       continue
            do 4803 nseg=1,nbcj0(nbl)
            if (abs(jbcinfo(nbl,nseg,1,1)).eq.2004 .or.
     .          abs(jbcinfo(nbl,nseg,1,1)).eq.1005 .or.
     .          abs(jbcinfo(nbl,nseg,1,1)).eq.1006 .or.
     .          abs(jbcinfo(nbl,nseg,1,1)).eq.2014 .or.
     .          abs(jbcinfo(nbl,nseg,1,1)).eq.2024 .or.
     .          abs(jbcinfo(nbl,nseg,1,1)).eq.2034 .or.
     .          abs(jbcinfo(nbl,nseg,1,1)).eq.2016) then
                jj = 1
                nplt3dtmp = nplt3dtmp + 1
                ipl3dtmp(1,nplt3dtmp)   = ng
                ipl3dtmp(2,nplt3dtmp)   = iptype
                ipl3dtmp(3,nplt3dtmp)   = jbcinfo(nbl,nseg,2,1)
                ipl3dtmp(4,nplt3dtmp)   = jbcinfo(nbl,nseg,3,1)
                ipl3dtmp(5,nplt3dtmp)   = 1
                ipl3dtmp(6,nplt3dtmp)   = jj
                ipl3dtmp(7,nplt3dtmp)   = jj
                ipl3dtmp(8,nplt3dtmp)   = 1
                ipl3dtmp(9,nplt3dtmp)   = jbcinfo(nbl,nseg,4,1)
                ipl3dtmp(10,nplt3dtmp)  = jbcinfo(nbl,nseg,5,1)
                ipl3dtmp(11,nplt3dtmp)  = 1
            end if
 4803       continue
            do 4804 nseg=1,nbcjdim(nbl)
            if (abs(jbcinfo(nbl,nseg,1,2)).eq.2004 .or.
     .          abs(jbcinfo(nbl,nseg,1,2)).eq.1005 .or.
     .          abs(jbcinfo(nbl,nseg,1,2)).eq.1006 .or.
     .          abs(jbcinfo(nbl,nseg,1,2)).eq.2014 .or.
     .          abs(jbcinfo(nbl,nseg,1,2)).eq.2024 .or.
     .          abs(jbcinfo(nbl,nseg,1,2)).eq.2034 .or.
     .          abs(jbcinfo(nbl,nseg,1,2)).eq.2016) then
                jj = jdim
                nplt3dtmp = nplt3dtmp + 1
                ipl3dtmp(1,nplt3dtmp)   = ng
                ipl3dtmp(2,nplt3dtmp)   = iptype
                ipl3dtmp(3,nplt3dtmp)   = jbcinfo(nbl,nseg,2,2)
                ipl3dtmp(4,nplt3dtmp)   = jbcinfo(nbl,nseg,3,2)
                ipl3dtmp(5,nplt3dtmp)   = 1
                ipl3dtmp(6,nplt3dtmp)   = jj
                ipl3dtmp(7,nplt3dtmp)   = jj
                ipl3dtmp(8,nplt3dtmp)   = 1
                ipl3dtmp(9,nplt3dtmp)   = jbcinfo(nbl,nseg,4,2)
                ipl3dtmp(10,nplt3dtmp)  = jbcinfo(nbl,nseg,5,2)
                ipl3dtmp(11,nplt3dtmp)  = 1
            end if
 4804       continue
            do 4805 nseg=1,nbck0(nbl)
            if (abs(kbcinfo(nbl,nseg,1,1)).eq.2004 .or.
     .          abs(kbcinfo(nbl,nseg,1,1)).eq.1005 .or.
     .          abs(kbcinfo(nbl,nseg,1,1)).eq.1006 .or.
     .          abs(kbcinfo(nbl,nseg,1,1)).eq.2014 .or.
     .          abs(kbcinfo(nbl,nseg,1,1)).eq.2024 .or.
     .          abs(kbcinfo(nbl,nseg,1,1)).eq.2034 .or.
     .          abs(kbcinfo(nbl,nseg,1,1)).eq.2016) then
                kk = 1
                nplt3dtmp = nplt3dtmp + 1
                ipl3dtmp(1,nplt3dtmp)   = ng
                ipl3dtmp(2,nplt3dtmp)   = iptype
                ipl3dtmp(3,nplt3dtmp)   = kbcinfo(nbl,nseg,2,1)
                ipl3dtmp(4,nplt3dtmp)   = kbcinfo(nbl,nseg,3,1)
                ipl3dtmp(5,nplt3dtmp)   = 1
                ipl3dtmp(6,nplt3dtmp)   = kbcinfo(nbl,nseg,4,1)
                ipl3dtmp(7,nplt3dtmp)   = kbcinfo(nbl,nseg,5,1)
                ipl3dtmp(8,nplt3dtmp)   = 1
                ipl3dtmp(9,nplt3dtmp)   = kk
                ipl3dtmp(10,nplt3dtmp)  = kk
                ipl3dtmp(11,nplt3dtmp)  = 1
            end if
 4805       continue
            do 4806 nseg=1,nbckdim(nbl)
            if (abs(kbcinfo(nbl,nseg,1,2)).eq.2004 .or.
     .          abs(kbcinfo(nbl,nseg,1,2)).eq.1005 .or.
     .          abs(kbcinfo(nbl,nseg,1,2)).eq.1006 .or.
     .          abs(kbcinfo(nbl,nseg,1,2)).eq.2014 .or.
     .          abs(kbcinfo(nbl,nseg,1,2)).eq.2024 .or.
     .          abs(kbcinfo(nbl,nseg,1,2)).eq.2034 .or.
     .          abs(kbcinfo(nbl,nseg,1,2)).eq.2016) then
                kk = kdim
                nplt3dtmp = nplt3dtmp + 1
                ipl3dtmp(1,nplt3dtmp)   = ng
                ipl3dtmp(2,nplt3dtmp)   = iptype
                ipl3dtmp(3,nplt3dtmp)   = kbcinfo(nbl,nseg,2,2)
                ipl3dtmp(4,nplt3dtmp)   = kbcinfo(nbl,nseg,3,2)
                ipl3dtmp(5,nplt3dtmp)   = 1
                ipl3dtmp(6,nplt3dtmp)   = kbcinfo(nbl,nseg,4,2)
                ipl3dtmp(7,nplt3dtmp)   = kbcinfo(nbl,nseg,5,2)
                ipl3dtmp(8,nplt3dtmp)   = 1
                ipl3dtmp(9,nplt3dtmp)   = kk
                ipl3dtmp(10,nplt3dtmp)  = kk
                ipl3dtmp(11,nplt3dtmp)  = 1
            end if
 4806       continue
         else
            nplt3dtmp = nplt3dtmp + 1
            ipl3dtmp(1,nplt3dtmp)   = ng
            ipl3dtmp(2,nplt3dtmp)   = iptype
            ipl3dtmp(3,nplt3dtmp)   = 1
            ipl3dtmp(4,nplt3dtmp)   = 1
            ipl3dtmp(5,nplt3dtmp)   = 1
            ipl3dtmp(6,nplt3dtmp)   = 1
            ipl3dtmp(7,nplt3dtmp)   = jdim
            ipl3dtmp(8,nplt3dtmp)   = 1
            ipl3dtmp(9,nplt3dtmp)   = 1
            ipl3dtmp(10,nplt3dtmp)  = kdim
            ipl3dtmp(11,nplt3dtmp)  = 1
         end if
c
 4800    continue
c
c        make sure nplots is large enough 
c
         if (nplt3dtmp.gt.nplots) then
            write(iunit11,1524)
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
c
      end if
c
c     mseq   - mesh sequencing flag for global grids (maximum 5)
c            = 1  single solution on finest grid
c            = 2  solution on second finest grid advanced ncyc(1) cycles
c                 followed by ncyc(2) cycles on finest grid.  The solu-
c                 tion on the finest grid is obtained by interpolation
c                 from the coarser grid.  If ncyc(2)=0, solution 
c                 terminated on second finest grid after ncyc(1) steps
c                 with restart file written for second finest grid at 
c                 that point.
c            > 2  sequencing from coarest to finest mesh as above
c     mgflag - multigrid flag
c            = 0  no multigrid
c            = 1  multigrid on coarser global meshes 
c            = 2  multigrid on coarser global meshes and on 
c                  embedded meshes
c     iconsf - conservation flag
c            = 0  nonconservative flux treatment for embedded grids
c            = 1  conservative flux treatment for embedded grids
c     mtt    = 0  no additional iterations on the "up" portion
c                 of the multigrid cycle
c            > 0  mtt additional iterations on the "up" portion
c                 of the multigrid cycle
c     ngam - multigrid cycle flag
c          = 1  V-cycle
c          = 2  W-cycle
c
      read(iunit5,10)
      read(iunit5,*) mseq,mgflag,iconsf,mtt,ngam
      write(iunit11,7832)
 7832 format(6x,4hmseq,4x,6hmgflag,4x,6hiconsf,7x,3hmtt,6x,4hngam)     
      write(iunit11,36) mseq,mgflag,iconsf,mtt,ngam
c
      if (mseq.gt.5) then
         write(iunit11,'('' Stopping... mseq cannot exceed 5'')')
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
      if (mseq.gt.1) then
         nbl = 0
         do 78 igrid=1,ngrid
         nbl = nbl+1
         if (iovrlp(nbl).gt.0) then
            write(iunit11,*)' Stopping!',
     .                 '  Mesh sequencing with Chimera not allowed!'
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
         ncg = ncgg(igrid)
         if (ncg.gt.0) nbl = nbl+ncg
   78    continue
      end if
      if (mseq.gt.1 .and. iemtot.ne.0 .and. iteravg .ne. 0) then
        write(iunit11,*)' Stopping!',
     .        '  Cannot use iteravg>0 with mseq>1 and embedded grids'
        call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (iexact_trunc .ne. 0 .and. iexact_disc .ne. 0) then
        write(iunit11,*)' Stopping!',
     .      '  cannot have both iexact_trunc and iexact_disc',
     .      ' .ne. 0'
        call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (iexact_trunc .ne. 0) then
        if (ncg .ne. 0 .or. mgflag .ne. 0 .or.
     +      ifullns .ne. 1 .or. i2d .ne. 1) then
          write(iunit11,*)' Stopping!',
     .        '  must have ncg=0,mgflag=0,ifullns=1,i2d=1 when',
     .        ' iexact_trunc.ne.0'
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
        end if
        if (ivmx .ne. 5) then
          write(iunit11,*)' Stopping!',
     .        '  currently must use ivisc=5 for iexact_trunc.ne.0'
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
        end if
        if (iexact_trunc .ne. 1 .and. iexact_trunc .ne. 2 .and.
     .      iexact_trunc .ne. 4) then
          write(iunit11,*)' Stopping!',
     .        '  currently must have iexact_trunc=1,2, or 4'
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
        end if
#   ifdef CMPLX
        write(iunit11,*)' Stopping!',
     .        '  currently cannot use iexact_trunc.ne.0 with complex'
        call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
#   endif
      end if
      if (iexact_disc .ne. 0) then
        if (ifullns .ne. 1 .or. i2d .ne. 1) then
          write(iunit11,*)' Stopping!',
     .        '  must have ifullns=1,i2d=1 when iexact_disc.ne.0'
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
        end if
        if (ivmx .ne. 5) then
          write(iunit11,*)' Stopping!',
     .        '  currently must use ivisc=5 for iexact_disc.ne.0'
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
        if (iexact_disc .ne. 1 .and. iexact_disc .ne. 2 .and.
     .      iexact_disc .ne. 4) then
          write(iunit11,*)' Stopping!',
     .        '  currently must have iexact_disc=1,2, or 4'
          call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
        end if
        end if
#   ifdef CMPLX
        write(iunit11,*)' Stopping!',
     .        '  currently cannot use iexact_disc.ne.0 with complex'
        call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
#   endif
      end if
      if (iexact_ring .ne. 0) then
        if (iexact_trunc .eq. 0 .and. iexact_disc .eq. 0) then
          write(iunit11,*)' iexact_ring must be used with',
     .      ' iexact_trunc=1 or iexact_disc=1'
        call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
        end if
      end if
c
c
c     issc      - correction smoothing flag
c               = 0  off
c               = 1  on
c     epsssc(m) - correction smoothing coefficient for I, J, 
c                 and K directions 
c                 typical values: 0.3, 0.3, 0.3
c     issr      - residual smoothing flag
c               = 0  off
c               = 1  on
c     epsssr(m) - residual smoothing coefficient for I, J,
c                 and K directions 
c                 typical values: 0.3, 0.3, 0.3  (usually not used)
c
      read(iunit5,10)
      read(iunit5,*) issc,(realval(n),n=1,3),issr,(realval(n),n=4,6)
      do n=1,3
         epsssc(n) = realval(n)
         epsssr(n) = realval(n+3)
      end do
      write(iunit11,7932)
 7932 format(6x,4hissc,1x,9hepsssc(1),1x,9hepsssc(2),1x,9hepsssc(3),
     .       6x,4hissr,1x,9hepsssr(1),1x,9hepsssr(2),1x,9hepsssr(3))
      write(iunit11,360) issc,(real(epsssc(n)),n=1,3),issr,
     .                   (real(epsssr(n)),n=1,3)
  360 format(2(i10,3f10.4))
c
c
c     for each sequence 1 through mseq (coarsest to finest):
c     ncyc1  - number of cycles
c     mglevg - number of grids to use in multigrid cycling for 
c              the global meshes
c            = 1 for single grid
c            = 2 for two levels
c            = m for m levels
c     nemgl  - number of embedded grid levels above the finest 
c              global grid (=0 for global grids coarser than the 
c              finest global grid)
c            = 0 no embedded grids
c            = 1 one embedded grid
c            = m m embedded grids
c     nitfo1 - number of first order iterations
c
      ncyctot = 0
      read(iunit5,10)
      write(iunit11,7833)
 7833 format(6x,4hncyc,4x,6hmglevg,5x,5hnemgl,5x,5hnitfo)
      do 12 m=1,mseq
      read(iunit5,*) ncyc1(m),mglevg(m),nemgl(m),nitfo1(m)
      idum = mglevg(m)+nemgl(m)
      if (m.gt.1 .and. ncyc1(m).ne.0 .and. nfreeze.ne.0) then
        write(iunit11,'('' not allowed to freeze turb model'',
     .   '' immediately upon mesh-sequencing-up to next level'')')
        call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
      if (m.gt.1 .and. ncyc1(m).ne.0 .and. iteravg.gt.0) then
        write(iunit11,'('' not allowed to use a running-average'',
     .   '' file (iteravg>0) with mesh sequencing'')')
        call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
c     check input parameters for consistency
c
      if (mglevg(m).le.0 .or. nemgl(m).lt.0 .or. idum.gt.maxbl) then
         write(iunit11,340)
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
  340 format(46h incorrect input value for mglevg and/or nemgl)
c
      if (nemgl(m).gt.0 .and. iemtot.eq.0) then
         write(iunit11,258)
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
  258 format(46h inconsistent values assigned to nemgl and iem)
c
      if (mglevg(m).gt.1 .and. mgflag.lt.1) then
         write(iunit11,268)
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
  268 format(50h inconsistent values assigned to mglevg and mgflag)
c
      if (m.lt.mseq .and. nemgl(m).gt.0) then
         write(iunit11,197)
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
  197 format(50h embedded grids only allowed on finest global grid)
c
      if (real(dt).lt.0.e0) then
         ncyctot = ncyctot+ncyc1(m)
      else
         if (ncyc1(m).gt.0) then
            ncyctot = ncyctot+ntstep
         end if
      end if
c
      if (real(dt).gt.0 .and. mgflag.ne.0 .and. ncyc1(m).eq.1) then
         write(iunit11,482)
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
  482 format(50h must do subits (ncyc>1) with MG and time-accurate)
   12 continue
c
c
      if (mseq.gt.1 .and. ivmx.ge.4) then
         iset=0
         do 2929 iseq=1,mseq
            if (ncyc1(iseq) .gt. 0) iset=iset+1
 2929    continue
         if (iset .gt. 1 .and. irest .ne. 0) then
            isminc=2
         end if
      end if
c
      do 2939 igrid = 1,ngrid
      if (inewgg(igrid).gt.0 .and. ivmx.ge.4) isminc = 2
 2939 continue
c
c     determine levelt and levelb
c
c     levelt - starting level for multigrid/time cycling
c     levelb - ending level for multigrid/time cycling
c
      do 17 m=1,mseq
      levelt(m) = ncgmax-(mseq-m)+nemgl(m)+1
      levelb(m) = levelt(m)-(mglevg(m)-1)-nemgl(m)
      if (levelb(m).lt.1) then
         write(iunit11,157)m,levelt(m),levelb(m),ncgmax
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
  157 format(1x,42herror in input, m, levelt, levelb, ncgmax=,4i5)
      write(iunit11,36)ncyc1(m),mglevg(m),nemgl(m),nitfo1(m)
      if (abs(ita).gt.1 .and. ncyc1(m).eq.1) then
         write(iunit11,198)
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
  198 format(43h you MUST use subiterations (ncyc>1) to get,
     .       31h second-order accuracy in time!)
   17 continue
c
      if (icall.gt.0) then
         if (ncyctot.gt.ncycmax) then
            write(iunit11,893)ncyctot,ncycmax
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
      end if
  893 format(19h number of cycles (,i5,19h) exceeds ncycmax (,i5,1h))
      if (ncyc1(1).le.0) then
         write(iunit11,199)
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
  199 format(47h you MUST do at least one cycle on the starting,
     .       7h level!)
c
c     mitL   - iterations on level L for each level L from coarsest
c               to finest (mitL=1 recommended)
c
      read(iunit5,10)
      write(iunit11,783)
  783 format(6x,4hmit1,6x,4hmit2,6x,4hmit3,6x,4hmit4,6x,9hmit5  ...)
      do 121 m=1,mseq
      read(iunit5,*)  (mit(i,m),i=1,(mglevg(m)+nemgl(m)))
      write(iunit11,36)(mit(i,m),i=1,(mglevg(m)+nemgl(m)))
  121 continue
c
c
c     nbli - number of 1:1 grid-point-connecting block interfaces 
c     NOTE: If nbli is negative, then information for checking the
c           1-1 blocking interfaces between grids is printed to
c           fort.91 and fort.92.
c
      read(iunit5,10)
      read(iunit5,10)
      write(iunit11,9098)
 9098 format(19h 1-1 blocking data:)
      write(iunit11,9099)
 9099 format(6x,4hnbli)
      read(iunit5,*) nbli
      ntest = abs(nbli)
      ntest = ntest + ncgmax*ntest
      if (ntest.gt.mxbli) then
         write(iunit11,9100)
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
 9100 format(47h too many block boundaries specified - increase,
     .6h mxbli)
c
      write(iunit11,36) nbli
c     nblon - block interface coupling flag for 1:1 point connections
c           = 0 on
c           < 0 off (not used)
      if (abs(nbli).gt.0) then
         do 9038 n=1,abs(nbli)
           nblon(n)=0
 9038   continue
      end if 
c
c
c     number     = 1-to-1 interface number; helpful to user only
c     nblk(ib,n) - grid/block identifer for the two blocks that contain
c                  the common interface m
c             ib = 1,2 for first or second block
c             n  = block number
c     limblk(ib,l,n) - i,j,k limits for 1:1 block interface
c                 ib = 1,2 for first or second block
c                 l  = 1-3 (i,j,k limits of the start of the block
c                           interface)
c                      4-6 (i,j,k limits of the end of the block
c                           interface)
c                 n  = block number
c     isva(ib,ind,n) - indicates which indices vary on 1:1 block
c                      interface
c                    = 1 for i varying
c                    = 2 for j varying
c                    = 3 for k varying
c             where:   isva(1,1,n) varies with isva(2,1,n)
c                      isva(1,2,n) varies with isva(2,2,n)
c                ib  = 1,2 for first or second block
c                ind = 1,2
c                n   = block number
c      NOTE:  On input, nblk refers to a grid index.
c             On exit from subroutine global, nblk refers to a block
c             index.
c      NOTE:  On input, limblk refers to grid point indices.
c             On exit from subroutine global, the limblk indices for 
c             the 1:1 patching index loops in the plane refer to 
c             cell-centers
c
c***** repeat above nbli times *****
c
      read(iunit5,10)
      write(iunit11,9047)
 9047 format(2x,6hnumber,4x,4hgrid,4x,4hista,4x,4hjsta,4x,4hksta,
     .       4x,4hiend,4x,4hjend,4x,4hkend,
     .       3x,5hisva1,3x,5hisva2)
      if (abs(nbli).gt.0) then
         do 9033 n=1,abs(nbli)
         read(iunit5,*)mdum,nblk(1,n),
     .     (limblk(1,l,n),l=1,6),(isva(1,ind,n),ind=1,2)
         write(iunit11,37)mdum,nblk(1,n),
     .     (limblk(1,l,n),l=1,6),(isva(1,ind,n),ind=1,2)
c        preliminary 1-1 interface check
         do 9790 m = 1, 3
         if (m + isva(1,1,n) + isva(1,2,n) .eq. 6) then
            if (limblk(1,m,n) .ne. limblk(1,m+3,n)) then
               if (m.eq.1) write(iunit11,9901) 'i'
               if (m.eq.2) write(iunit11,9901) 'j'
               if (m.eq.3) write(iunit11,9901) 'k'
               call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
            end if
         end if
 9790    continue
 9033    continue
      end if
 37   format(10i8)
c
      read(iunit5,10)
      write(iunit11,9147)
 9147 format(2x,6hnumber,4x,4hgrid,4x,4hista,4x,4hjsta,4x,4hksta,
     .       4x,4hiend,4x,4hjend,4x,4hkend,
     .       3x,5hisva1,3x,5hisva2)
      if (abs(nbli).gt.0) then
         do 9133 n=1,abs(nbli)
         read(iunit5,*)mdum,nblk(2,n),
     .     (limblk(2,l,n),l=1,6),(isva(2,ind,n),ind=1,2)
         write(iunit11,37)mdum,nblk(2,n),
     .     (limblk(2,l,n),l=1,6),(isva(2,ind,n),ind=1,2)
c        preliminary 1-1 interface check
         do 9791 m = 1, 3
         if (m + isva(2,1,n) + isva(2,2,n) .eq. 6) then
            if (limblk(2,m,n) .ne. limblk(2,m+3,n)) then
               if (m.eq.1) write(iunit11,9901) 'i'
               if (m.eq.2) write(iunit11,9901) 'j'
               if (m.eq.3) write(iunit11,9901) 'k'
               call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
            end if
         end if
 9791    continue
 9133    continue
      end if
c
 9901 format(34h problem with the above 1-1 input:,
     .       1x,a1,28h not a constant at interface)
c
c
c     ninter - total number of block interpolations
c            = 0  no patching
c            < 0  read interpolation coefficients previously calculated
c                 (input from file patch.bin)
c            > 0  dynamic mesh patching
c
c
      read(iunit5,10)
      read(iunit5,10)
      write(iunit11,1260)
 1260 format(22h patch interface data:)
      write(iunit11,1261)
 1261 format(4x,6hninter)
      read(iunit5,36)ninter
      write(iunit11,36)ninter
c
      if (iunst.eq.0 .and. ninter.gt.0) then
c
c        no longer an option to do generalized patching 
c        internally....need to run preprocessor "ronnie"
c        only the less general case of dynamic patching
c        is done internally
c
         write(iunit11,*)' coefficients for patching must be ',
     .   'precalculated and available in file patch.bin'
         call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
c     determine size requirements for patched applications
c
      if (ninter.eq.0) then
         lnsub1 = 1
         lintmax = 1
         lmaxxe  = 1
      else
         call rpatch0(intmax,nsub1,iindex,ninter)
         rewind(22)
      end if
c
c     for each block, set level and grid from which the block originated
c
      lfem   = 0
      lfgm   = ncgmax+1
      lcgm   = 1
      do 205 igrid=1,ngrid
      iem         = iemg(igrid)
      nbl         = nblg(igrid)
      ncg         = ncgg(igrid)
      igridg(nbl) = igrid
      levelg(nbl) = iem+lfgm
      if (iem.eq.0) lfgm = max(lfgm,levelg(nbl))
      if (iem.gt.0) lfem = max(lfem,levelg(nbl))
      if (ncg.gt.0 .and. iem.eq.0) then
      do 204 n=1,ncg
      nbl         = nbl+1
      igridg(nbl) = igrid
      levelg(nbl) = levelg(nbl-1)-1
  204 continue
      end if
  205 continue
c
c     determine highest global level (levgg) and highest top
c     level (levtt) based on input value of mseq 
c     (top level > global level for embedded meshes)
c
      levtt = levelt(mseq)
      levgg = levtt
      if (nemgl(mseq).gt.0) levgg = levtt - nemgl(mseq)
c
c     determine highest global level actually computed (levglb)
c     and highest top level actually computed (levtop)
c     levels associated with mesh sequence ms are not computed
c     if ncyc(ms) = 0
c
      msmax = 0
      do 177 ms=1,mseq
      if(ncyc1(ms).gt.0) msmax=ms
  177 continue
      levtop = levelt(msmax)
      levglb = levtop
      if(nemgl(msmax).gt.0) levglb = levtop - nemgl(msmax)
c
c***** plot3d output data *****
c
c    igrid  - designated grid number for output
c    iptype - type of plot3d output
c           = 0 grid point type;  grid and q files output
c           = 1 cell center type; grid and q files output
c           = 2 cell center type; grid and turb files output
c     abs() > 2 cell center type; grid and function files output
c             =  3, output smin   (distance to nearest viscous surface)
c             =  4, output vist3d (turbulent eddy viscosity)
c             =  5...cp at cell centers
c             = -5...cp at grid points
c             =  6...p/pinf at grid points
c             = -6...p/pinf at cell centers
c    istart - starting location in I-direction
c    iend   - ending location in I-direction
c    iinc   - increment factor in I-direction
c    jstart - starting location in J-direction
c    jend   - ending location in J-direction
c    jinc   - increment factor in J-direction
c    kstart - starting location in K-direction
c    kend   - ending location in K-direction
c    kinc   - increment factor in K-direction
c
c***** repeated nplot3d times *****
c
      read(iunit5,10)
      read(iunit5,10)
c
      if (nplot3d.ne.0) then
c
         iflag = 0
         ip3dsurf = 0
         if (nplot3d.lt.0) then
            iflag   = 1
            if (i2d.eq.0) ip3dsurf = 1
         end if
         write(iunit11,*) 'plot3d output:'
         write(iunit11,43)
   43    format(2x,4hgrid,1x,5hblock,1x,5hiptyp,2x,4hista,
     .   2x,4hiend,2x,4hiinc,2x,4hjsta,2x,4hjend,2x,
     .   4hjinc,2x,4hksta,2x,4hkend,2x,4hkinc)
c
         np3dhold = abs(nplot3d)
         if (nplot3d.lt.0) then
            nplot3d = nplt3dtmp
         end if
         n1      = 0
c
         np3dread = max(nplot3d,np3dhold)
         do 136 n=1,np3dread
c
         if (iflag.eq.0) then
            read(iunit5,*)(im(l),l=1,11)
            if (im(2) .eq. -1) then
             write(iunit11,*)'stopping...iptyp should not be -1',
     .       ' unless nplot3d=-1 (surface only)'
             call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
            end if
            if (im(2) .eq. 7) then
             write(iunit11,*)'stopping...iptyp should not be 7',
     .       ' unless nplot3d=-1 (surface only)'
             call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
            end if
            if (ivmx .eq. 72 .and. im(2) .eq. 2) then
               write(iunit11,*)'stopping...iptyp=2 not implemented',
     .         ' for ivisc=72... turbulence variables are ouput',
     .         ' automatically to turre...plt'
               call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
            end if
         else
            if (n.le.np3dhold) then
               read(iunit5,*) idum,iptyp,(idum,l=3,11)
            end if
            if (iptyp .eq. 7 .and. i2d .ne. 0) then
             write(iunit11,*)'stopping...iptyp 7',
     .       ' does not work for 2D'
             call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
            end if
            if (ivmx .eq. 72 .and. iptyp .eq. 2) then
               write(iunit11,*)'stopping...iptyp=2 not implemented',
     .         ' for ivisc=72... turbulence variables are ouput',
     .         ' automatically to turre...plt'
               call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
            end if
         end if
c
         if (n.le.nplot3d) then
c
         if (iflag.ne.0) then
            im(1) = ipl3dtmp(1,n)
            im(2) = iptyp
            do 137 l=3,11
            im(l) = ipl3dtmp(l,n)
  137       continue
c           fixup range if cell center, since ipl3dtmp was set
c           assuming grid point type
            if (iptyp.gt.0) then
               idimm = idimg(nblg(im(1)))
               jdimm = jdimg(nblg(im(1)))
               kdimm = kdimg(nblg(im(1)))
               if (im(4).gt.im(3)) then
                  im(4) = im(4) - 1
               else if (im(3).eq.idimm) then
                  im(3) = idimm - 1
                  im(4) = idimm - 1
               end if
               if (im(7).gt.im(6)) then
                  im(7) = im(7) - 1
               else if (im(6).eq.jdimm) then
                  im(6) = jdimm - 1
                  im(7) = jdimm - 1
               end if
               if (im(10).gt.im(9)) then
                  im(10) = im(10) - 1
               else if (im(9).eq.kdimm) then
                  im(9)  = kdimm - 1
                  im(10) = kdimm - 1
               end if
            end if
         end if
c
c        check for inapropriate iptyp
c
         iptyp = im(2)
         if ((iptyp-0)*(iptyp-1)*(iptyp-1)*(iptyp-3)*
     .       (iptyp-4)*(iptyp-5)*(iptyp+5)*(iptyp-6)*
     .       (iptyp+6)*(iptyp-2)*(iptyp+4)*(iptyp+1)*
     .       (iptyp-7) .ne. 0) then
             write(iunit11,*)'stopping...inappropriate value',
     .       ' for iptyp for plot3d data'
             call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
c
         if (ifort50write .ne. 0 .and. 
     +       (ngrid .gt. 1 .or. iptyp .ne. 2)) then
            write(iunit11,*)'stopping...ifort50write must',
     .         ' have ngrid=+-1 and iptyp=2'
               call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
c
         igrid = im(1)
         nbl =  nblg(igrid)
         jd  = jdimg(nbl)
         kd  = kdimg(nbl)
         id  = idimg(nbl)
c
c        set entire index range if all zeros are input
c
         if (im(3).eq.0 .and. im(4).eq.0) then
            im(3) = 1
            im(4) = id
            im(5) = 1
         end if
         if (im(6).eq.0 .and. im(7).eq.0) then
            im(6) = 1
            im(7) = jd
            im(8) = 1
         end if
         if (im(9).eq.0 .and. im(10).eq.0) then
            im(9)  = 1
            im(10) = kd
            im(11) = 1
         end if
c
c        determine block number corresponding to igrid on highest 
c        global level actually computed on (and on embedded levels,
c        if applicable)
c
         nbl = 0
         iflg = 0
         do 178 nnn = 1,nblock
c        levglb: highest global level actually computed on
         if (levelg(nnn).eq.levglb) then
            if (igridg(nnn).eq.igrid) nbl = nnn
         end if
         if (iemg(igrid) .gt.0) then
c           levgg: highest global level in sequence; embedding allowed
c           only above highest global level in sequence
            if (levelg(nnn).gt.levgg) then
               if (igridg(nnn).eq.igrid) nbl = nnn
               if (levtop.lt.levtt) then
                  iflg    = 1
               end if
            end if
         end if
  178    continue
c
         if (nbl.eq.0) then
            write(iunit11,75) igrid
   75       format(3x,39hstopping: no block corresponds to grid ,i3,
     .      33h at the highest level computed on)
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
c
         if (iflg .eq. 0) then
c
            n1 = n1 + 1
            do 657 l = 1,11
            inpl3d(n1,l) = im(l)
657         continue
c
            inpl3d(n1,1) = nbl
            jd  = jdimg(nbl)
            kd  = kdimg(nbl)
            id  = idimg(nbl)
c
c           if necessary, adjust indicies on plot3d input to correspond
c           to the block number consistent with the highest level 
c           actually computed 
c
            if (levglb.lt.levgg) then
c           (if max i-index=2, don't reduce further)
               do 179 mm=3,11
               if(mm.eq.4 .and. inpl3d(n1,mm).eq.2) then
               continue
               else
               inpl3d(n1,mm) = (inpl3d(n1,mm)-1)/2**(levgg-levglb) + 1
               end if
  179          continue
            end if
c
            if (inpl3d(n1,3).le.0)   inpl3d(n1,3)  = 1
            if (inpl3d(n1,3).gt.id)  inpl3d(n1,3)  = id
            if (inpl3d(n1,4).le.0)   inpl3d(n1,4)  = id
            if (inpl3d(n1,4).gt.id)  inpl3d(n1,4)  = id
            if (inpl3d(n1,5).eq.0)   inpl3d(n1,5)  = 1
            if (inpl3d(n1,6).le.0)   inpl3d(n1,6)  = 1
            if (inpl3d(n1,6).gt.jd)  inpl3d(n1,6)  = jd
            if (inpl3d(n1,7).le.0)   inpl3d(n1,7)  = jd
            if (inpl3d(n1,7).gt.jd)  inpl3d(n1,7)  = jd
            if (inpl3d(n1,8).eq.0)   inpl3d(n1,8)  = 1
            if (inpl3d(n1,9).le.0)   inpl3d(n1,9)  = 1
            if (inpl3d(n1,9).gt.kd)  inpl3d(n1,9)  = kd
            if (inpl3d(n1,10).le.0)  inpl3d(n1,10) = kd
            if (inpl3d(n1,10).gt.kd) inpl3d(n1,10) = kd
            if (inpl3d(n1,11).eq.0)  inpl3d(n1,11) = 1
c
            if (i2d.eq.1) then
               inpl3d(n1,3) = 1
               inpl3d(n1,4) = 1
               inpl3d(n1,5) = 1
            end if
c
            write(iunit11,76)igrid,(inpl3d(n1,l),l=1,11)
   76       format(i6,11i6)
c
         else
c
            write(iunit11,1788) igrid
1788        format(4x,46h stopping: embedded grids are computed only on,
     .      13h finest level,/,4x,33h the way you are mesh sequencing,,
     .      30h can't output solution on grid,i3)
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
c
         end if
c
         end if
c
  136    continue
c
c        if any zones have function file output instead of q file
c        output, *all* zones must have function file output
c
         ifunc = 0
         if (nplot3d.gt.0) then
            do nnn = 1,nplot3d
               if (inpl3d(nnn,2).gt.2) ifunc = inpl3d(nnn,2)
            enddo
            if (ifunc.gt.0) then
               do nnn = 1,nplot3d
                  if(inpl3d(nnn,2).ne.ifunc) then
                  write(iunit11,*)'stopping...must have iptype same ',
     .            'for all zones if any zone has iptype > 2'
                  call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
                  endif
               enddo
            end if 
         end if
c
      end if
c
c     movie - animation output flag
c            = 0 plot3d file written only at end of computation
c            > 0 plot3d file written every movie iterations
c            < 0 plot3d file written at start of computation
c                AND plot3d file written every movie iterations 
c
      read(iunit5,*)
      read(iunit5,*) movie
      write(iunit11,1302)
 1302 format(1x,5hmovie)
      write(iunit11,76) movie
c
c***** printed output data *****
c
c    igrid  - designated grid number for output
c    iptype - type of print out output
c           = 0 grid point type
c           = 1 cell center type
c    istart - starting location in I-direction
c    iend   - ending location in I-direction
c    iinc   - increment factor in I-direction
c    jstart - starting location in J-direction
c    jend   - ending location in J-direction
c    jinc   - increment factor in J-direction
c    kstart - starting location in K-direction
c    kend   - ending location in K-direction
c    kinc   - increment factor in K-direction
c
c***** repeated nprint times *****
c
      read(iunit5,10)
      read(iunit5,10)
c
      if (nprint.ne.0) then
c
         iflag = 0
         iprnsurf = 0
         if (nprint.lt.0) then
            iflag  = 1
            if (i2d.eq.0) iprnsurf = 1
         end if
         write(iunit11,*) 'print out:'
         write(iunit11,43)
c
         if (nprint.lt.0) then
            nprnthold = abs(nprint)
            nprint = nplt3dtmp
         end if
         n1     = 0
c
         do 236 n=1,nprint
c
         if (iflag.eq.0) then
            read(iunit5,*)(im(l),l=1,11)
         else
            if (n.le.nprnthold) then
               read(iunit5,*) idum,iptyp,(idum,l=3,11)
            end if
         end if
c
         if (im(5).lt.0 .or. im(8).lt.0 .or. im(11).lt.0) then
           write(iunit11,'('' input error: inc < 0'')')
           call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
         if (im(3) .gt. im(4) .or. im(6) .gt. im(7) .or.
     .       im(9) .gt. im(10)) then
           write(iunit11,'('' input error: end > start'')')
           call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
c
         if (iflag.ne.0) then
            im(1) = ipl3dtmp(1,n)
            im(2) = iptyp
            do 237 l=3,11
            im(l) = ipl3dtmp(l,n)
  237       continue
c           fixup range if cell center, since ipl3dtmp was set
c           assuming grid point type
            if (iptyp.gt.0) then
               idimm = idimg(nblg(im(1)))
               jdimm = jdimg(nblg(im(1)))
               kdimm = kdimg(nblg(im(1)))
               if (im(4).gt.im(3)) then
                  im(4) = im(4) - 1
               else if (im(3).eq.idimm) then
                  im(3) = idimm - 1
                  im(4) = idimm - 1
               end if
               if (im(7).gt.im(6)) then
                  im(7) = im(7) - 1
               else if (im(6).eq.jdimm) then
                  im(6) = jdimm - 1
                  im(7) = jdimm - 1
               end if
               if (im(10).gt.im(9)) then
                  im(10) = im(10) - 1
               else if (im(9).eq.kdimm) then
                  im(9)  = kdimm - 1
                  im(10) = kdimm - 1
               end if
            end if
         end if
c
c        check for inapropriate iptyp
c
         iptyp = im(2)
         if ((iptyp-0)*(iptyp-1) .ne. 0) then
             write(iunit11,*)'stopping...inappropriate value',
     .       ' for iptyp for printout data'
             call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
c
         igrid = im(1)
         nbl =  nblg(igrid)
         jd  = jdimg(nbl)
         kd  = kdimg(nbl)
         id  = idimg(nbl)
c
c        set entire index range if all zeros are input
c
         if (im(3).eq.0 .and. im(4).eq.0) then
            im(3) = 1
            im(4) = id
            im(5) = 1
         end if
         if (im(6).eq.0 .and. im(7).eq.0) then
            im(6) = 1
            im(7) = jd
            im(8) = 1
         end if
         if (im(9).eq.0 .and. im(10).eq.0) then
            im(9)  = 1
            im(10) = kd
            im(11) = 1
         end if
c
c        determine block number corresponding to igrid on highest 
c        global level actually computed on (and on embedded levels,
c        if applicable)
c
         nbl = 0
         iflg = 0
         do 278 nnn = 1,nblock
c        levglb: highest global level actually computed on
         if (levelg(nnn).eq.levglb) then
            if (igridg(nnn).eq.igrid) nbl = nnn
         end if
         if (iemg(igrid) .gt.0) then
c           levgg: highest global level in sequence; embedding allowed
c           only above highest global level in sequence
            if (levelg(nnn).gt.levgg) then
               if (igridg(nnn).eq.igrid) nbl = nnn
               if (levtop.lt.levtt) then
                  iflg    = 1
               end if
            end if
         end if
  278    continue
c
         if (nbl.eq.0) then
            write(iunit11,75) igrid
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
c
         if (iflg .eq. 0) then
c
            n1 = n1 + 1
            do 658 l = 1,11
            inpr(n1,l) = im(l)
658         continue
c
            inpr(n1,1) = nbl
            jd  = jdimg(nbl)
            kd  = kdimg(nbl)
            id  = idimg(nbl)
c
c           if necessary, adjust indicies on print input to correspond
c           to the block number consistent with the highest level
c           actually computed
c
            if (nemgl(mseq).gt.0) levgg = levtt - nemgl(mseq)
            if (levglb.lt.levgg) then
c           (if max i-index=2, don't reduce further)
               do 279 mm=3,11
               if(mm.eq.4 .and. inpr(n1,mm).eq.2) then
               continue
               else
               inpr(n1,mm) = (inpr(n1,mm)-1)/2**(levgg-levglb) + 1
               end if
  279          continue
            end if
c
            if (inpr(n1,3).le.0)   inpr(n1,3)  = 1
            if (inpr(n1,3).gt.id)  inpr(n1,3)  = id
            if (inpr(n1,4).le.0)   inpr(n1,4)  = id
            if (inpr(n1,4).gt.id)  inpr(n1,4)  = id
            if (inpr(n1,5).eq.0)   inpr(n1,5)  = 1
            if (inpr(n1,6).le.0)   inpr(n1,6)  = 1
            if (inpr(n1,6).gt.jd)  inpr(n1,6)  = jd
            if (inpr(n1,7).le.0)   inpr(n1,7)  = jd
            if (inpr(n1,7).gt.jd)  inpr(n1,7)  = jd
            if (inpr(n1,8).eq.0)   inpr(n1,8)  = 1
            if (inpr(n1,9).le.0)   inpr(n1,9)  = 1
            if (inpr(n1,9).gt.kd)  inpr(n1,9)  = kd
            if (inpr(n1,10).le.0)  inpr(n1,10) = kd
            if (inpr(n1,10).gt.kd) inpr(n1,10) = kd
            if (inpr(n1,11).eq.0)  inpr(n1,11) = 1
c
            if (i2d.eq.1) then
               inpr(n1,3) = 1
               inpr(n1,4) = 1
               inpr(n1,5) = 1
            end if
c
            write(iunit11,76)igrid,(inpr(n1,l),l=1,11)
c
         else
c
            write(iunit11,1788) igrid
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
c
         end if
c
  236    continue
c
      end if
c
c***** contol surface data  *****
c
c    ncs    - the number of control surfaces
c
c    igrid  - designated grid number for the control surface
c    istart - starting location in I-direction
c    iend   - ending location in I-direction
c    jstart - starting location in J-direction
c    jend   - ending location in J-direction
c    kstart - starting location in K-direction
c    kend   - ending location in K-direction
c    iwall  - 0 for a flow surface
c             1 for a solid wall
c    inorm  - defines the direction of the control surface outward normal
c             relative to the grid surface normal. This is only required
c             in the calculation of the total forces for all control
c             surfaces to obtain the correct signs
c             1 - the control volume outward normal is in the same 
c                 direction as the grid surface normal
c            -1 - the control volume outward normal is in the opposite
c                 direction to the grid surface normal
c             0 - do not add this surface into the summation of the 
c                 total forces.  Note, if all the surfaces have inorm=0,
c                 the totals are not calculated.
c
c***** repeated ncs times *****
c
      read(iunit5,10)
      read(iunit5,10)
      read(iunit5,*)ncs
      read(iunit5,10)
c
      ncs2    = 0
c
      if (ncs.ne.0) then
c
         write(iunit11,*) 'control surfaces:'
         write(iunit11,*) '  ncs'
         write(iunit11,820)ncs
 820     format(i6)
         write(iunit11,830)
 830     format(2x,4hgrid,1x,5hblock,2x,4hista,2x,4hiend,2x,4hjsta,2x,
     .          4hjend,2x,4hksta,2x,4hkend,1x,5hiwall,1x,5hinorm)
c
         ncs1    = ncs
         ncs     = 0
         ntestcs = 0
         do 800 n=1,ncs1
         read(iunit5,*)(im(l),l=1,9)
c
         igrid = im(1)
c
c
c        set entire index range if all zeros are input
c
         nbl = nblg(igrid)
         jd  = jdimg(nbl)
         kd  = kdimg(nbl)
         id  = idimg(nbl)
c
         if (im(2).eq.0 .and. im(3).eq.0) then
            im(2) = 1
            im(3) = id
         end if
         if (im(4).eq.0 .and. im(5).eq.0) then
            im(4) = 1
            im(5) = jd
         end if
         if (im(6).eq.0 .and. im(7).eq.0) then
            im(6)  = 1
            im(7) = kd
         end if
c
         ntestcs = ntestcs + (ncgg(igrid)+1)
c
         if (icall .gt. 0) then
            if (ntestcs .gt. maxcs-1) then
               write(iunit11,*)
     .         ' stopping...to accomodate fine and coarse',
     .         ' levels, for an input ncs = ',ncs1
               write(iunit11,*)' the parameter maxcs must be at least ',
     .         ncs1*(ncgg(igrid)+1) + 1
               call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
            end if
         end if
c
         do 805 ncg=1,ncgg(igrid)+1
c
         nbl = nblg(igrid) + ncg - 1
         ncs  = ncs + 1
         ncs2 = ncs2 + 1
         if (icall.eq.0) then
            ncs = min(maxcs,ncs)
         end if
         icsinfo(ncs,1) = nbl
         do 379 mm=2,7
  379    icsinfo(ncs,mm) = (im(mm)-1)/2**(ncg-1) + 1
         icsinfo(ncs,8)  = im(8)
         icsinfo(ncs,9)  = im(9)
c
         jd  = jdimg(nbl)
         kd  = kdimg(nbl)
         id  = idimg(nbl)
c
         if (icsinfo(ncs,2).le.0)  icsinfo(ncs,2)  = 1
         if (icsinfo(ncs,2).gt.id) icsinfo(ncs,2)  = id
         if (icsinfo(ncs,3).le.0)  icsinfo(ncs,3)  = id
         if (icsinfo(ncs,3).gt.id) icsinfo(ncs,3)  = id
         if (icsinfo(ncs,4).le.0)  icsinfo(ncs,4)  = 1
         if (icsinfo(ncs,4).gt.jd) icsinfo(ncs,4)  = jd
         if (icsinfo(ncs,5).le.0)  icsinfo(ncs,5)  = jd
         if (icsinfo(ncs,5).gt.jd) icsinfo(ncs,5)  = jd
         if (icsinfo(ncs,6).le.0)  icsinfo(ncs,6)  = 1
         if (icsinfo(ncs,6).gt.kd) icsinfo(ncs,6)  = kd
         if (icsinfo(ncs,7).le.0)  icsinfo(ncs,7)  = kd
         if (icsinfo(ncs,7).gt.kd) icsinfo(ncs,7)  = kd
c
         if (i2d.eq.1) then
            icsinfo(ncs,2) = 1
            icsinfo(ncs,3) = 2
         end if
c
         if (ncg.eq.1) then
            write(iunit11,840)igrid,(icsinfo(ncs,l),l=1,9)
  840       format(10i6)
         end if
c
  805    continue
c
  800    continue
c
c
      end if
c
c***** input for rigid grid motion  *****
c
c     grid    - grid number
c             > 0...input corresponds to grid
c             = 0...input corresponds to moment center
c
c     itrans  - type of modulation for translation
c             = 0...no translation
c             = 1...constant velocity
c             = 2...sinusoidal variation displacement
c             = 3...smooth increase in displacement, asypmtotically
c                   reaching a fixed displacement
c
c     irotat  - type of modulation for rotation
c             = 0...no rotation
c             = 1...constant rotation speed
c             = 2...sinusoidal variation of rotational displacement
c             = 3...smooth inrease in rotational displacement,
c                   asypmtotically reaching a fixed rotational
c                   displacement
c     idefrm  - deforming mesh flag
c             = 0...rigid mesh
c
c      rfreq  = reduced frequency for sinusoidal variation 
c
c      lref   = grid-equivalent of the dimensional reference length that
c               was used to define the reduced frequency, or the 
c               non-dimensional rotation speed input with irotat=1 
c
c      utrans = if itrans=1, x-component of translational velocity
c               if itrans>1, maximum displacement in y-direction
c
c      vtrans = if itrans=1, y-component of translational velocity
c               if itrans>1, maximum displacement in y-direction 
c
c      wtrans = if itrans=1, z-component of translational velocity
c               if itrans>1, maximum displacement in z-direction
c
c      dxmx   = maximum x-displacement of grid before being reset to 
c               t=0 position; if zero, then x-position of grid is 
c               never reset 
c
c      dymx   = maximum y-displacement of grid before being reset to 
c               t=0 position; if zero, then y-position of grid is
c               never reset 
c
c      dzmx   = maximum z-displacement of grid before being reset to 
c               t=0 position; if zero, then z-position of grid is
c               never reset 
c
c      omegax = if irotat=1, x-component of rotational velocity
c               if irotat>1, maximum angular rotation about the x-axis
c
c      omegay = if irotat=1, y-component of rotational velocity
c               if irotat>1, maximum angular rotation about the y-axis
c
c      omegaz = if irotat=1, y-component of rotational velocity
c               if irotat>1, maximum angular rotation about the z-axis
c
c      dthxmx = maximum angular displacement of grid about x-axis before
c               being reset to t=0 position; if zero, then angular 
c               position about the x-axis is never reset
c
c      dthymx = maximum angular displacement of grid about y-axis before
c               being reset to t=0 position; if zero, then angular 
c               position about the y-axis is never reset
c
c      dthzmx = maximum angular displacement of grid about z-axis before
c               being reset to t=0 position; if zero, then angular 
c               position about the z-axis is never reset
c
c      itransmc,irotatmc,utransmc...dthymxmc,dthzmxmc: 
c               translational/rotational parameters for the motion
c               of the moment center, with definitions analogous
c               to itrans,irotat,utrans...above. The moment center 
c               motion is input in the same way as grid motion, with 
c               grid=0 in the input deck. NOTE: moment center is 
c               stationary if itransmc = irotatmc = 0. The moment center
c               parameters are global (not block) values.
c
c NOTE: at the end (after both rigid and deforming/aeroelastic mesh data
c       are set), if ialph = 1 we swap y and z grid motion parameters
c       so that they remain consitant with the grid orientation as it
c       also gets swapped around later in the code if ialph=1
c
c****************************************
c
c     set default values to no motion
c
      itransmc = 0
      irotatmc = 0
      rfreqtmc = 0.
      rfreqrmc = 0.
      utransmc = 0.
      vtransmc = 0.
      wtransmc = 0.
      omegaxmc = 0.
      omegaymc = 0.
      omegazmc = 0.
      thetaxmc = 0.
      thetaymc = 0.
      thetazmc = 0.
      xorigmc  = 0.
      yorigmc  = 0.
      zorigmc  = 0.
      xorig0mc = 0.
      yorig0mc = 0.
      zorig0mc = 0.
      dxmxmc   = 0.
      dymxmc   = 0.
      dzmxmc   = 0.
      dthxmxmc = 0.
      dthymxmc = 0.
      dthzmxmc = 0.
      time2mc  = 0.
      xmc0     = 0.
      ymc0     = 0.
      zmc0     = 0.
c
      do 611 igrid=1,ngrid
      nbl = nblg(igrid)
      ncg = ncgg(igrid)
      itrans(nbl)  = 0
      utrans(nbl)  = 0.
      vtrans(nbl)  = 0.
      wtrans(nbl)  = 0.
      dxmx(nbl)    = 0.
      dymx(nbl)    = 0.
      dzmx(nbl)    = 0.
      irotat(nbl)  = 0
      omegax(nbl)  = 0.
      omegay(nbl)  = 0.
      omegaz(nbl)  = 0.
      xorig(nbl)   = 0.
      yorig(nbl)   = 0.
      zorig(nbl)   = 0.
      xorig0(nbl)  = 0.
      yorig0(nbl)  = 0.
      zorig0(nbl)  = 0.
      thetax(nbl)  = 0.
      thetay(nbl)  = 0.
      thetaz(nbl)  = 0.
      dthxmx(nbl)  = 0.
      dthymx(nbl)  = 0.
      dthzmx(nbl)  = 0.
      rfreqt(nbl)  = 0.
      rfreqr(nbl)  = 0.
      time2(nbl)   = 0.
      thetaxl(nbl) = 0.
      thetayl(nbl) = 0.
      thetazl(nbl) = 0.
c
      if (ncg.gt.0) then
         do 650 n=1,ncg
         nbl          = nbl+1
         itrans(nbl)  = itrans(nbl-1) 
         utrans(nbl)  = utrans(nbl-1)
         vtrans(nbl)  = vtrans(nbl-1)
         wtrans(nbl)  = wtrans(nbl-1)
         irotat(nbl)  = irotat(nbl-1)
         omegax(nbl)  = omegax(nbl-1)
         omegay(nbl)  = omegay(nbl-1)
         omegaz(nbl)  = omegaz(nbl-1)
         xorig(nbl)   = xorig(nbl-1)
         yorig(nbl)   = yorig(nbl-1)
         zorig(nbl)   = zorig(nbl-1)
         xorig0(nbl)  = xorig0(nbl-1)
         yorig0(nbl)  = yorig0(nbl-1)
         zorig0(nbl)  = zorig0(nbl-1)
         thetax(nbl)  = thetax(nbl-1)
         thetay(nbl)  = thetay(nbl-1)
         thetaz(nbl)  = thetaz(nbl-1)
         rfreqt(nbl)  = rfreqt(nbl-1)
         rfreqr(nbl)  = rfreqr(nbl-1)
         dxmx(nbl)    = dxmx(nbl-1)
         dymx(nbl)    = dymx(nbl-1)
         dzmx(nbl)    = dzmx(nbl-1)
         dthxmx(nbl)  = dthxmx(nbl-1)
         dthymx(nbl)  = dthymx(nbl-1)
         dthzmx(nbl)  = dthzmx(nbl-1)
         time2(nbl)   = time2(nbl-1)
         thetaxl(nbl) = thetaxl(nbl-1)
         thetayl(nbl) = thetayl(nbl-1)
         thetazl(nbl) = thetazl(nbl-1)
  650    continue
      end if
  611 continue
c
      if (iunst.eq.1 .or. iunst.eq.3) then
c
c        translation data
c
         read(iunit5,*)
         read(iunit5,*)
         read(iunit5,*) ntrans
         write(iunit11,600)
         write(iunit11,602)
         write(iunit11,604) ntrans
         read(iunit5,*)
         write(iunit11,606)
         if (ntrans .gt. 0) then
            read(iunit5,*) realval(1)
            lreftra = realval(1)
            write(iunit11,608) real(lreftra)
         end if
         read(iunit5,*)
         write(iunit11,612)
         if (ntrans .gt. 0) then
            ichk1 = 0
            if (ntrans .gt. ngrid) then
               if (ntrans .eq. ngrid+1) then
                  ichk1 = 1
               else
                  write(iunit11,*)
     .            ' stopping...must have ntrans <= ngrid + 1'
                  call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
               end if
            end if
            do 610 igrid=1,ntrans
            read(iunit5,*) igr,itran,(realval(i),i=1,4)
            rfreq = realval(1)
            utran = realval(2)
            vtran = realval(3)
            wtran = realval(4)
            if (igr .gt. 0) then
               write(iunit11,614) igr,itran,real(rfreq),real(utran),
     .                            real(vtran),real(wtran)
               nbl = nblg(igr)
               ncg = ncgg(igr)
               itrans(nbl) = itran
               utrans(nbl) = utran
               vtrans(nbl) = vtran
               wtrans(nbl) = wtran
               if (itran.eq.1) lreftra = 1.
               rfreqt(nbl) = 2.*pi*rfreq/lreftra
               if (ncg.gt.0) then
                  do 630 n=1,ncg
                  nbl         = nbl+1
                  itrans(nbl) = itrans(nbl-1)
                  utrans(nbl) = utrans(nbl-1)
                  vtrans(nbl) = vtrans(nbl-1)
                  wtrans(nbl) = wtrans(nbl-1)
                  rfreqt(nbl) = rfreqt(nbl-1)
  630             continue
               end if
            end if
            if (igr .eq. 0) then
               write(iunit11,616) itran,rfreq,utran,vtran,wtran
               itransmc = itran
               utransmc = utran
               vtransmc = vtran
               wtransmc = wtran
               rfreqtmc = 2.*pi*rfreq/lreftra
               ichk1    = 0
            end if
  610       continue
            if (ichk1 .ne. 0) then
               write(iunit11,*) ' stopping...since ntrans = ngrid + 1,',
     .                     ' must have grid = 0 (moment center) in'
               write(iunit11,*) ' one input line of translation data'
               call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
            end if
         end if
         read(iunit5,*)
         write(iunit11,618)
         if (ntrans .gt. 0) then
            do 710 igrid=1,ntrans
            read(iunit5,*) igr,(realval(i),i=1,3)
            dxmax = realval(1)
            dymax = realval(2)
            dzmax = realval(3)
            if (igr .gt. 0) then
               write(iunit11,604) igr,real(dxmax),real(dymax),
     .                            real(dzmax)
               nbl = nblg(igr)
               ncg = ncgg(igr)
               dxmx(nbl) = dxmax
               dymx(nbl) = dymax
               dzmx(nbl) = dzmax
               if (ncg.gt.0) then
                  do 730 n=1,ncg
                  nbl         = nbl+1
                  dxmx(nbl) = dxmx(nbl-1)
                  dymx(nbl) = dymx(nbl-1)
                  dzmx(nbl) = dzmx(nbl-1)
  730             continue
               end if
            end if
            if (igr .eq. 0) then
               write(iunit11,622) real(dxmax),real(dymax),real(dzmax)
               dxmxmc = dxmax
               dymxmc = dymax
               dzmxmc = dzmax
            end if
  710       continue
         end if
c
c        rotation data
c
         read(iunit5,*)
         read(iunit5,*)
         read(iunit5,*) nrotat
         write(iunit11,624)
         write(iunit11,626) 
         write(iunit11,604) nrotat
         read(iunit5,*)
         write(iunit11,606)
         if (nrotat .gt. 0) then
            read(iunit5,*) realval(1)
            lrefrot = realval(1)
            write(iunit11,608) real(lrefrot)
         end if
         read(iunit5,*)
         write(iunit11,628)
         if (nrotat .gt. 0) then
            ichk1 = 0
            if (nrotat .gt. ngrid) then
               if (nrotat .eq. ngrid+1) then
                  ichk1 = 1
               else
                  write(iunit11,*)
     .            ' stopping...must have nrotat <= ngrid + 1'
                  call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
               end if
            end if
            do 620 igrid=1,nrotat
            read(iunit5,*) igr,irota,(realval(i),i=1,7)
            rfreq = realval(1)
            omegx = realval(2) 
            omegy = realval(3)
            omegz = realval(4)
            xorg  = realval(5)
            yorg  = realval(6)
            zorg  = realval(7)
            if (igr .gt. 0) then
               write(iunit11,614) igr,irota,real(rfreq),real(omegx),
     .                            real(omegy),real(omegz),real(xorg),
     .                            real(yorg),real(zorg)
               nbl = nblg(igr)
               ncg = ncgg(igr)
               irotat(nbl) = irota
               xorig(nbl)  = xorg
               yorig(nbl)  = yorg
               zorig(nbl)  = zorg
               if (irota.eq.1) then
                  omegax(nbl) = 2.*pi*omegx/lrefrot
                  omegay(nbl) = 2.*pi*omegy/lrefrot
                  omegaz(nbl) = 2.*pi*omegz/lrefrot
               end if
               if (irota.eq.2 .or. irota.eq.3) then
                  omegax(nbl) = omegx/radtodeg
                  omegay(nbl) = omegy/radtodeg
                  omegaz(nbl) = omegz/radtodeg
               end if
               rfreqr(nbl) = 2.*pi*rfreq/lrefrot
               if (ncg.gt.0) then
                  do 640 n=1,ncg
                  nbl         = nbl+1
                  irotat(nbl) = irotat(nbl-1)
                  omegax(nbl) = omegax(nbl-1)
                  omegay(nbl) = omegay(nbl-1)
                  omegaz(nbl) = omegaz(nbl-1)
                  xorig(nbl)  = xorig(nbl-1)
                  yorig(nbl)  = yorig(nbl-1)
                  zorig(nbl)  = zorig(nbl-1)
                  thetax(nbl) = thetax(nbl-1)
                  thetay(nbl) = thetay(nbl-1)
                  thetaz(nbl) = thetaz(nbl-1)
                  rfreqr(nbl) = rfreqr(nbl-1)
  640             continue
               end if
            end if
            if (igr .eq. 0) then
               write(iunit11,616) irota,real(rfreq),real(omegx),
     .                            real(omegy),real(omegz),real(xorg),
     .                            real(yorg),real(zorg)
               irotatmc = irota
               xorigmc  = xorg
               yorigmc  = yorg
               zorigmc  = zorg
               if (irota.eq.1) then
                  omegaxmc = 2.*pi*omegx/lrefrot
                  omegaymc = 2.*pi*omegy/lrefrot
                  omegazmc = 2.*pi*omegz/lrefrot
               end if
               if (irota.eq.2 .or. irota.eq.3) then
                  omegaxmc = omegx/radtodeg
                  omegaymc = omegy/radtodeg
                  omegazmc = omegz/radtodeg
               end if
               rfreqrmc = 2.*pi*rfreq/lrefrot
               ichk1    = 0
            end if
  620       continue
            if (ichk1 .ne. 0) then
               write(iunit11,*) ' stopping...since nrotat = ngrid + 1,',
     .                     ' must have grid = 0 (moment center) in'
               write(iunit11,*) ' one input line of rotation data'
               call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
            end if
         end if
         read(iunit5,*)
         write(iunit11,632)
         if (nrotat .gt. 0) then
            do 720 igrid=1,nrotat
            read(iunit5,*) igr,(realval(i),i=1,3)
            dtxmax = realval(1)
            dtymax = realval(2)
            dtzmax = realval(3)
            if (igr .gt. 0) then
               write(iunit11,604) igr,real(dtxmax),real(dtymax),
     .                            real(dtzmax)
               nbl = nblg(igr)
               ncg = ncgg(igr)
               dthxmx(nbl) = dtxmax/radtodeg
               dthymx(nbl) = dtymax/radtodeg
               dthzmx(nbl) = dtzmax/radtodeg
               if (ncg.gt.0) then
                  do 740 n=1,ncg
                  nbl         = nbl+1
                  dthxmx(nbl) = dthxmx(nbl-1)
                  dthymx(nbl) = dthymx(nbl-1)
                  dthzmx(nbl) = dthzmx(nbl-1)
  740             continue
               end if
            end if
            if (igr .eq. 0) then
               write(iunit11,622) real(dtxmax),real(dtymax),real(dtzmax)
               dthxmxmc = dtxmax/radtodeg
               dthymxmc = dtymax/radtodeg
               dthzmxmc = dtzmax/radtodeg
            end if
  720       continue
         end if
c
c        store off t=0 values of xorig,yorig,zorig, and also xmc,ymc,zmc
c        
         xmc0 = xmc
         ymc0 = ymc
         zmc0 = zmc
         xorig0mc = xorigmc
         yorig0mc = yorigmc
         zorig0mc = zorigmc
         do 760 n=1,nblock
         xorig0(n) = xorig(n)
         yorig0(n) = yorig(n)
         zorig0(n) = zorig(n)
  760    continue
c
  600    format(33h dynamic grid data - translation:)
  602    format(8h  ntrans)
  604    format(5x,i3,8(1x,f8.4))
  606    format(8h    lref)
  608    format(f8.4)
  612    format(4x,4hgrid,3x,6hitrans,3x,6hrfreqt,3x,6hutrans,
     .          3x,6hvtrans,3x,6hwtrans)
  614    format(5x,i3,7x,i2,8(1x,f8.4))
  616    format(6x,2hMC,7x,i2,8(1x,f8.4))
  618    format(4x,4hgrid,5x,4hdxmx,5x,4hdymx,5x,4hdzmx)
  622    format(6x,2hMC,8(1x,f8.4))
  624    format(30h dynamic grid data - rotation:)
  626    format(8h  nrotat)
  628    format(4x,4hgrid,3x,6hirotat,3x,6hrfreqr,3x,6homegax,
     .          3x,6homegay,3x,6homegaz,4x,5hxorig,
     .          4x,5hyorig,4x,5hzorig)
  632    format(4x,4hgrid,3x,6hdthxmx,3x,6hdthymx,3x,6hdthzmx)
c
      end if
c
c***** input for elastic grid motion  *****
c
c     idefrm  - deforming mesh flag
c             = 0...rigid mesh
c             > 0...deforming mesh 
c               (99 indicates deformation 
c               due to aeroelastic effects)
c
c******************************************
c
      mxmds    = 1
      mxaes    = 1
      mxdefseg = 1
      nsprgit  = 0
c
c     default to no deformation, no aeroelasticty
c
      naesrf = 0
      ndefrm = 0
c
      do iaes=1,maxaes
         do ll=1,5
            aesrfdat(ll,iaes) = 0.
         end do
         do n = 1,nmds
            nm1  = 2*n-1
            nm2  = 2*n
            freq(n,iaes)  = 0.
            gmass(n,iaes) = 0.
            damp(n,iaes)  = 0.
            x0(nm1,iaes)  = 0.
            x0(nm2,iaes)  = 0.
            gf0(nm1,iaes) = 0.
            gf0(nm2,iaes) = 0.
            do ll=1,4
               perturb(n,iaes,ll) = 0.
            end do
         end do
      end do
c
      do igrid=1,ngrid
         nbl = nblg(igrid)
         ncg = ncgg(igrid)
         idefrm(nbl)   = 0
         nsegdfrm(nbl) = 0
         do iseg=1,maxsegdg
            utrnsae(nbl,iseg)  = 0.
            vtrnsae(nbl,iseg)  = 0.
            wtrnsae(nbl,iseg)  = 0.
            omgxae(nbl,iseg)   = 0.
            omgyae(nbl,iseg)   = 0.
            omgzae(nbl,iseg)   = 0.
            xorgae(nbl,iseg)   = 0.
            yorgae(nbl,iseg)   = 0.
            zorgae(nbl,iseg)   = 0.
            xorgae0(nbl,iseg)  = 0.
            yorgae0(nbl,iseg)  = 0.
            zorgae0(nbl,iseg)  = 0.
            icouple(nbl,iseg)  = 0
            thtxae(nbl,iseg)   = 0.
            thtyae(nbl,iseg)   = 0.
            thtzae(nbl,iseg)   = 0.
            rfrqtae(nbl,iseg)  = 0.
            rfrqrae(nbl,iseg)  = 0.
            icsi(nbl,iseg)     = 0
            icsf(nbl,iseg)     = 0
            jcsi(nbl,iseg)     = 0
            jcsf(nbl,iseg)     = 0
            kcsi(nbl,iseg)     = 0
            kcsf(nbl,iseg)     = 0
            idfrmseg(nbl,iseg) = 0
            iaesurf(nbl,iseg)  = 0
         end do
         if (ncg.gt.0) then
            do n=1,ncg
               nbl           = nbl+1
               idefrm(nbl)   = idefrm(nbl-1)
               nsegdfrm(nbl) = nsegdfrm(nbl-1)
               do iseg=1,maxsegdg
                  utrnsae(nbl,iseg)  = utrnsae(nbl-1,iseg)
                  vtrnsae(nbl,iseg)  = vtrnsae(nbl-1,iseg)
                  wtrnsae(nbl,iseg)  = wtrnsae(nbl-1,iseg)
                  omgxae(nbl,iseg)   = omgxae(nbl-1,iseg)
                  omgyae(nbl,iseg)   = omgyae(nbl-1,iseg)
                  omgzae(nbl,iseg)   = omgzae(nbl-1,iseg)
                  xorgae(nbl,iseg)   = xorgae(nbl-1,iseg)
                  yorgae(nbl,iseg)   = yorgae(nbl-1,iseg)
                  zorgae(nbl,iseg)   = zorgae(nbl-1,iseg)
                  xorgae0(nbl,iseg)  = xorgae0(nbl-1,iseg)
                  yorgae0(nbl,iseg)  = yorgae0(nbl-1,iseg)
                  zorgae0(nbl,iseg)  = zorgae0(nbl-1,iseg)
                  icouple(nbl,iseg)  = icouple(nbl-1,iseg)
                  thtxae(nbl,iseg)   = thtxae(nbl-1,iseg)
                  thtyae(nbl,iseg)   = thtyae(nbl-1,iseg)
                  thtzae(nbl,iseg)   = thtzae(nbl-1,iseg)
                  rfrqtae(nbl,iseg)  = rfrqtae(nbl-1,iseg)
                  rfrqrae(nbl,iseg)  = rfrqrae(nbl-1,iseg)
                  icsi(nbl,iseg)     = icsi(nbl-1,iseg)
                  icsf(nbl,iseg)     = icsf(nbl-1,iseg)
                  jcsi(nbl,iseg)     = jcsi(nbl-1,iseg)
                  jcsf(nbl,iseg)     = jcsf(nbl-1,iseg)
                  kcsi(nbl,iseg)     = kcsi(nbl-1,iseg)
                  kcsf(nbl,iseg)     = kcsf(nbl-1,iseg)
                  idfrmseg(nbl,iseg) = idfrmseg(nbl-1,iseg)
                  iaesurf(nbl,iseg)  = iaesurf(nbl-1,iseg)
               end do
            end do
         end if
      end do
c
c     deforming mesh data
c
      if (iunst.gt.1 .or. idef_ss.gt.0) then
         allocate( riskp(maxbl,500), stat=stats )
         allocate( rjskp(maxbl,500), stat=stats )
         allocate( rkskp(maxbl,500), stat=stats )
         allocate( iskipt(maxbl,500), stat=stats )
         allocate( jskipt(maxbl,500), stat=stats )
         allocate( kskipt(maxbl,500), stat=stats )
c
c        set default skip values for block
c
         do nb = 1,maxbl
           do ii = 1,500
             iskip(nb,ii) = 0
             jskip(nb,ii) = 0
             kskip(nb,ii) = 0
           enddo
         enddo
         do ig=1,ngrid
            nbl     = nblg(ig) 
            jdim1   = jdimg(nbl)-1
            kdim1   = kdimg(nbl)-1
            idim1   = idimg(nbl)-1
            ntimes  = 10
c
            jrngmin = jdim1
            krngmin = kdim1
            do iseg=1,nbci0(nbl)
               jrange = abs(ibcinfo(nbl,iseg,3,1)-ibcinfo(nbl,iseg,2,1))
               krange = abs(ibcinfo(nbl,iseg,5,1)-ibcinfo(nbl,iseg,4,1))
               if (jrange.lt.jrngmin) jrngmin = jrange
               if (krange.lt.krngmin) krngmin = krange
               do jexp=1,ntimes
                  power2 = 2**(jexp-1)
                  fact   = jrange/power2
                  if (real(int(fact)).eq.real(fact)) jskip_i0 = power2
               end do
               do kexp=1,ntimes
                  power2 = 2**(kexp-1)
                  fact   = krange/power2
                  if (real(int(fact)).eq.real(fact)) kskip_i0 = power2
               end do
            end do
            do iseg=1,nbcidim(nbl)
               jrange = abs(ibcinfo(nbl,iseg,3,2)-ibcinfo(nbl,iseg,2,2))
               krange = abs(ibcinfo(nbl,iseg,5,2)-ibcinfo(nbl,iseg,4,2))
               if (jrange.lt.jrngmin) jrngmin = jrange
               if (krange.lt.krngmin) krngmin = krange
               do jexp=1,ntimes
                  power2 = 2**(jexp-1)
                  fact   = jrange/power2
                  if (real(int(fact)).eq.real(fact)) jskip_id = power2
               end do
               do kexp=1,ntimes
                  power2 = 2**(kexp-1)
                  fact   = krange/power2
                  if (real(int(fact)).eq.real(fact)) kskip_id = power2
               end do
            end do
            if (jdim1/jrngmin*jrngmin .eq. jdim1) then
               jskip_i = jrngmin
            else
               jskip_i = min(jskip_i0,jskip_id)
            end if
            if (kdim1/krngmin*krngmin .eq. kdim1) then
               kskip_i = krngmin
            else
               kskip_i = min(kskip_i0,kskip_id)
            end if
c
            irngmin = idim1
            krngmin = kdim1
            do jseg=1,nbcj0(nbl)
               irange = abs(jbcinfo(nbl,jseg,3,1)-jbcinfo(nbl,jseg,2,1))
               krange = abs(jbcinfo(nbl,jseg,5,1)-jbcinfo(nbl,jseg,4,1))
               if (irange.lt.irngmin) irngmin = irange
               if (krange.lt.krngmin) krngmin = krange
               do iexp=1,ntimes
                  power2 = 2**(iexp-1)
                  fact   = irange/power2
                  if (real(int(fact)).eq.real(fact)) iskip_j0 = power2
               end do
               do kexp=1,ntimes
                  power2 = 2**(kexp-1)
                  fact   = krange/power2
                  if (real(int(fact)).eq.real(fact)) kskip_j0 = power2
               end do
            end do
            do jseg=1,nbcjdim(nbl)
               irange = abs(jbcinfo(nbl,jseg,3,2)-jbcinfo(nbl,jseg,2,2))
               krange = abs(jbcinfo(nbl,jseg,5,2)-jbcinfo(nbl,jseg,4,2))
               if (irange.lt.irngmin) irngmin = irange
               if (krange.lt.krngmin) krngmin = krange
               do iexp=1,ntimes
                  power2 = 2**(iexp-1)
                  fact   = irange/power2
                  if (real(int(fact)).eq.real(fact)) iskip_jd = power2
               end do
               do kexp=1,ntimes
                  power2 = 2**(kexp-1)
                  fact   = krange/power2
                  if (real(int(fact)).eq.real(fact)) kskip_jd = power2
               end do
            end do
            if (idim1/irngmin*irngmin .eq. idim1) then
               iskip_j = irngmin
            else
               iskip_j = min(iskip_j0,iskip_jd)
            end if
            if (kdim1/krngmin*krngmin .eq. kdim1) then
               kskip_j = krngmin
            else
               kskip_j = min(kskip_j0,kskip_jd)
            end if
c
            irngmin = idim1
            jrngmin = jdim1
            do kseg=1,nbck0(nbl)
               irange = abs(kbcinfo(nbl,kseg,3,1)-kbcinfo(nbl,kseg,2,1))
               jrange = abs(kbcinfo(nbl,kseg,5,1)-kbcinfo(nbl,kseg,4,1))
               if (irange.lt.irngmin) irngmin = irange
               if (jrange.lt.jrngmin) jrngmin = jrange
               do iexp=1,ntimes
                  power2 = 2**(iexp-1)
                  fact   = irange/power2
                  if (real(int(fact)).eq.real(fact)) iskip_k0 = power2
               end do
               do jexp=1,ntimes
                  power2 = 2**(jexp-1)
                  fact   = jrange/power2
                  if (real(int(fact)).eq.real(fact)) jskip_k0 = power2
               end do
            end do
            do kseg=1,nbckdim(nbl)
               irange = abs(kbcinfo(nbl,kseg,3,2)-kbcinfo(nbl,kseg,2,2))
               jrange = abs(kbcinfo(nbl,kseg,5,2)-kbcinfo(nbl,kseg,4,2))
               if (irange.lt.irngmin) irngmin = irange
               if (jrange.lt.jrngmin) jrngmin = jrange
               do iexp=1,ntimes
                  power2 = 2**(iexp-1)
                  fact   = irange/power2
                  if (real(int(fact)).eq.real(fact)) iskip_kd = power2
               end do
               do jexp=1,ntimes
                  power2 = 2**(jexp-1)
                  fact   = jrange/power2
                  if (real(int(fact)).eq.real(fact)) jskip_kd = power2
               end do
            end do
            if (idim1/irngmin*irngmin .eq. idim1) then
               iskip_k = irngmin
            else
               iskip_k = min(iskip_k0,iskip_kd)
            end if
            if (jdim1/jrngmin*jrngmin .eq. jdim1) then
               jskip_k = jrngmin
            else
               jskip_k = min(jskip_k0,jskip_kd)
            end if
c
            iskip(nbl,1) = min(iskip_j,iskip_k)
            jskip(nbl,1) = min(jskip_i,jskip_k)
            kskip(nbl,1) = min(kskip_i,kskip_j)
c
         end do 
c
c        account for 1-1 interfaces
c
         if (abs(nbli).gt.0) then
           do it = 1,ngrid
            nskpdif = 0
            do n=1,abs(nbli)
               igr1 = nblk(1,n)
               igr2 = nblk(2,n)
               nbl1 = nblg(igr1)
               nbl2 = nblg(igr2)
               iskipt(nbl1,1) = iskip(nbl1,1)
               iskipt(nbl2,1) = iskip(nbl2,1)
               jskipt(nbl1,1) = jskip(nbl1,1)
               jskipt(nbl2,1) = jskip(nbl2,1)
               kskipt(nbl1,1) = kskip(nbl1,1)
               kskipt(nbl2,1) = kskip(nbl2,1)
               if (isva(1,1,n) .eq. 1) then
                  if (isva(2,1,n) .eq. 1) then
                    call blocking_skip(iskip(nbl1,1) ,iskip(nbl2,1),
     .                                 iskipt(nbl1,1),iskipt(nbl2,1))
                  else if (isva(2,1,n) .eq. 2) then
                    call blocking_skip(iskip(nbl1,1) ,jskip(nbl2,1),
     .                                 iskipt(nbl1,1),jskipt(nbl2,1))
                  else if (isva(2,1,n) .eq. 3) then
                    call blocking_skip(iskip(nbl1,1) ,kskip(nbl2,1),
     .                                 iskipt(nbl1,1),kskipt(nbl2,1))
                  end if
               else if (isva(1,1,n) .eq. 2) then
                  if (isva(2,1,n) .eq. 1) then
                    call blocking_skip(jskip(nbl1,1) ,iskip(nbl2,1),
     .                                 jskipt(nbl1,1),iskipt(nbl2,1))
                  else if (isva(2,1,n) .eq. 2) then
                    call blocking_skip(jskip(nbl1,1) ,jskip(nbl2,1),
     .                                 jskipt(nbl1,1),jskipt(nbl2,1))
                  else if (isva(2,1,n) .eq. 3) then
                    call blocking_skip(jskip(nbl1,1) ,kskip(nbl2,1),
     .                                 jskipt(nbl1,1),kskipt(nbl2,1))
                  end if
               else if (isva(1,1,n) .eq. 3) then
                  if (isva(2,1,n) .eq. 1) then
                    call blocking_skip(kskip(nbl1,1) ,iskip(nbl2,1),
     .                                 kskipt(nbl1,1),iskipt(nbl2,1))
                  else if (isva(2,1,n) .eq. 2) then
                    call blocking_skip(kskip(nbl1,1) ,jskip(nbl2,1),
     .                                 kskipt(nbl1,1),jskipt(nbl2,1))
                  else if (isva(2,1,n) .eq. 3) then
                    call blocking_skip(kskip(nbl1,1) ,kskip(nbl2,1),
     .                                 kskipt(nbl1,1),kskipt(nbl2,1))
                  end if
               end if
c
               if (isva(1,2,n) .eq. 1) then
                  if (isva(2,2,n) .eq. 1) then
                    call blocking_skip(iskip(nbl1,1) ,iskip(nbl2,1),
     .                                 iskipt(nbl1,1),iskipt(nbl2,1))
                  else if (isva(2,2,n) .eq. 2) then
                    call blocking_skip(iskip(nbl1,1) ,jskip(nbl2,1),
     .                                 iskipt(nbl1,1),jskipt(nbl2,1))
                  else if (isva(2,2,n) .eq. 3) then
                    call blocking_skip(iskip(nbl1,1) ,kskip(nbl2,1),
     .                                 iskipt(nbl1,1),kskipt(nbl2,1))
                  end if
               else if (isva(1,2,n) .eq. 2) then
                  if (isva(2,2,n) .eq. 1) then
                    call blocking_skip(jskip(nbl1,1) ,iskip(nbl2,1),
     .                                 jskipt(nbl1,1),iskipt(nbl2,1))
                  else if (isva(2,2,n) .eq. 2) then
                    call blocking_skip(jskip(nbl1,1) ,jskip(nbl2,1),
     .                                 jskipt(nbl1,1),jskipt(nbl2,1))
                  else if (isva(2,2,n) .eq. 3) then
                    call blocking_skip(jskip(nbl1,1) ,kskip(nbl2,1),
     .                                 jskipt(nbl1,1),kskipt(nbl2,1))
                  end if
               else if (isva(1,2,n) .eq. 3) then
                  if (isva(2,2,n) .eq. 1) then
                    call blocking_skip(kskip(nbl1,1) ,iskip(nbl2,1),
     .                                 kskipt(nbl1,1),iskipt(nbl2,1))
                  else if (isva(2,2,n) .eq. 2) then
                    call blocking_skip(kskip(nbl1,1) ,jskip(nbl2,1),
     .                                 kskipt(nbl1,1),jskipt(nbl2,1))
                  else if (isva(2,2,n) .eq. 3) then
                    call blocking_skip(kskip(nbl1,1) ,kskip(nbl2,1),
     .                                 kskipt(nbl1,1),kskipt(nbl2,1))
                  end if
               end if
               nskpdif = nskpdif + abs(iskip(nbl1,1) - iskipt(nbl1,1))
               nskpdif = nskpdif + abs(iskip(nbl2,1) - iskipt(nbl2,1))
               nskpdif = nskpdif + abs(jskip(nbl1,1) - jskipt(nbl1,1))
               nskpdif = nskpdif + abs(jskip(nbl2,1) - jskipt(nbl2,1))
               nskpdif = nskpdif + abs(kskip(nbl1,1) - kskipt(nbl1,1))
               nskpdif = nskpdif + abs(kskip(nbl2,1) - kskipt(nbl2,1))
               iskip(nbl1,1) = iskipt(nbl1,1)
               iskip(nbl2,1) = iskipt(nbl2,1)
               jskip(nbl1,1) = jskipt(nbl1,1)
               jskip(nbl2,1) = jskipt(nbl2,1)
               kskip(nbl1,1) = kskipt(nbl1,1)
               kskip(nbl2,1) = kskipt(nbl2,1)
            end do
            if(nskpdif.eq.0) goto 850
           end do
         end if
850      continue
c
         read(iunit5,*)
         read(iunit5,*)
         read(iunit5,*) ndefrm
         write(iunit11,'('' moving grid data - deforming surface'',
     .                   '' (forced motion):'')')
         write(iunit11,'(''  ndefrm'')')
         if (ndefrm.ge.0) then
            write(iunit11,'(i8)') ndefrm
         end if
         if (ndefrm.gt.0) then
            read(iunit5,*)
            read(iunit5,*) realval(1)
            lrefdef = realval(1)
            read(iunit5,*)
            write(iunit11,'(''    lref'')')
            write(iunit11,'(f8.4)') real(lrefdef)
            write(iunit11,'(''    grid   idefrm    rfreq u/omegax'',
     .                      '' v/omegay w/omegaz    xorig'',
     .                      ''    yorig    zorig'')')
            do ndef=1,ndefrm
               read(iunit5,*) igr,idef,(realval(i),i=1,7)
               rfr   = realval(1)
               omgx  = realval(2)
               omgy  = realval(3)
               omgz  = realval(4)
               origx = realval(5)
               origy = realval(6)
               origz = realval(7)
               nbl                = nblg(igr)
               ncg                = ncgg(igr)
               idefrm(nbl)        = 1
               nsegdfrm(nbl)      = nsegdfrm(nbl) + 1
               iseg               = nsegdfrm(nbl)
               idfrmseg(nbl,iseg) = idef
               if (idef.eq.2 .or. idef.eq.3) then
c                 deformation via rotation
                  rfrqrae(nbl,iseg) = 2.*pi*rfr/lrefdef
                  omgxae(nbl,iseg)  = omgx/radtodeg
                  omgyae(nbl,iseg)  = omgy/radtodeg
                  omgzae(nbl,iseg)  = omgz/radtodeg
                  xorgae(nbl,iseg)  = origx
                  yorgae(nbl,iseg)  = origy
                  zorgae(nbl,iseg)  = origz
                  xorgae0(nbl,iseg) = origx
                  yorgae0(nbl,iseg) = origy
                  zorgae0(nbl,iseg) = origz
                  write(iunit11,'(i8,i9,7f9.4)')igr,idef,real(rfr),
     .                  real(omgx),real(omgy),real(omgz),real(origx),
     .                  real(origy),real(origz)
               else if (idef.eq.1) then
c                 deformation via translation
                  rfrqtae(nbl,iseg) = 2.*pi*rfr/lrefdef
                  utrnsae(nbl,iseg) = omgx
                  vtrnsae(nbl,iseg) = omgy
                  wtrnsae(nbl,iseg) = omgz
                  xorgae(nbl,iseg)  = origx
                  yorgae(nbl,iseg)  = origy
                  zorgae(nbl,iseg)  = origz
                  xorgae0(nbl,iseg) = origx
                  yorgae0(nbl,iseg) = origy
                  zorgae0(nbl,iseg) = origz
                  write(iunit11,'(i8,i9,7f9.4)')igr,idef,real(rfr),
     .                  real(omgx),real(omgy),real(omgz),real(origx),
     .                  real(origy),real(origz)
               else
c                 no deformation...output (and use) default values
                  write(iunit11,'(i8,i9,7f9.4)')igr,idfrmseg(nbl,iseg),
     .            real(rfrqrae(nbl,iseg)),real(omgxae(nbl,iseg)),
     .            real(omgyae(nbl,iseg)),real(omgzae(nbl,iseg)),
     .            real(xorgae(nbl,iseg)),real(yorgae(nbl,iseg)),
     .            real(zorgae(nbl,iseg))
               end if
               if (ncg.gt.0) then
                  do n=1,ncg
                     nbl = nbl + 1
                     idefrm(nbl)        = idefrm(nbl-1)
                     nsegdfrm(nbl)      = nsegdfrm(nbl-1)
                     idfrmseg(nbl,iseg) = idfrmseg(nbl-1,iseg)
                     rfrqrae(nbl,iseg)  = rfrqrae(nbl-1,iseg)
                     omgxae(nbl,iseg)   = omgxae(nbl-1,iseg)
                     omgyae(nbl,iseg)   = omgyae(nbl-1,iseg)
                     omgzae(nbl,iseg)   = omgzae(nbl-1,iseg)
                     rfrqtae(nbl,iseg)  = rfrqtae(nbl-1,iseg)
                     utrnsae(nbl,iseg)  = utrnsae(nbl-1,iseg)
                     vtrnsae(nbl,iseg)  = vtrnsae(nbl-1,iseg)
                     wtrnsae(nbl,iseg)  = wtrnsae(nbl-1,iseg)
                     xorgae(nbl,iseg)   = xorgae(nbl-1,iseg)
                     yorgae(nbl,iseg)   = yorgae(nbl-1,iseg)
                     zorgae(nbl,iseg)   = zorgae(nbl-1,iseg)
                     xorgae0(nbl,iseg)  = xorgae0(nbl-1,iseg)
                     yorgae0(nbl,iseg)  = yorgae0(nbl-1,iseg)
                     zorgae0(nbl,iseg)  = zorgae0(nbl-1,iseg)
                  end do
               end if
            end do
            read(iunit5,*)
            write(iunit11,'(''    grid     icsi     icsf     jcsi'',
     .                      ''     jcsf     kcsi     kcsf'')')
c
c           reset nsegdfrm so it can be used again for setting the 
c           index ranges. NOTE: this assumes that the motion data
c           set above and the index range data set below are in the
c           same order on input
c
            do nbl=1,nblock
               nsegdfrm(nbl) = 0
            end do
c
            do ndef=1,ndefrm
               read(iunit5,*) igr,is,ie,js,je,ks,ke
               write(iunit11,'(i8,6i9)') igr,is,ie,js,je,ks,ke
               nbl = nblg(igr)
               ncg = ncgg(igr)
               nsegdfrm(nbl) = nsegdfrm(nbl) + 1
               iseg = nsegdfrm(nbl)
               icsi(nbl,iseg) = is
               icsf(nbl,iseg) = ie
               jcsi(nbl,iseg) = js
               jcsf(nbl,iseg) = je
               kcsi(nbl,iseg) = ks
               kcsf(nbl,iseg) = ke
               if (ncg.gt.0) then
                  do n=1,ncg
                     nbl = nbl + 1
                     nsegdfrm(nbl)  = nsegdfrm(nbl-1)
                     icsi(nbl,iseg) = icsi(nbl-1,iseg)/2 + 1
                     icsf(nbl,iseg) = icsf(nbl-1,iseg)/2 + 1
                     jcsi(nbl,iseg) = jcsi(nbl-1,iseg)/2 + 1
                     jcsf(nbl,iseg) = jcsf(nbl-1,iseg)/2 + 1
                     kcsi(nbl,iseg) = kcsi(nbl-1,iseg)/2 + 1
                     kcsf(nbl,iseg) = kcsf(nbl-1,iseg)/2 + 1
                  end do
               end if
            end do
c
         else if (ndefrm.lt.0) then
c
c           set segment counter and index ranges to correspond to all 
c           solid surfaces.
c
            ndefrm0 = 0
            call setseg(maxgr,maxbl,maxseg,nblg,ncgg,idimg,jdimg,
     .                  kdimg,ibcinfo,jbcinfo,kbcinfo,nbci0,
     .                  nbcidim,nbcj0,nbcjdim,nbck0,nbckdim,
     .                  ndefrm0,idefrm,nsegdfrm,icsi,icsf,
     .                  jcsi,jcsf,kcsi,kcsf,maxsegdg,ngrid)
c
            write(iunit11,'(i8)') ndefrm0
            read(iunit5,*)
            read(iunit5,*) realval(1)
            lrefdef = realval(1)
            read(iunit5,*)
            write(iunit11,'(''    lref'')')
            write(iunit11,'(f8.4)') real(lrefdef)
            write(iunit11,'(''    grid   idefrm    rfreq u/omegax'',
     .                      '' v/omegay w/omegaz    xorig'',
     .                      ''    yorig    zorig'')')
            do ndef=1,abs(ndefrm)
               read(iunit5,*) igr,idef,(realval(i),i=1,7)
               rfr   = realval(1)
               omgx  = realval(2)
               omgy  = realval(3)
               omgz  = realval(4)
               origx = realval(5)
               origy = realval(6)
               origz = realval(7)
            end do
            do ig=1,ngrid
               nbl = nblg(ig)
               ncg = ncgg(ig)
               do iseg=1,nsegdfrm(nbl)
                  idfrmseg(nbl,iseg) = idef           
                  if (idef.eq.2 .or. idef.eq.3) then
c                    deformation via rotation
                     rfrqrae(nbl,iseg) = 2.*pi*rfr/lrefdef
                     omgxae(nbl,iseg)  = omgx/radtodeg
                     omgyae(nbl,iseg)  = omgy/radtodeg
                     omgzae(nbl,iseg)  = omgz/radtodeg
                     xorgae(nbl,iseg)  = origx
                     yorgae(nbl,iseg)  = origy
                     zorgae(nbl,iseg)  = origz
                     xorgae0(nbl,iseg) = origx
                     yorgae0(nbl,iseg) = origy
                     zorgae0(nbl,iseg) = origz
                     write(iunit11,'(i8,i9,7f9.4)')ig,idef,real(rfr),
     .                     real(omgx),real(omgy),real(omgz),real(origx),
     .                     real(origy),real(origz)
                  else if (idef.eq.1) then
c                    deformation via translation
                     rfrqtae(nbl,iseg) = 2.*pi*rfr/lrefdef
                     utrnsae(nbl,iseg) = omgx
                     vtrnsae(nbl,iseg) = omgy
                     wtrnsae(nbl,iseg) = omgz
                     xorgae(nbl,iseg)  = origx
                     yorgae(nbl,iseg)  = origy
                     zorgae(nbl,iseg)  = origz
                     xorgae0(nbl,iseg) = origx
                     yorgae0(nbl,iseg) = origy
                     zorgae0(nbl,iseg) = origz
                     write(iunit11,'(i8,i9,7f9.4)')ig,idef,real(rfr),
     .                     real(omgx),real(omgy),real(omgz),real(origx),
     .                     real(origy),real(origz)
                  else
c                    no deformation...output (and use) default values
                     write(iunit11,'(i8,i9,7f9.4)')
     .               ig,idfrmseg(nbl,iseg),
     .               real(rfrqrae(nbl,iseg)),real(omgxae(nbl,iseg)),
     .               real(omgyae(nbl,iseg)),real(omgzae(nbl,iseg)),
     .               real(xorgae(nbl,iseg)),real(yorgae(nbl,iseg)),
     .               real(zorgae(nbl,iseg))
                  end if
               end do
               if (ncg.gt.0) then
                  do n=1,ncg
                     nbl = nbl + 1
                     idefrm(nbl)        = idefrm(nbl-1)
                     nsegdfrm(nbl)      = nsegdfrm(nbl-1)
                     do iseg=1,nsegdfrm(nbl)
                        idfrmseg(nbl,iseg) = idfrmseg(nbl-1,iseg)
                        rfrqrae(nbl,iseg)  = rfrqrae(nbl-1,iseg)
                        omgxae(nbl,iseg)   = omgxae(nbl-1,iseg)
                        omgyae(nbl,iseg)   = omgyae(nbl-1,iseg)
                        omgzae(nbl,iseg)   = omgzae(nbl-1,iseg)
                        rfrqtae(nbl,iseg)  = rfrqtae(nbl-1,iseg)
                        utrnsae(nbl,iseg)  = utrnsae(nbl-1,iseg)
                        vtrnsae(nbl,iseg)  = vtrnsae(nbl-1,iseg)
                        wtrnsae(nbl,iseg)  = wtrnsae(nbl-1,iseg)
                        xorgae(nbl,iseg)   = xorgae(nbl-1,iseg)
                        yorgae(nbl,iseg)   = yorgae(nbl-1,iseg)
                        zorgae(nbl,iseg)   = zorgae(nbl-1,iseg)
                        xorgae0(nbl,iseg)  = xorgae0(nbl-1,iseg)
                        yorgae0(nbl,iseg)  = yorgae0(nbl-1,iseg)
                        zorgae0(nbl,iseg)  = zorgae0(nbl-1,iseg)
                     end do
                  end do
               end if
            end do
c
            read(iunit5,*)
            write(iunit11,'(''    grid     icsi     icsf     jcsi'',
     .                      ''     jcsf     kcsi     kcsf'')')
            do ndef=1,abs(ndefrm)
               read(iunit5,*) igr,is,ie,js,je,ks,ke
            end do
            do ig=1,ngrid
               nbl = nblg(ig)
               do iseg=1,nsegdfrm(nbl)
                  write(iunit11,'(i8,6i9)') ig,icsi(nbl,iseg),
     .            icsf(nbl,iseg),jcsi(nbl,iseg),jcsf(nbl,iseg),
     .            kcsi(nbl,iseg),kcsf(nbl,iseg)
               end do
            end do
c
            ndefrm = ndefrm0
c
         else
            read(iunit5,*)
            write(iunit11,'(''    lref'')')
            read(iunit5,*)
            write(iunit11,'(''    grid  idefrm    rfreq u/omegax'',
     .                      '' v/omegay w/omegaz    xorig'',
     .                      ''    yorig    zorig'')')
            read(iunit5,*)
            write(iunit11,'(''    grid    icsi     icsf     jcsi'',
     .                      ''     jcsf     kcsi     kcsf'')')
         end if
c
c        aeroelastic mesh data
c
         read(iunit5,*)
         read(iunit5,*)
         read(iunit5,*) naesrf
         write(iunit11,'('' moving grid data - aeroelastic surface'',
     .                   '' (aeroelastic motion):'')')
         write(iunit11,'(''  naesrf'')')
         write(iunit11,'(i8)') naesrf
         if (naesrf.gt.0) then
            if (naesrf.gt.maxaes) then
               write(iunit11,'('' input error: naesrf > maxaes'')')
               call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
            end if
            do naes=1,naesrf
               read(iunit5,*)
               write(iunit11,'(''  iaesrf    ngrid    grefl     uinf'',
     .                         ''     qinf   nmodes iskyhook'')')
               read(iunit5,*) iaes,ngd,(realval(i),i=1,3),nmodes,iskyhk
               grefl = realval(1)
               uinf  = realval(2)
               qinf  = realval(3)
               aesrfdat(1,iaes) = iskyhk
               aesrfdat(2,iaes) = grefl
               aesrfdat(3,iaes) = uinf
               aesrfdat(4,iaes) = qinf
               aesrfdat(5,iaes) = nmodes
               write(iunit11,'(i8,i9,3f9.4,2i9)') iaes,ngd,real(grefl),
     .                             real(uinf),real(qinf),nmodes,iskyhk
               if (real(grefl).le.0.) then
                  write(iunit11,'('' input error: invalid grefl'')')
                  call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
               end if
               if (nmodes.gt.nmds) then
                  write(iunit11,'('' input error: number of modes'',
     .                            '' exceeds nmds'')')
                  call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
               end if
               read(iunit5,*)
               write(iunit11,'(''    freq    gmass     damp x0(2n-1)'',
     .                         ''   x0(2n)  gf0(2n)'')')
               do nm = 1,nmodes
                  nm1  = 2*nm-1
                  nm2  = 2*nm
                  read(iunit5,*) (realval(i),i=1,6)
                  freq1  = realval(1)
                  gmass1 = realval(2)
                  damp1  = realval(3)
                  x01    = realval(4)
                  x02    = realval(5)
                  gf02   = realval(6)
                  freq(nm,iaes)  = freq1
                  gmass(nm,iaes) = gmass1
                  damp(nm,iaes)  = damp1
                  x0(nm1,iaes)   = x01
                  x0(nm2,iaes)   = x02
                  gf0(nm1,iaes)  = gf02
                  gf0(nm2,iaes)  = gf02
                  write(iunit11,'(f8.4,5f9.4)') real(freq1),
     .                  real(gmass1),real(damp1),real(x01),
     .                  real(x02),real(gf02)
                  if (real(freq1).le.0.) then
                     write(iunit11,'('' input error: invalid'',
     .                               '' frequency'')')
                     call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
                  end if
                  if (real(damp1).ge.1.) then
                     write(iunit11,'('' input error: invalid'',
     .                               ''  damping'')')
                     call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
                  end if
               end do
               read(iunit5,*)
               write(iunit11,'(''  moddfl      amp     freq'',
     .                         ''       t0'')')
               do nm = 1,nmodes
                  read(iunit5,*)  moddfl,(realval(i),i=1,3)
                  amp   = realval(1)
                  freqp = realval(2)
                  t0    = realval(3)
                  perturb(nm,iaes,1) = float(moddfl)
                  perturb(nm,iaes,2) = amp
                  perturb(nm,iaes,3) = freqp
                  perturb(nm,iaes,4) = t0
                  write(iunit11,'(i8,4f9.4)') moddfl,real(amp),
     .                                        real(freqp),real(t0)
               end do
               read(iunit5,*)
               write(iunit11,'(''    grid     iaei     iaef     jaei'',
     .                         ''     jaef     kaei     kaef'')')
               if (ngd.gt.0) then
                  do ng = 1,ngd
                     read(iunit5,*) igr,is,ie,js,je,ks,ke
                     write(iunit11,'(i8,6i9)') igr,is,ie,js,je,ks,ke
                     nbl                = nblg(igr)
                     ncg                = ncgg(igr)
                     nsegdfrm(nbl)      = nsegdfrm(nbl)+1
                     iseg               = nsegdfrm(nbl)
                     icsi(nbl,iseg)     = is
                     icsf(nbl,iseg)     = ie
                     jcsi(nbl,iseg)     = js
                     jcsf(nbl,iseg)     = je
                     kcsi(nbl,iseg)     = ks
                     kcsf(nbl,iseg)     = ke
                     idfrmseg(nbl,iseg) = 99
                     iaesurf(nbl,iseg)  = iaes
                     idefrm(nbl)        = 1
                     if (ncg.gt.0) then
                        do n=1,ncg
                           nbl = nbl + 1
                           idefrm(nbl)        = idefrm(nbl-1)
                           nsegdfrm(nbl)      = nsegdfrm(nbl-1)
                           idfrmseg(nbl,iseg) = idfrmseg(nbl-1,iseg)
                           iaesurf(nbl,iseg)  = iaesurf(nbl-1,iseg)
                           icsi(nbl,iseg)     = icsi(nbl-1,iseg)/2 + 1
                           icsf(nbl,iseg)     = icsf(nbl-1,iseg)/2 + 1
                           jcsi(nbl,iseg)     = jcsi(nbl-1,iseg)/2 + 1
                           jcsf(nbl,iseg)     = jcsf(nbl-1,iseg)/2 + 1
                           kcsi(nbl,iseg)     = kcsi(nbl-1,iseg)/2 + 1
                           kcsf(nbl,iseg)     = kcsf(nbl-1,iseg)/2 + 1
                        end do
                     end if
                  end do
               else if (ngd.lt.0) then
c
c                 set segment counter and index ranges to correspond to
c                 all solid surfaces.
c
                  ndefrm0 = 0
                  call setseg(maxgr,maxbl,maxseg,nblg,ncgg,idimg,jdimg,
     .                        kdimg,ibcinfo,jbcinfo,kbcinfo,nbci0,
     .                        nbcidim,nbcj0,nbcjdim,nbck0,nbckdim,
     .                        ndefrm0,idefrm,nsegdfrm,icsi,icsf,
     .                        jcsi,jcsf,kcsi,kcsf,maxsegdg,ngrid)
c
                  do ng = 1,abs(ngd)
                     read(iunit5,*) igr,is,ie,js,je,ks,ke
                  end do
                  do ig=1,ngrid
                     nbl = nblg(ig)
                     ncg = ncgg(ig)
                     do iseg=1,nsegdfrm(nbl)
                        idfrmseg(nbl,iseg) = 99
                        iaesurf(nbl,iseg)  = iaes
                        write(iunit11,'(i8,6i9)') ig,icsi(nbl,iseg),
     .                  icsf(nbl,iseg), jcsi(nbl,iseg),jcsf(nbl,iseg),
     .                  kcsi(nbl,iseg),kcsf(nbl,iseg)
                     end do
                     if (ncg.gt.0) then
                         do n=1,ncg
                            nbl = nbl + 1
                            do iseg=1,nsegdfrm(nbl)
                               idfrmseg(nbl,iseg) = idfrmseg(nbl-1,iseg)
                               iaesurf(nbl,iseg)  = iaesurf(nbl-1,iseg)
                            end do
                         end do
                     end if
                  end do
               end if
            end do
c
            mxmds = max(mxmds,nmodes)
c
         else
            read(iunit5,*)
            read(iunit5,*)
            read(iunit5,*)
            read(iunit5,*)
            write(iunit11,'(''  iaesrf    ngrid    grefl     uinf'',
     .                      ''     qinf   nmodes iskyhook'')')
            write(iunit11,'(''    freq    gmass     damp x0(2n-1)'',
     .                      ''   x0(2n)  gf0(2n)'')')
            write(iunit11,'(''  moddfl      amp     freq       t0'')')
            write(iunit11,'(''    grid     iaei     iaef     jaei'',
     .                      ''     jaef     kaei     kaef'')')
         end if
c
         mxaes  = max(mxaes,naesrf)
c
c        deforming mesh skip data - used to determine slave points
c        for offbody/multiblock mesh movement
c
c        Set default skip type so that standard skip values are read:
c
         isktyp = 1
c
c        Set counter for number of control point updates to ensure
c        1-1 blocking interface continuity:
c
         itcpadd = 0
c
c        If isktyp = 2 then an array of points in the i,j,k directions
c        are read in.  These points are the control points rather than
c        those computed using standard skip values
c
         read(iunit5,*)
         read(iunit5,*)
         read(iunit5,*) nskip,isktyp,realval(1),realval(2),
     .                  realval(3),realval(4),nsprgit 
         beta1 = realval(1)
         alpha1= realval(2)
         beta2 = realval(3)
         alpha2= realval(4)
         write(iunit11,'('' moving grid data - data for field/'',
     .                   ''multiblock mesh movement'')')
         write(iunit11,
     .     '(''   nskip   isktyp    beta1   alpha1    beta2   alpha2'',
     .                   '' nsprngit'')')
         if(abs(isktyp).lt.1.or.abs(isktyp).gt.2) then
           write(iunit11,'(i8,i9,4f9.6,i9)') nskip,isktyp,beta1,alpha1,
     .          beta2,alpha2,nsprgit
           write(iunit11,'('' input error: Invalid isktyp value'')')
           call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
         if(abs(isktyp).eq.1) then
           if (nskip.gt.0) then
              if (nskip.gt.ngrid) then
               write(iunit11,'(i8,i9,4f9.6,i9)') nskip,isktyp,beta1,
     .              alpha1,beta2,alpha2,nsprgit
               write(iunit11,'('' input error: nskip > ngrid'')')
               call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
              end if
              write(iunit11,'(i8,i9,4f9.6,i9)') nskip,isktyp,
     .             real(beta1),real(alpha1),real(beta2),real(alpha2),
     .             nsprgit
              if (real(beta1).lt.0.2) then
                 write(iunit11,'(/,'' WARNING: beta this small may '',
     .        ''severely compromise fidelity of surface geometry'',/)')
              end if
 
              read(iunit5,*)
              do ng = 1,nskip
                 read(iunit5,*) igr,iskp,jskp,kskp
                 nbl = nblg(igr)
                 ncg = ncgg(igr)
c                zero is shortcut for max skip value
                 if (iskp.eq.0) then
                    iskp = idimg(nbl)-1
                 end if
                 if (jskp.eq.0) then
                    jskp = jdimg(nbl)-1
                 end if
                 if (kskp.eq.0) then
                    kskp = kdimg(nbl)-1
                 end if
                 if (i2d.ne.0)  iskp = idimg(nbl)-1
c                skip value must divide evenly into dim-1 value
                 istop = 0
                 if ((idimg(nbl)-1)/iskp*iskp.ne.(idimg(nbl)-1)) then
                    istop = 1
                    write(iunit11,'('' stopping...iskip must divide'',
     .                  '' evenly into idim-1 for grid'',i4)') igr
                 end if
                 if ((jdimg(nbl)-1)/jskp*jskp.ne.(jdimg(nbl)-1)) then
                    istop = 1
                    write(iunit11,'('' stopping...jskip must divide'',
     .                    '' evenly into jdim-1 for grid'',i4)') igr
                 end if
                 if ((kdimg(nbl)-1)/kskp*kskp.ne.(kdimg(nbl)-1)) then
                    istop = 1
                    write(iunit11,'('' stopping...kskip must divide'',
     .                  '' evenly into kdim-1 for grid'',i4)') igr
                 end if
                 if (istop .eq. 1) then
                    call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
                 end if
                 iskip(nbl,1) = iskp
                 jskip(nbl,1) = jskp
                 kskip(nbl,1) = kskp
              end do
           else if (nskip.lt.0) then
               write(iunit11,'(i8,i9,4f9.6,i9)') ngrid,isktyp,
     .              real(beta1),real(alpha1),real(beta2),real(alpha2),
     .              nsprgit
              if (real(beta1).lt.0.2) then
                 write(iunit11,'(/,'' WARNING: beta this small may '',
     .        ''severely compromise fidelity of surface geometry'',/)')
              end if
              read(iunit5,*)
              do ng = 1,abs(nskip)
                 read(iunit5,*) igr,iskp0,jskp0,kskp0
              end do
              do igr=1,ngrid
                 nbl = nblg(igr)
                 ncg = ncgg(igr)
                 iskp = iskp0
                 jskp = jskp0
                 kskp = kskp0
c                zero is shortcut for max skip value
                 if (iskp0.eq.0) iskp = idimg(nbl)-1
                 if (jskp0.eq.0) jskp = jdimg(nbl)-1
                 if (kskp0.eq.0) kskp = kdimg(nbl)-1
c                skip value must divide evenly into dim-1 value
                 istop = 0
                 if ((idimg(nbl)-1)/iskp*iskp.ne.(idimg(nbl)-1)) then
                    istop = 1
                    write(iunit11,'('' stopping...iskip must divide'',
     .                    '' evenly into idim-1 for grid'',i4)') igr
                 end if
                 if ((jdimg(nbl)-1)/jskp*jskp.ne.(jdimg(nbl)-1)) then
                    istop = 1
                    write(iunit11,'('' stopping...jskip must divide'',
     .                    '' evenly into jdim-1 for grid'',i4)') igr
                 end if
                 if ((kdimg(nbl)-1)/kskp*kskp.ne.(kdimg(nbl)-1)) then
                    istop = 1
                    write(iunit11,'('' stopping...kskip must divide'',
     .                    '' evenly into kdim-1 for grid'',i4)') igr
                 end if
                 if (istop .eq. 1) then
                    call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
                 end if
                 iskip(nbl,1) = iskp
                 jskip(nbl,1) = jskp
                 kskip(nbl,1) = kskp
              end do
           else
               write(iunit11,'(i8,i9,4f9.6,i9)') ngrid,isktyp,
     .              real(beta1),real(alpha1),real(beta2),real(alpha2),
     .              nsprgit
              if (real(beta1).lt.0.2) then
                 write(iunit11,'(/,'' WARNING: beta this small may '',
     .        ''severely compromise fidelity of surface geometry'',/)')
              end if
              read(iunit5,*)
           end if
         else
            if (nskip.ne.ngrid.and.nskip.ne.0) then
             write(iunit11,'(i8,i9,4f9.6,i9)') nskip,isktyp,beta1,
     .              alpha1,beta2,alpha2,nsprgit
             write(iunit11,'('' input error: nskip .ne. ngrid.'',
     .        '' May mean improper isktyp value.'',/)')
             call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
            end if
            write(iunit11,'(i8,i9,4f9.6,i9)') ngrid,isktyp,
     .              real(beta1),real(alpha1),real(beta2),real(alpha2),
     .              nsprgit
            if (real(beta1).lt.0.2) then
               write(iunit11,'(/,'' WARNING: beta this small may '',
     .      ''severely compromise fidelity of surface geometry'',/)')
            end if
            read(iunit5,*)
            if(nskip.eq.ngrid) then
             do ng = 1,nskip
               read(iunit5,*)
               read(iunit5,*) igr,nskpi1(ng),nskpj1(ng),nskpk1(ng)
               if (igr.ne.ng) then
                  write(iunit11,'(/,'' ERROR: control point input '',
     .            ''must include all blocks, and be entered in '',/,
     .            ''ascending block order'',/)')
                  call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
               end if
               if (nskpi1(ng).gt.500) then
                  write(iunit11,'(/,'' ERROR: control point input '',
     .            ''number of control points must be <= 500'',/)')
                  call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
               end if
               if (nskpj1(ng).gt.500) then
                  write(iunit11,'(/,'' ERROR: control point input '',
     .            ''number of control points must be <= 500'',/)')
                  call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
               end if
               if (nskpk1(ng).gt.500) then
                  write(iunit11,'(/,'' ERROR: control point input '',
     .            ''number of control points must be <= 500'',/)')
                  call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
               end if
               nbl = nblg(igr)
               ncg = ncgg(igr)
               read(iunit5,*)
               iis = -9
               iie =  0
               do jj = 1,50
                 iis = iis + 10
                 iie = iie + 10
                 if(iie.gt.nskpi1(ng)) iie = nskpi1(ng)
                 read(iunit5,*) (iskip(nbl,ii),ii=iis,iie)
                 if(iie.eq.nskpi1(ng)) goto 325
               enddo
325            continue
               read(iunit5,*)
               iis = -9
               iie =  0
               do jj = 1,50
                 iis = iis + 10
                 iie = iie + 10
                 if(iie.gt.nskpj1(ng)) iie = nskpj1(ng)
                 read(iunit5,*) (jskip(nbl,ii),ii=iis,iie)
                 if(iie.eq.nskpj1(ng)) goto 350
               enddo
350            continue
               read(iunit5,*)
               iis = -9
               iie =  0
               do jj = 1,50
                 iis = iis + 10
                 iie = iie + 10
                 if(iie.gt.nskpk1(ng)) iie = nskpk1(ng)
                 read(iunit5,*) (kskip(nbl,ii),ii=iis,iie)
                 if(iie.eq.nskpk1(ng)) goto 375
               enddo
375            continue
               jdim   = jdimg(nbl)
               kdim   = kdimg(nbl)
               idim   = idimg(nbl)
               if(iskip(nbl,nskpi1(ng)).ne.idim) then
               write(iunit11,'(/,'' STOPPING: ic(iskip) ne idim'',/)')
                  stop
               end if
               if(jskip(nbl,nskpj1(ng)).ne.jdim) then
               write(iunit11,'(/,'' STOPPING: ic(jskip) ne jdim'',/)')
                  stop
               end if
               if(kskip(nbl,nskpk1(ng)).ne.kdim) then
               write(iunit11,'(/,'' STOPPING: ic(kskip) ne kdim'',/)')
                  stop
               end if
               do jj = 1,nskpi1(ng)
                 iskipt(nbl,jj) = iskip(nbl,jj)
               enddo
               do jj = 1,nskpj1(ng)
                 jskipt(nbl,jj) = jskip(nbl,jj)
               enddo
               do jj = 1,nskpk1(ng)
                 kskipt(nbl,jj) = kskip(nbl,jj)
               enddo
             enddo
            end if
            if(nskip.eq.0) then 
             do ng = 1,ngrid    
               nbl          = nblg(ng)
               ncg          = ncgg(ng)
               jdim   = jdimg(nbl)
               kdim   = kdimg(nbl)
               idim   = idimg(nbl)
               nskpi1(ng)   = 2
               nskpj1(ng)   = 2
               nskpk1(ng)   = 2
               iskip(nbl,1) = 1
               iskip(nbl,2) = idim 
               jskip(nbl,1) = 1
               jskip(nbl,2) = jdim 
               kskip(nbl,1) = 1
               kskip(nbl,2) = kdim 
             enddo  
             do n = 1,abs(nbli)
               do ii = 1,2
                nbl    = nblg(nblk(ii,n))
                jdim   = jdimg(nbl)
                kdim   = kdimg(nbl)
                idim   = idimg(nbl)
                do jj = 1,2
                 if(isva(ii,jj,n).eq.1) then
                   if(limblk(ii,1,n).gt.1.and.limblk(ii,1,n).lt.idim)
     .                call skordr(nbl,nskpi1(nblk(ii,n)),
     .                            limblk(ii,1,n),iskip,maxbl) 
                   if(limblk(ii,4,n).gt.1.and.limblk(ii,4,n).lt.idim)
     .                call skordr(nbl,nskpi1(nblk(ii,n)),
     .                            limblk(ii,4,n),iskip,maxbl) 
                 else if(isva(ii,jj,n).eq.2) then
                   if(limblk(ii,2,n).gt.1.and.limblk(ii,2,n).lt.jdim)
     .                call skordr(nbl,nskpj1(nblk(ii,n)),
     .                            limblk(ii,2,n),jskip,maxbl) 
                   if(limblk(ii,5,n).gt.1.and.limblk(ii,5,n).lt.jdim)
     .                call skordr(nbl,nskpj1(nblk(ii,n)),
     .                              limblk(ii,5,n),jskip,maxbl) 
                 else if(isva(ii,jj,n).eq.3) then
                   if(limblk(ii,3,n).gt.1.and.limblk(ii,3,n).lt.kdim)
     .                call skordr(nbl,nskpk1(nblk(ii,n)),
     .                            limblk(ii,3,n),kskip,maxbl) 
                   if(limblk(ii,6,n).gt.1.and.limblk(ii,6,n).lt.kdim)
     .                call skordr(nbl,nskpk1(nblk(ii,n)),
     .                            limblk(ii,6,n),kskip,maxbl) 
                 end if
                enddo
               enddo
             enddo
             do n=1,ngrid
               nbl = nblg(n)
               jdim   = jdimg(nbl)
               kdim   = kdimg(nbl)
               idim   = idimg(nbl)
c
               do iseg = 1,nsegdfrm(nbl)
c
c icsi surface boundary
c
                 if(icsi(nbl,iseg).gt.1.and.icsi(nbl,iseg).lt.idim) 
     .              then 
                  call skordr(nbl,nskpi1(n),icsi(nbl,iseg)
     .                        ,iskip,maxbl) 
c
c          Check J0 boundary
c
                  if(jcsi(nbl,iseg).eq.jcsf(nbl,iseg).and.
     .               jcsi(nbl,iseg).eq.1) then
                    call bndchk(n,nbl,iseg,2,3,1,maxbl,maxsegdg,maxseg,
     .                          iskip,nskpi1,jbcinfo,nbcj0,icsi) 
                  end if
c
c          Check JDIM boundary
c
                  if(jcsi(nbl,iseg).eq.jcsf(nbl,iseg).and.
     .               jcsi(nbl,iseg).eq.jdim) then
                    call bndchk(n,nbl,iseg,2,3,2,maxbl,maxsegdg,maxseg,
     .                          iskip,nskpi1,jbcinfo,nbcjdim,icsi) 
                  end if
c
c          Check K0 boundary
c
                  if(kcsi(nbl,iseg).eq.kcsf(nbl,iseg).and.
     .               kcsi(nbl,iseg).eq.1) then
                    call bndchk(n,nbl,iseg,2,3,1,maxbl,maxsegdg,maxseg,
     .                          iskip,nskpi1,kbcinfo,nbck0,icsi) 
                  end if
c
c          Check KDIM boundary
c
                  if(kcsi(nbl,iseg).eq.kcsf(nbl,iseg).and.
     .               kcsi(nbl,iseg).eq.kdim) then
                    call bndchk(n,nbl,iseg,2,3,2,maxbl,maxsegdg,maxseg,
     .                          iskip,nskpi1,kbcinfo,nbckdim,icsi) 
                  end if
                 end if
c
c icsf surface boundary
c
                 if(icsf(nbl,iseg).gt.1.and.icsf(nbl,iseg).lt.idim)
     .              then 
                  call skordr(nbl,nskpi1(n),icsf(nbl,iseg)
     .                        ,iskip,maxbl) 
c
c          Check J0 boundary
c
                  if(jcsi(nbl,iseg).eq.jcsf(nbl,iseg).and.
     .               jcsi(nbl,iseg).eq.1) then
                    call bndchk(n,nbl,iseg,2,3,1,maxbl,maxsegdg,maxseg,
     .                          iskip,nskpi1,jbcinfo,nbcj0,icsf) 
                  end if
c
c          Check JDIM boundary
c
                  if(jcsi(nbl,iseg).eq.jcsf(nbl,iseg).and.
     .               jcsi(nbl,iseg).eq.jdim) then
                    call bndchk(n,nbl,iseg,2,3,2,maxbl,maxsegdg,maxseg,
     .                          iskip,nskpi1,jbcinfo,nbcjdim,icsf) 
                  end if
c
c          Check K0 boundary
c
                  if(kcsi(nbl,iseg).eq.kcsf(nbl,iseg).and.
     .               kcsi(nbl,iseg).eq.1) then
                    call bndchk(n,nbl,iseg,2,3,1,maxbl,maxsegdg,maxseg,
     .                          iskip,nskpi1,kbcinfo,nbck0,icsf) 
                  end if
c
c          Check KDIM boundary
c
                  if(kcsi(nbl,iseg).eq.kcsf(nbl,iseg).and.
     .               kcsi(nbl,iseg).eq.kdim) then
                    call bndchk(n,nbl,iseg,2,3,2,maxbl,maxsegdg,maxseg,
     .                          iskip,nskpi1,kbcinfo,nbckdim,icsf) 
                  end if
                 end if
c
c jcsi surface boundary
c
                 if(jcsi(nbl,iseg).gt.1.and.jcsi(nbl,iseg).lt.jdim)
     .              then 
                  call skordr(nbl,nskpj1(n),jcsi(nbl,iseg)
     .                        ,jskip,maxbl) 
c
c          Check K0 boundary
c
                  if(kcsi(nbl,iseg).eq.kcsf(nbl,iseg).and.
     .               kcsi(nbl,iseg).eq.1) then
                    call bndchk(n,nbl,iseg,4,5,1,maxbl,maxsegdg,maxseg,
     .                          jskip,nskpj1,kbcinfo,nbck0,jcsi) 
                  end if
c
c          Check KDIM boundary
c
                  if(kcsi(nbl,iseg).eq.kcsf(nbl,iseg).and.
     .               kcsi(nbl,iseg).eq.kdim) then
                    call bndchk(n,nbl,iseg,4,5,2,maxbl,maxsegdg,maxseg,
     .                          jskip,nskpj1,kbcinfo,nbckdim,jcsi) 
                  end if
c
c          Check I0 boundary
c
                  if(icsi(nbl,iseg).eq.icsf(nbl,iseg).and.
     .               icsi(nbl,iseg).eq.1) then
                    call bndchk(n,nbl,iseg,2,3,1,maxbl,maxsegdg,maxseg,
     .                          jskip,nskpj1,ibcinfo,nbci0,jcsi) 
                  end if
c
c          Check IDIM boundary
c
                  if(icsi(nbl,iseg).eq.icsf(nbl,iseg).and.
     .               icsi(nbl,iseg).eq.idim) then
                    call bndchk(n,nbl,iseg,2,3,2,maxbl,maxsegdg,maxseg,
     .                          jskip,nskpj1,ibcinfo,nbcidim,jcsi) 
                  end if
                 end if
c
c jcsf surface boundary
c
                 if(jcsf(nbl,iseg).gt.1.and.jcsf(nbl,iseg).lt.jdim)
     .              then 
                  call skordr(nbl,nskpj1(n),jcsf(nbl,iseg)
     .                        ,jskip,maxbl) 
c
c          Check K0 boundary
c
                  if(kcsi(nbl,iseg).eq.kcsf(nbl,iseg).and.
     .               kcsi(nbl,iseg).eq.1) then
                    call bndchk(n,nbl,iseg,4,5,1,maxbl,maxsegdg,maxseg,
     .                          jskip,nskpj1,kbcinfo,nbck0,jcsf) 
                  end if
c
c          Check KDIM boundary
c
                  if(kcsi(nbl,iseg).eq.kcsf(nbl,iseg).and.
     .               kcsi(nbl,iseg).eq.kdim) then
                    call bndchk(n,nbl,iseg,4,5,2,maxbl,maxsegdg,maxseg,
     .                          jskip,nskpj1,kbcinfo,nbckdim,jcsf) 
                  end if
c
c          Check I0 boundary
c
                  if(icsi(nbl,iseg).eq.icsf(nbl,iseg).and.
     .               icsi(nbl,iseg).eq.1) then
                    call bndchk(n,nbl,iseg,2,3,1,maxbl,maxsegdg,maxseg,
     .                          jskip,nskpj1,ibcinfo,nbci0,jcsf) 
                  end if
c
c          Check IDIM boundary
c
                  if(icsi(nbl,iseg).eq.icsf(nbl,iseg).and.
     .               icsi(nbl,iseg).eq.idim) then
                    call bndchk(n,nbl,iseg,2,3,2,maxbl,maxsegdg,maxseg,
     .                          jskip,nskpj1,ibcinfo,nbcidim,jcsf) 
                  end if
                 end if
c
c kcsi surface boundary
c
                 if(kcsi(nbl,iseg).gt.1.and.kcsi(nbl,iseg).lt.kdim)
     .              then 
                  call skordr(nbl,nskpk1(n),kcsi(nbl,iseg)
     .                        ,kskip,maxbl) 
c
c          Check I0 boundary
c
                  if(icsi(nbl,iseg).eq.icsf(nbl,iseg).and.
     .               icsi(nbl,iseg).eq.1) then
                    call bndchk(n,nbl,iseg,4,5,1,maxbl,maxsegdg,maxseg,
     .                          kskip,nskpk1,ibcinfo,nbci0,kcsi) 
                  end if
c
c          Check IDIM boundary
c
                  if(icsi(nbl,iseg).eq.icsf(nbl,iseg).and.
     .               icsi(nbl,iseg).eq.idim) then
                    call bndchk(n,nbl,iseg,4,5,2,maxbl,maxsegdg,maxseg,
     .                          kskip,nskpk1,ibcinfo,nbcidim,kcsi) 
                  end if
c
c          Check J0 boundary
c
                  if(jcsi(nbl,iseg).eq.jcsf(nbl,iseg).and.
     .               jcsi(nbl,iseg).eq.1) then
                    call bndchk(n,nbl,iseg,4,5,1,maxbl,maxsegdg,maxseg,
     .                          kskip,nskpk1,jbcinfo,nbcj0,kcsi) 
                  end if
c
c          Check JDIM boundary
c
                  if(jcsi(nbl,iseg).eq.jcsf(nbl,iseg).and.
     .               jcsi(nbl,iseg).eq.jdim) then
                    call bndchk(n,nbl,iseg,4,5,2,maxbl,maxsegdg,maxseg,
     .                          kskip,nskpk1,jbcinfo,nbcjdim,kcsi) 
                  end if
                 end if
c
c kcsf surface boundary
c
                 if(kcsf(nbl,iseg).gt.1.and.kcsf(nbl,iseg).lt.kdim)
     .              then 
                  call skordr(nbl,nskpk1(n),kcsf(nbl,iseg)
     .                        ,kskip,maxbl) 
c
c          Check I0 boundary
c
                  if(icsi(nbl,iseg).eq.icsf(nbl,iseg).and.
     .               icsi(nbl,iseg).eq.1) then
                    call bndchk(n,nbl,iseg,4,5,1,maxbl,maxsegdg,maxseg,
     .                          kskip,nskpk1,ibcinfo,nbci0,kcsf) 
                  end if
c
c          Check IDIM boundary
c
                  if(icsi(nbl,iseg).eq.icsf(nbl,iseg).and.
     .               icsi(nbl,iseg).eq.idim) then
                    call bndchk(n,nbl,iseg,4,5,2,maxbl,maxsegdg,maxseg,
     .                          kskip,nskpk1,ibcinfo,nbcidim,kcsf) 
                  end if
c
c          Check J0 boundary
c
                  if(jcsi(nbl,iseg).eq.jcsf(nbl,iseg).and.
     .               jcsi(nbl,iseg).eq.1) then
                    call bndchk(n,nbl,iseg,4,5,1,maxbl,maxsegdg,maxseg,
     .                          kskip,nskpk1,jbcinfo,nbcj0,kcsf) 
                  end if
c
c          Check JDIM boundary
c
                  if(jcsi(nbl,iseg).eq.jcsf(nbl,iseg).and.
     .               jcsi(nbl,iseg).eq.jdim) then
                    call bndchk(n,nbl,iseg,4,5,2,maxbl,maxsegdg,maxseg,
     .                          kskip,nskpk1,jbcinfo,nbcjdim,kcsf) 
                  end if
                 end if
               enddo
             enddo
             do n = 1,ngrid
               nbl = nblg(n)
               jdim   = jdimg(nbl)
               kdim   = kdimg(nbl)
               idim   = idimg(nbl)
               do iseg=1,nbci0(nbl)
                 if   (ibcinfo(nbl,iseg,2,1).gt.1.
     .             and.ibcinfo(nbl,iseg,2,1).lt.jdim) 
     .             call skordr(nbl,nskpj1(n),ibcinfo(nbl,iseg,2,1)
     .                        ,jskip,maxbl) 
                 if   (ibcinfo(nbl,iseg,3,1).gt.1.
     .             and.ibcinfo(nbl,iseg,3,1).lt.jdim) 
     .             call skordr(nbl,nskpj1(n),ibcinfo(nbl,iseg,3,1)
     .                        ,jskip,maxbl) 
                 if   (ibcinfo(nbl,iseg,4,1).gt.1.
     .             and.ibcinfo(nbl,iseg,4,1).lt.kdim) 
     .             call skordr(nbl,nskpk1(n),ibcinfo(nbl,iseg,4,1)
     .                        ,kskip,maxbl) 
                 if   (ibcinfo(nbl,iseg,5,1).gt.1.
     .             and.ibcinfo(nbl,iseg,5,1).lt.kdim) 
     .             call skordr(nbl,nskpk1(n),ibcinfo(nbl,iseg,5,1)
     .                        ,kskip,maxbl) 
               enddo
               do iseg=1,nbcidim(nbl)
                 if   (ibcinfo(nbl,iseg,2,2).gt.1.
     .             and.ibcinfo(nbl,iseg,2,2).lt.jdim) 
     .             call skordr(nbl,nskpj1(n),ibcinfo(nbl,iseg,2,2)
     .                        ,jskip,maxbl) 
                 if   (ibcinfo(nbl,iseg,3,2).gt.1.
     .             and.ibcinfo(nbl,iseg,3,2).lt.jdim) 
     .             call skordr(nbl,nskpj1(n),ibcinfo(nbl,iseg,3,2)
     .                        ,jskip,maxbl) 
                 if   (ibcinfo(nbl,iseg,4,2).gt.1.
     .             and.ibcinfo(nbl,iseg,4,2).lt.kdim) 
     .             call skordr(nbl,nskpk1(n),ibcinfo(nbl,iseg,4,2)
     .                        ,kskip,maxbl) 
                 if   (ibcinfo(nbl,iseg,5,2).gt.1.
     .             and.ibcinfo(nbl,iseg,5,2).lt.kdim) 
     .             call skordr(nbl,nskpk1(n),ibcinfo(nbl,iseg,5,2)
     .                        ,kskip,maxbl) 
               enddo
               do iseg=1,nbcj0(nbl)
                 if   (jbcinfo(nbl,iseg,2,1).gt.1.
     .             and.jbcinfo(nbl,iseg,2,1).lt.idim) 
     .             call skordr(nbl,nskpi1(n),jbcinfo(nbl,iseg,2,1)
     .                        ,iskip,maxbl) 
                 if   (jbcinfo(nbl,iseg,3,1).gt.1.
     .             and.jbcinfo(nbl,iseg,3,1).lt.idim) 
     .             call skordr(nbl,nskpi1(n),jbcinfo(nbl,iseg,3,1)
     .                        ,iskip,maxbl) 
                 if   (jbcinfo(nbl,iseg,4,1).gt.1.
     .             and.jbcinfo(nbl,iseg,4,1).lt.kdim) 
     .             call skordr(nbl,nskpk1(n),jbcinfo(nbl,iseg,4,1)
     .                        ,kskip,maxbl) 
                 if   (jbcinfo(nbl,iseg,5,1).gt.1.
     .             and.jbcinfo(nbl,iseg,5,1).lt.kdim) 
     .             call skordr(nbl,nskpk1(n),jbcinfo(nbl,iseg,5,1)
     .                        ,kskip,maxbl) 
               enddo
               do iseg=1,nbcjdim(nbl)
                 if   (jbcinfo(nbl,iseg,2,2).gt.1.
     .             and.jbcinfo(nbl,iseg,2,2).lt.idim) 
     .             call skordr(nbl,nskpi1(n),jbcinfo(nbl,iseg,2,2)
     .                        ,iskip,maxbl) 
                 if   (jbcinfo(nbl,iseg,3,2).gt.1.
     .             and.jbcinfo(nbl,iseg,3,2).lt.idim) 
     .             call skordr(nbl,nskpi1(n),jbcinfo(nbl,iseg,3,2)
     .                        ,iskip,maxbl) 
                 if   (jbcinfo(nbl,iseg,4,2).gt.1.
     .             and.jbcinfo(nbl,iseg,4,2).lt.kdim) 
     .             call skordr(nbl,nskpk1(n),jbcinfo(nbl,iseg,4,2)
     .                        ,kskip,maxbl) 
                 if   (jbcinfo(nbl,iseg,5,2).gt.1.
     .             and.jbcinfo(nbl,iseg,5,2).lt.kdim) 
     .             call skordr(nbl,nskpk1(n),jbcinfo(nbl,iseg,5,2)
     .                        ,kskip,maxbl) 
               enddo
               do iseg=1,nbck0(nbl)
                 if   (kbcinfo(nbl,iseg,2,1).gt.1.
     .             and.kbcinfo(nbl,iseg,2,1).lt.idim) 
     .             call skordr(nbl,nskpi1(n),kbcinfo(nbl,iseg,2,1)
     .                        ,iskip,maxbl) 
                 if   (kbcinfo(nbl,iseg,3,1).gt.1.
     .             and.kbcinfo(nbl,iseg,3,1).lt.idim) 
     .             call skordr(nbl,nskpi1(n),kbcinfo(nbl,iseg,3,1)
     .                        ,iskip,maxbl) 
                 if   (kbcinfo(nbl,iseg,4,1).gt.1.
     .             and.kbcinfo(nbl,iseg,4,1).lt.jdim) 
     .             call skordr(nbl,nskpj1(n),kbcinfo(nbl,iseg,4,1)
     .                        ,jskip,maxbl) 
                 if   (kbcinfo(nbl,iseg,5,1).gt.1.
     .             and.kbcinfo(nbl,iseg,5,1).lt.jdim) 
     .             call skordr(nbl,nskpj1(n),kbcinfo(nbl,iseg,5,1)
     .                        ,jskip,maxbl) 
               enddo
               do iseg=1,nbckdim(nbl)
                 if   (kbcinfo(nbl,iseg,2,2).gt.1.
     .             and.kbcinfo(nbl,iseg,2,2).lt.idim) 
     .             call skordr(nbl,nskpi1(n),kbcinfo(nbl,iseg,2,2)
     .                        ,iskip,maxbl) 
                 if   (kbcinfo(nbl,iseg,3,2).gt.1.
     .             and.kbcinfo(nbl,iseg,3,2).lt.idim) 
     .             call skordr(nbl,nskpi1(n),kbcinfo(nbl,iseg,3,2)
     .                        ,iskip,maxbl) 
                 if   (kbcinfo(nbl,iseg,4,2).gt.1.
     .             and.kbcinfo(nbl,iseg,4,2).lt.jdim) 
     .             call skordr(nbl,nskpj1(n),kbcinfo(nbl,iseg,4,2)
     .                        ,jskip,maxbl) 
                 if   (kbcinfo(nbl,iseg,5,2).gt.1.
     .             and.kbcinfo(nbl,iseg,5,2).lt.jdim) 
     .             call skordr(nbl,nskpj1(n),kbcinfo(nbl,iseg,5,2)
     .                        ,jskip,maxbl) 
               enddo
             enddo
            end if
            do iter = 1,ngrid
              nskpsmo = 0
              nskpsmn = 0
              itcpadd = itcpadd + 1
              do n = 1,ngrid
                nskpsmo = nskpsmo+nskpi1(n)+nskpj1(n)+nskpk1(n) 
              enddo
              do n = 1,abs(nbli)
               nbl1   = nblg(nblk(1,n))
               jdim1  = jdimg(nbl1)
               kdim1  = kdimg(nbl1)
               idim1  = idimg(nbl1)
               nbl2   = nblg(nblk(2,n))
               jdim2  = jdimg(nbl2)
               kdim2  = kdimg(nbl2)
               idim2  = idimg(nbl2) 
               do jj = 1,2    
                if(isva(1,jj,n).eq.1.and.isva(2,jj,n).eq.1) then 
                  call cpadd(nbl1,nbl2,nblk(1,n),nblk(2,n),
     .                       limblk(1,1,n),limblk(1,4,n),limblk(2,1,n),
     .                       limblk(2,4,n),nskpi1(nblk(1,n)),
     .                       iskip(nbl1,1:500),nskpi1(nblk(2,n)), 
     .                       iskip(nbl2,1:500),maxbl,idim1,idim2,1,1)
                else          
     .          if(isva(1,jj,n).eq.1.and.isva(2,jj,n).eq.2) then 
                  call cpadd(nbl1,nbl2,nblk(1,n),nblk(2,n),
     .                       limblk(1,1,n),limblk(1,4,n),limblk(2,2,n),
     .                       limblk(2,5,n),nskpi1(nblk(1,n)),
     .                       iskip(nbl1,1:500),nskpj1(nblk(2,n)), 
     .                       jskip(nbl2,1:500),maxbl,idim1,jdim2,1,2)
                else          
     .          if(isva(1,jj,n).eq.2.and.isva(2,jj,n).eq.1) then 
                  call cpadd(nbl1,nbl2,nblk(1,n),nblk(2,n),
     .                       limblk(1,2,n),limblk(1,5,n),limblk(2,1,n),
     .                       limblk(2,4,n),nskpj1(nblk(1,n)),
     .                       jskip(nbl1,1:500),nskpi1(nblk(2,n)),
     .                       iskip(nbl2,1:500),maxbl,jdim1,idim2,2,1)
                else
     .          if(isva(1,jj,n).eq.2.and.isva(2,jj,n).eq.2) then
                  call cpadd(nbl1,nbl2,nblk(1,n),nblk(2,n),
     .                       limblk(1,2,n),limblk(1,5,n),limblk(2,2,n),
     .                       limblk(2,5,n),nskpj1(nblk(1,n)),
     .                       jskip(nbl1,1:500),nskpj1(nblk(2,n)),
     .                       jskip(nbl2,1:500),maxbl,jdim1,jdim2,2,2)
                else
     .          if(isva(1,jj,n).eq.3.and.isva(2,jj,n).eq.2) then
                  call cpadd(nbl1,nbl2,nblk(1,n),nblk(2,n),
     .                       limblk(1,3,n),limblk(1,6,n),limblk(2,2,n),
     .                       limblk(2,5,n),nskpk1(nblk(1,n)),
     .                       kskip(nbl1,1:500),nskpj1(nblk(2,n)),
     .                       jskip(nbl2,1:500),maxbl,kdim1,jdim2,3,2)
                else
     .          if(isva(1,jj,n).eq.2.and.isva(2,jj,n).eq.3) then
                  call cpadd(nbl1,nbl2,nblk(1,n),nblk(2,n),
     .                       limblk(1,2,n),limblk(1,5,n),limblk(2,3,n),
     .                       limblk(2,6,n),nskpj1(nblk(1,n)),
     .                       jskip(nbl1,1:500),nskpk1(nblk(2,n)),
     .                       kskip(nbl2,1:500),maxbl,jdim1,kdim2,2,3)
                else
     .          if(isva(1,jj,n).eq.3.and.isva(2,jj,n).eq.3) then
                  call cpadd(nbl1,nbl2,nblk(1,n),nblk(2,n),
     .                       limblk(1,3,n),limblk(1,6,n),limblk(2,3,n),
     .                       limblk(2,6,n),nskpk1(nblk(1,n)),
     .                       kskip(nbl1,1:500),nskpk1(nblk(2,n)),
     .                       kskip(nbl2,1:500),maxbl,kdim1,kdim2,3,3)
                else
     .          if(isva(1,jj,n).eq.3.and.isva(2,jj,n).eq.1) then
                  call cpadd(nbl1,nbl2,nblk(1,n),nblk(2,n),
     .                       limblk(1,3,n),limblk(1,6,n),limblk(2,1,n),
     .                       limblk(2,4,n),nskpk1(nblk(1,n)),
     .                       kskip(nbl1,1:500),nskpi1(nblk(2,n)),
     .                       iskip(nbl2,1:500),maxbl,kdim1,idim2,3,1)
                else
     .          if(isva(1,jj,n).eq.1.and.isva(2,jj,n).eq.3) then
                  call cpadd(nbl1,nbl2,nblk(1,n),nblk(2,n),
     .                       limblk(1,1,n),limblk(1,4,n),limblk(2,3,n),
     .                       limblk(2,6,n),nskpi1(nblk(1,n)),
     .                       iskip(nbl1,1:500),nskpk1(nblk(2,n)),
     .                       kskip(nbl2,1:500),maxbl,idim1,kdim2,1,3)
                end if
               enddo
              enddo
              do n = 1,ngrid
                nskpsmn = nskpsmn+nskpi1(n)+nskpj1(n)+nskpk1(n) 
              enddo
              if(nskpsmn.eq.nskpsmo) goto 378
            enddo
378         continue 
            do n = 1,ngrid    
              nbl   = nblg(n)
              do jj = 1,nskpi1(n)
                iskipt(nbl,jj) = iskip(nbl,jj)
              enddo
              do jj = 1,nskpj1(n)
                jskipt(nbl,jj) = jskip(nbl,jj)
              enddo
              do jj = 1,nskpk1(n)
                kskipt(nbl,jj) = kskip(nbl,jj)
              enddo
            enddo
         end if
c
         if(abs(isktyp).eq.1) then
           write(iunit11,'(''    grid    iskip    jskip    kskip'')')
           do ng=1,ngrid
             nbl = nblg(ng)
             ncg = ncgg(ng)
             write(iunit11,'(i8,6i9)') ng,iskip(nbl,1),jskip(nbl,1),
     .                                  kskip(nbl,1)
             if (ncg.gt.0) then
                do n=1,ncg
                   nbl = nbl + 1
                   if(idim.gt.2) then
                     iskip(nbl,1) = iskip(nbl-1,1)/2
                   else
                     iskip(nbl,1) = iskip(nbl-1,1)
                   end if
                   jskip(nbl,1) = jskip(nbl-1,1)/2
                   kskip(nbl,1) = kskip(nbl-1,1)/2
                end do
                nbl = nblg(ng)
             end if
             ncgt = ncg
             if(idim.eq.2) then
               nskpi1(ng) = 2
             else
               nskpi1(ng) = (idimg(nbl)-1)/iskip(nbl,1) + 1
             end if
             nskpj1(ng) = (jdimg(nbl)-1)/jskip(nbl,1) + 1
             nskpk1(ng) = (kdimg(nbl)-1)/kskip(nbl,1) + 1
             iskipt(nbl,1) = 1
             jskipt(nbl,1) = 1
             kskipt(nbl,1) = 1
             if (nskpi1(ng).gt.500) then
                write(iunit11,'(/,'' ERROR: nskpi1 > 500 '',
     .          '' Too many i-control points (iskip too small)'',/)')
                call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
             end if
             if (nskpj1(ng).gt.500) then
                write(iunit11,'(/,'' ERROR: nskpj1 > 500 '',
     .          '' Too many j-control points (jskip too small)'',/)')
                call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
             end if
             if (nskpk1(ng).gt.500) then
                write(iunit11,'(/,'' ERROR: nskpk1 > 500 '',
     .          '' Too many k-control points (kskip too small)'',/)')
                call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
             end if
             do jj = 2,nskpi1(ng)
               iskipt(nbl,jj) = iskipt(nbl,jj-1)+iskip(nbl,1)
             enddo
             do jj = 2,nskpj1(ng)
               jskipt(nbl,jj) = jskipt(nbl,jj-1)+jskip(nbl,1)
             enddo
             do jj = 2,nskpk1(ng)
               kskipt(nbl,jj) = kskipt(nbl,jj-1)+kskip(nbl,1)
             enddo
           enddo
         else
           do ng=1,ngrid
              nbl = nblg(ng)
              ncg = ncgg(ng)
              if (ncg.gt.0) then
                nbl = nblg(ng)
                do ii = 1,nskpi1(ng)
                  iskipt(nbl,ii)= iskip(nbl,ii)
                enddo
                if(idim.gt.2) then
                  do n=1,ncg
                    nbl = nbl + 1
                    do ii = 1,nskpi1(ng)
                      iskip(nbl,ii) = (iskip(nbl-1,ii)-1)/2 + 1
                      iskipt(nbl,ii)= iskip(nbl,ii)
                      riskp(nbl,ii) = .5*real(iskip(nbl-1,ii)-1) + 1.0  
                    end do
                    do ii = 2,nskpi1(ng)-1  
                      if(iskip(nbl,ii).ge.iskip(nbl,ii+1).and.
     .                   real(iskip(nbl,ii)).ne.riskp(nbl,ii)) 
     .                    iskip(nbl,ii) = iskip(nbl,ii)-1
                      if(iskip(nbl,ii).le.iskip(nbl,ii-1).and.
     .                   real(iskip(nbl,ii)).ne.riskp(nbl,ii))  
     .                    iskip(nbl,ii) = iskip(nbl,ii)+1
                    enddo
                    do ii = 1,nskpi1(ng)
                      iskipt(nbl,ii)= iskip(nbl,ii)
                    enddo
                  enddo
                else
                  do n=1,ncg
                    nbl = nbl + 1
                    do ii = 1,nskpi1(ng)
                      iskip(nbl,ii) = iskip(nbl-1,ii)
                    end do
                    do ii = 1,nskpi1(ng)
                      iskipt(nbl,ii)= iskip(nbl,ii)
                    enddo
                  enddo
                end if 
                nbl = nblg(ng)
                do ii = 1,nskpj1(ng)
                  jskipt(nbl,ii)= jskip(nbl,ii)
                enddo
                do n=1,ncg
                  nbl = nbl + 1
                  do ii = 1,nskpj1(ng)
                    jskip(nbl,ii) = (jskip(nbl-1,ii)-1)/2 + 1
                    rjskp(nbl,ii) = .5*real(jskip(nbl-1,ii)-1) + 1.0
                  end do
                  do ii = 2,nskpj1(ng)-1
                    if(jskip(nbl,ii).ge.jskip(nbl,ii+1).and.
     .                 real(jskip(nbl,ii)).ne.rjskp(nbl,ii))
     .                  jskip(nbl,ii) = jskip(nbl,ii)-1
                    if(jskip(nbl,ii).le.jskip(nbl,ii-1).and.
     .                 real(jskip(nbl,ii)).ne.rjskp(nbl,ii))
     .                  jskip(nbl,ii) = jskip(nbl,ii)+1
                  enddo
                  do ii = 1,nskpj1(ng)
                    jskipt(nbl,ii)= jskip(nbl,ii)
                  enddo
                enddo
                nbl = nblg(ng)
                do ii = 1,nskpk1(ng)
                  kskipt(nbl,ii)= kskip(nbl,ii)
                enddo
                do n=1,ncg
                  nbl = nbl + 1
                  do ii = 1,nskpk1(ng)
                    kskip(nbl,ii) = (kskip(nbl-1,ii)-1)/2 + 1
                    rkskp(nbl,ii) = .5*real(kskip(nbl-1,ii)-1) + 1.0
                  end do
                  do ii = 2,nskpk1(ng)-1
                    if(kskip(nbl,ii).ge.kskip(nbl,ii+1).and.
     .                 real(kskip(nbl,ii)).ne.rkskp(nbl,ii))
     .                  kskip(nbl,ii) = kskip(nbl,ii)-1
                    if(kskip(nbl,ii).le.kskip(nbl,ii-1).and.
     .                 real(kskip(nbl,ii)).ne.rkskp(nbl,ii))
     .                  kskip(nbl,ii) = kskip(nbl,ii)+1
                  enddo
                  do ii = 1,nskpk1(ng)
                    kskipt(nbl,ii)= kskip(nbl,ii)
                  enddo
                enddo
              end if
           enddo
         end if
         if(abs(isktyp).eq.1.or.(abs(isktyp).eq.2.and.nskip.eq.0)
     .     .or.(itcpadd.gt.1)) then
          if(meshdef.eq.1) then
           open(196,file='meshdef.inp',status='unknown',form=
     .             'formatted')
           do ng=1,ngrid
             nbl = nblg(ng)
             ncg = ncgg(ng)
             write(196,19019)
19019        format('   GRID   NIND   NJND   NKND ')
             write(196,19020) ng,nskpi1(ng),nskpj1(ng),nskpk1(ng)
19020        format(10i7)
             write(196,19018)
19018        format('************************** I NODE INDICES ******'
     .             ,'**********************')
             iis = -9
             iie =  0
             do jj = 1,50
               iis = iis + 10
               iie = iie + 10
               if(iie.gt.nskpi1(ng)) iie = nskpi1(ng)
               write(196,19020) (iskipt(nbl,ii),ii=iis,iie)
               if(iie.eq.nskpi1(ng)) goto 525
             enddo
525          continue
             write(196,19017)
19017        format('************************** J NODE INDICES ******'
     .             ,'**********************')
             iis = -9
             iie =  0
             do jj = 1,50
               iis = iis + 10
               iie = iie + 10
               if(iie.gt.nskpj1(ng)) iie = nskpj1(ng)
               write(196,19020) (jskipt(nbl,ii),ii=iis,iie)
               if(iie.eq.nskpj1(ng)) goto 550
             enddo
550          continue
             write(196,19016)
19016        format('************************** K NODE INDICES ******'
     .             ,'**********************')
             iis = -9
             iie =  0
             do jj = 1,50
               iis = iis + 10
               iie = iie + 10
               if(iie.gt.nskpk1(ng)) iie = nskpk1(ng)
               write(196,19020) (kskipt(nbl,ii),ii=iis,iie)
               if(iie.eq.nskpk1(ng)) goto 575
             enddo
575          continue
           enddo
           close(196)
          end if
         end if
         if(abs(isktyp).eq.2) then
          do ng=1,ngrid
            nbl = nblg(ng)
            ncg = ncgg(ng)
            ncgt = ncg
            write(iunit11,
     .      '(''      ng     nipt     njpt     nkpt  '')')
            write(iunit11,'(i8,4i9)') ng,nskpi1(ng),nskpj1(ng)
     .                          ,nskpk1(ng)
            nbl = nblg(ng)
            write(iunit11,
     .    '(''    control point i-indices for grid levels  '',5i4)')
     .            (nbl+i-1,i=1,ncgt+1)
            do i = 1,nskpi1(ng)
              write(iunit11,'(i8,5i9)') (iskipt(nbl+ii-1,i),ii=1,ncgt+1)
            enddo
            nbl = nblg(ng)
            write(iunit11,
     .    '(''    control point j-indices for grid levels  '',5i4)')
     .            (nbl+i-1,i=1,ncgt+1)
            do i = 1,nskpj1(ng)
              write(iunit11,'(i8,5i9)') (jskipt(nbl+ii-1,i),ii=1,ncgt+1)
            enddo
            nbl = nblg(ng)
            write(iunit11,
     .    '(''    control point k-indices for grid levels  '',5i4)')
     .            (nbl+i-1,i=1,ncgt+1)
            do i = 1,nskpk1(ng)
              write(iunit11,'(i8,5i9)') (kskipt(nbl+ii-1,i),ii=1,ncgt+1)
            enddo
          enddo
         end if
c
c        for all blocks that have not yet been identified as 
c        deforming, set idefrm flag to 999...indicating the mesh 
c        deforms, but only through interaction with other blocks 
c        containing moving solid surfaces. also set index ranges
c        to include entire block
c
         do nbl=1,nblock
            if (idefrm(nbl).eq.0) then
               idefrm(nbl)  = 999
               nsegdfrm(nbl) = 1
               do iseg=1,nsegdfrm(nbl)
                  icsi(nbl,iseg) = 1
                  icsf(nbl,iseg) = idimg(nbl)
                  jcsi(nbl,iseg) = 1
                  jcsf(nbl,iseg) = jdimg(nbl)
                  kcsi(nbl,iseg) = 1
                  kcsf(nbl,iseg) = kdimg(nbl)
               end do
            end if
         end do
c
c        multi-mode grid movement data - used to couple rigid and
c        deforming grid motions, or deforming translation plus 
c        deforming rotation.
c
         read(iunit5,*)
         read(iunit5,*)
         read(iunit5,*) ncoupl
         write(iunit11,'('' moving grid data - multi-motion '',
     .                   ''coupling'')')
         write(iunit11,'(''  ncoupl'')')
         if (ncoupl.gt.0) then
            write(iunit11,'(i8)') ncoupl
            if (ncoupl.gt.ngrid) then
               write(iunit11,'('' input error: ncoupl > ngrid'')')
               call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
            end if
            read(iunit5,*)
            write(iunit11,'(''  slave   master   xorig   yorig'',
     .      ''   zorig'')')
            do nn=1,ncoupl
               read(iunit5,*) igslav,igmast,xorg,yorg,zorg
               write(iunit11,'(i8,i9,3f9.4)') igslav,igmast,xorg,
     .         yorg,zorg
               nblslav = nblg(igslav)
               if (igmast.gt.0) then
                  nblmast = nblg(igmast)
                  xorig(nblmast) = xorg
                  yorig(nblmast) = yorg
                  zorig(nblmast) = zorg
                  xorig0(nblmast) = xorg
                  yorig0(nblmast) = yorg
                  zorig0(nblmast) = zorg
               else
                  nblmast = 0
               end if
               ncg = ncgg(igslav)
               do iseg=1,nsegdfrm(nblslav)
                  icouple(nblslav,iseg) = nblmast
               end do
               if (ncg.gt.0) then
                   do n=1,ncg
                      nblslav = nblslav + 1
                      if (nblmast.gt.0) then
                         nblmast = nblmast + 1 
                         xorig(nblmast) = xorig(nblmast-1)
                         yorig(nblmast) = yorig(nblmast-1)
                         zorig(nblmast) = zorig(nblmast-1)
                         xorig0(nblmast) = xorig0(nblmast-1)
                         yorig0(nblmast) = yorig0(nblmast-1)
                         zorig0(nblmast) = zorig0(nblmast-1)
                      end if
                      do iseg=1,nsegdfrm(nblslav)
                         icouple(nblslav,iseg) = nblmast
                      end do
                   end do
               end if
            end do
         else if (ncoupl.lt.0) then
            read(iunit5,*) igslav,igmast,xorg,yorg,zorg
            ncoupl = ngrid
            write(iunit11,'(i8)') ncoupl
            read(iunit5,*)
            write(iunit11,'(''  slave   master   xorig   yorig'',
     .      ''   zorig'')')
            do nn=1,ncoupl
               igslav = nn
               write(iunit11,'(i8,i9,3f9.4)') igslav,igmast,xorg,
     .         yorg,zorg
               nblslav = nblg(igslav)
               if (igmast.gt.0) then
                  nblmast = nblg(igmast)
                  xorig(nblmast) = xorg
                  yorig(nblmast) = yorg
                  zorig(nblmast) = zorg
                  xorig0(nblmast) = xorg
                  yorig0(nblmast) = yorg
                  zorig0(nblmast) = zorg
               else
                  nblmast = 0
               end if
               ncg = ncgg(igslav)
               do iseg=1,nsegdfrm(nblslav)
                  icouple(nblslav,iseg) = nblmast
               end do
               if (ncg.gt.0) then
                   do n=1,ncg
                      nblslav = nblslav + 1
                      if (nblmast.gt.0) then
                         nblmast = nblmast + 1 
                         xorig(nblmast) = xorig(nblmast-1)
                         yorig(nblmast) = yorig(nblmast-1)
                         zorig(nblmast) = zorig(nblmast-1)
                         xorig0(nblmast) = xorig0(nblmast-1)
                         yorig0(nblmast) = yorig0(nblmast-1)
                         zorig0(nblmast) = zorig0(nblmast-1)
                      end if
                      do iseg=1,nsegdfrm(nblslav)
                         icouple(nblslav,iseg) = nblmast
                      end do
                   end do
               end if
            end do
         else
            write(iunit11,'(i8)') ncoupl
            read(iunit5,*)
            write(iunit11,'(''  slave   master   xorig   yorig'',
     .      ''   zorig'')')
         end if
c
         deallocate(riskp,rjskp,rkskp)
         deallocate(iskipt,jskipt,kskipt)
      end if
c
      do nbl=1,nblock
         mxdefseg = max(mxdefseg,nsegdfrm(nbl))
      end do
c
c     check that the deforming mesh segments are on the
c     appropriate block face (min or max of the direction
c     for which one of the indicies is constant
c
      do nbl=1,nblock
         jdim = jdimg(nbl)
         kdim = kdimg(nbl)
         idim = idimg(nbl)
         if (idefrm(nbl).lt.99) then
            do iseg=1,nsegdfrm(nbl)
               ist = icsi(nbl,iseg)
               ifn = icsf(nbl,iseg)
               jst = jcsi(nbl,iseg)
               jfn = jcsf(nbl,iseg)
               kst = kcsi(nbl,iseg)
               kfn = kcsf(nbl,iseg)
               iok = 1
               jok = 1
               kok = 1
               lconst = 1
               if (ist.eq.ifn) then
                  if (ist.eq.1 .or. ist.eq.idim) then
                     iok = 1
                  else
                     iok = 0
                  end if
               else if (jst.eq.jfn) then
                  if (jst.eq.1 .or. jst.eq.jdim) then
                     jok = 1
                  else
                     jok = 0
                  end if
               else if (kst.eq.kfn) then
                  if (kst.eq.1 .or. kst.eq.kdim) then
                     kok = 1
                  else
                     kok = 0
                  end if
               else
                  lconst = 0
               end if
               if (iok.eq.0 .or. jok.eq.0 .or. kok.eq.0) then
                  write(iunit11,'('' stopping...deforming surface'',
     .                         '' does not lie on a block face'')')
               end if
               if (lconst.eq.0) then
                  write(iunit11,'('' stopping...at least one index'',
     .                           '' must be constant on a deforming'',
     .                           '' surface'')')
               end if 
               isum = iok + jok + kok + lconst
               if (isum.ne.4) then
                  write(iunit11,'(''   segment'',i3,'' block '',
     .                            i4,'':'')') iseg,nbl
                  write(iunit11,'(''   icsi,icsf,jcsi,jcsf,kcsi,kcsf '',
     .                            6i4)') ist,ifn,jst,jfn,kst,kfn
                  call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
               end if
            end do
         end if
      end do
c
c     Now swap y and z grid motion parameters as needed for ialph > 0
c
      if (ialph > 0) then
c
c       rigid mesh movement data
c
        temp     = vtransmc
        vtransmc = -wtransmc
        wtransmc = temp
c
        temp     = omegaymc
        omegaymc = -omegazmc
        omegazmc = temp
c
        temp     = thetaymc
        thetaymc = -thetazmc
        thetazmc = temp
c
        temp     = yorigmc
        yorigmc  = -zorigmc
        zorigmc  = temp
c
        temp     = yorig0mc
        yorig0mc = -zorig0mc
        zorig0mc = temp
c
        temp     = dymxmc
        dymxmc   = -dzmxmc
        dzmxmc   = temp
c
        temp     = dthymxmc
        dthymxmc = -dthzmxmc
        dthzmxmc = temp
c
        temp     = ymc0
        ymc0     = -zmc0
        zmc0     = temp
c
        do nbl=1,nblock
c
          temp         = vtrans(nbl)
          vtrans(nbl)  = -wtrans(nbl)
          wtrans(nbl)  = temp
c
          temp         = dymx(nbl)
          dymx(nbl)    = -dzmx(nbl)
          dzmx(nbl)    = temp
c
          temp         = omegay(nbl)
          omegay(nbl)  = -omegaz(nbl)
          omegaz(nbl)  = temp
c
          temp         = yorig(nbl)
          yorig(nbl)   = -zorig(nbl)
          zorig(nbl)   = temp
c
          temp         = yorig0(nbl)
          yorig0(nbl)  = -zorig0(nbl)
          zorig0(nbl)  = temp
c
          temp         = thetay(nbl)
          thetay(nbl)  = -thetaz(nbl)
          thetaz(nbl)  = temp
c
          temp         = dthymx(nbl)
          dthymx(nbl)  = -dthzmx(nbl)
          dthzmx(nbl)  = temp
c
          temp         = thetayl(nbl)
          thetayl(nbl) = -thetazl(nbl)
          thetazl(nbl) = temp
c
        end do
c
c       deforming mesh movement data
c
        do nbl=1,nblock
c
          do iseg=1,maxsegdg
c
            temp               = vtrnsae(nbl,iseg)
            vtrnsae(nbl,iseg)  = -wtrnsae(nbl,iseg)
            wtrnsae(nbl,iseg)  = temp
c
            temp               = omgyae(nbl,iseg)
            omgyae(nbl,iseg)   = -omgzae(nbl,iseg)
            omgzae(nbl,iseg)   = temp
c
            temp               = yorgae(nbl,iseg)
            yorgae(nbl,iseg)   = zorgae(nbl,iseg)
            zorgae(nbl,iseg)   = temp
c
            temp               = yorgae0(nbl,iseg)
            yorgae0(nbl,iseg)  = -zorgae0(nbl,iseg)
            zorgae0(nbl,iseg)  = temp
c
            temp               = thtyae(nbl,iseg)
            thtyae(nbl,iseg)   = -thtyae(nbl,iseg)
            thtzae(nbl,iseg)   = temp
c
          end do
c
        end do
c
      end if
c
c     print summary by grids
c
      write(iunit11,1919)
 1919 format(/,17h SUMMARY BY GRIDS)
      write(iunit11,117)
  117 format(6x,4hgrid,5x,5hlevel,5x,5hblock,
     .       6x,4hjdim,6x,4hkdim,6x,4hidim,3x,9hgrid pts.)
c
      igptot = 0
      do 105 igrid=1,ngrid
      iem         = iemg(igrid)
      nbl         = nblg(igrid)
      ncg         = ncgg(igrid)
      igpts       = jdimg(nbl)*kdimg(nbl)*idimg(nbl)
      igptot      = igptot+igpts
      write(iunit11,7)igrid,levelg(nbl),nbl,
     .           jdimg(nbl),kdimg(nbl),idimg(nbl),igpts
      if (ncg.gt.0 .and. iem.eq.0) then
      do 104 n=1,ncg
      nbl         = nbl+1
      igpts       = jdimg(nbl)*kdimg(nbl)*idimg(nbl)
      igptot      = igptot+igpts
      write(iunit11,7)igrid,levelg(nbl),nbl,
     .           jdimg(nbl),kdimg(nbl),idimg(nbl),igpts
  104 continue
      end if
  105 continue
    7 format(6i10,i12)
      write(iunit11,2021)igptot
 2021 format(/,55x,5hTOTAL,i12)
c
c     print summary by levels
c
      write(iunit11,2304)
 2304 format(/,18h SUMMARY BY LEVELS)
      write(iunit11,110)
  110 format(5x,5hlevel,6x,4hgrid,5x,5hblock)
      lf = lfem
      if (lfem.eq.0) lf = lfgm
      do 115 levelc=lf,lcgm,-1
      ncell(levelc) = 0
      do 114 nbl=1,nblock
      if (levelc.ne.levelg(nbl)) go to 114
      igrid = igridg(nbl)
      ncell(levelc) = ncell(levelc) 
     .              + (jdimg(nbl)-1)*(kdimg(nbl)-1)*(idimg(nbl)-1)
      write(iunit11,7)levelg(nbl),igrid,nbl
  114 continue
  115 continue

c
      write(iunit11,8)
    8 format(/,5x,5hlevel,5x,5hcells)
      do levelc=lf,lcgm,-1
        write(iunit11,7) levelc,ncell(levelc)
      enddo
c
      if (mseq.gt.(ncgmax+1)) then
      write(iunit11,153)mseq,ncgmax
  153 format(1x,12hmseq,ncgmax=,2i4)
      write(iunit11,155)mseq,ncgmax
  155 format(1x,28herror in input, mseq, ncgmax,2i5)
      call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
      end if
c
      if (lfem.ne.0) write(iunit11,4009)lfem
 4009 format(/,31h level of finest embedded mesh=,i3)
      write(iunit11,4019)lcgm
 4019 format(/,31h level of coarsest global mesh=,i3)
      write(iunit11,4029)lfgm
 4029 format(/,31h level of finest global mesh  =,i3)
c
c     print summary of grid sequences
c
      write(iunit11,7227)
 7227 format(/,26h SUMMARY OF GRID SEQUENCES)
      write(iunit11,9)
    9 format(2x,8hsequence,3x,14hstarting level,5x,12hending level,
     .       2x,28hcells(finest global+embeded)) 
      do 19 m=1,mseq
      if (nemgl(m).eq.0) then
         ncellseq = ncell(levelt(m))
      else
         ncellseq = ncell(levelt(m)-nemgl(m)) + ncell(levelt(m))
      end if
      write(iunit11,29)m,levelt(m),levelb(m),ncellseq
   19 continue
   29 format(i10,3i17)
c
      if (abs(nbli).gt.0) then
         nblict = abs(nbli)
         do 5000 n=1,abs(nbli)
c
c    check to see if the block values are valid
c
         if (nblk(1,n).lt.1 .or. nblk(1,n).gt.ngrid) then
            write(iunit11,9350) '1',n
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
         if (nblk(2,n).lt.1 .or. nblk(2,n).gt.ngrid) then
            write(iunit11,9350) '2',n
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
 9350    format(13h side number ,a1,16h is out of range,
     .          14h for 1:1 plane,i3)
c
c    check to see if the isva values are valid for block 1
c
         if (isva(1,1,n).lt.1 .or. isva(1,1,n).gt.3) then
            write(iunit11,9400) '1',n
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
 9400    format(15h ind1 for side ,a1,16h is out of range,
     .          14h for 1:1 plane,i3)
         if (isva(1,2,n).lt.1 .or. isva(1,2,n).gt.3) then
            write(iunit11,9450) '1',n
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
 9450    format(15h ind2 for side ,a1,16h is out of range)
         if (isva(1,1,n).eq.isva(1,2,n)) then
            write(iunit11,9500) '1',n
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
 9500    format(23h ind1 = ind2 for block ,a1,
     .          14h for 1:1 plane,i3)
c
c    check to see if the isva values are valid for block 2
c
         if (isva(2,1,n).lt.1 .or. isva(2,1,n).gt.3) then
            write(iunit11,9400) '2',n
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
         if (isva(2,2,n).lt.1 .or. isva(2,2,n).gt.3) then
            write(iunit11,9450) '2',n
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
         if (isva(2,1,n).eq.isva(2,2,n)) then
            write(iunit11,9500) '2',n
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
c
c    check to see if the limblk values are valid
c
         do 9025 ib = 1, 2
         ntblk = nblg(nblk(ib,n)) 
         iflag = 0
         if (limblk(ib,1,n).lt.1 .or. 
     .       limblk(ib,1,n).gt.idimg(ntblk)) iflag = 1
         if (limblk(ib,2,n).lt.1 .or. 
     .       limblk(ib,2,n).gt.jdimg(ntblk)) iflag = 1
         if (limblk(ib,3,n).lt.1 .or. 
     .       limblk(ib,3,n).gt.kdimg(ntblk)) iflag = 1
         if (limblk(ib,4,n).lt.1 .or. 
     .       limblk(ib,4,n).gt.idimg(ntblk)) iflag = 1
         if (limblk(ib,5,n).lt.1 .or. 
     .       limblk(ib,5,n).gt.jdimg(ntblk)) iflag = 1
         if (limblk(ib,6,n).lt.1 .or. 
     .       limblk(ib,6,n).gt.kdimg(ntblk)) iflag = 1
         if (iflag.eq.1) then
            write(iunit11,9550) ib,n
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
 9550    format(29h limits out of range for side,i4,
     .          14h for 1:1 plane,i3)
 9025    continue
         if (isva(1,1,n).eq.1) i1b1 = 1
         if (isva(1,1,n).eq.2) i1b1 = 2
         if (isva(1,1,n).eq.3) i1b1 = 3
         if (isva(1,2,n).eq.1) i2b1 = 1
         if (isva(1,2,n).eq.2) i2b1 = 2
         if (isva(1,2,n).eq.3) i2b1 = 3
         if (isva(2,1,n).eq.1) i1b2 = 1
         if (isva(2,1,n).eq.2) i1b2 = 2
         if (isva(2,1,n).eq.3) i1b2 = 3
         if (isva(2,2,n).eq.1) i2b2 = 1
         if (isva(2,2,n).eq.2) i2b2 = 2
         if (isva(2,2,n).eq.3) i2b2 = 3
         idif1 = abs(limblk(1,i1b1,n) - limblk(1,i1b1+3,n))
         idif2 = abs(limblk(1,i2b1,n) - limblk(1,i2b1+3,n))
         idif3 = abs(limblk(2,i1b2,n) - limblk(2,i1b2+3,n))
         idif4 = abs(limblk(2,i2b2,n) - limblk(2,i2b2+3,n))
c
c    check to see if there is one to one communication between blocks
c
         if (idif1.ne.idif3) then
            write(iunit11,9600) '1',n
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
 9600    format(18h the limits of ind,a1,
     .          32h are not the same for both sides,
     .          14h for 1:1 plane,i3)
         if (idif2.ne.idif4) then
            write(iunit11,9600) '2',n
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
         n1   = nblk(1,n)
         n2   = nblg(nblk(1,n))  
c        if (i1b1+i2b1.eq.3) write(15,9650) '1','  k',n,n1,n2
c        if (i1b1+i2b1.eq.4) write(15,9650) '1','  j',n,n1,n2
c        if (i1b1+i2b1.eq.5) write(15,9650) '1','n i',n,n1,n2
         n1   = nblk(2,n)
         n2   = nblg(nblk(2,n))  
c        if (i1b2+i2b2.eq.3) write(15,9650) '2','  k',n,n1,n2
c        if (i1b2+i2b2.eq.4) write(15,9650) '2','  j',n,n1,n2
c        if (i1b2+i2b2.eq.5) write(15,9650) '2','n i',n,n1,n2
 9650    format(6h side ,a1,5h is a,a3,17h = constant plane,
     .          18h for 1:1 interface,i3,13h (grid/block=,i3,1h/,i3,1h))
         if (ncgg(nblk(1,n)) .ne. ncgg(nblk(2,n))) then
            write(iunit11,9700)
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
 9700    format(51h both sides are not at the same level of coarseness)
c
c    add additional 1:1 interfaces due to coarse grid definition
c
         if (ncgg(nblk(1,n)) .gt. 0) then
c           write(15,9750) ncgg(nblk(1,n))
c           write(15,9800)
            do 9050 m = 1, ncgg(nblk(1,n))
            nblict = nblict + 1
            nblon(nblict) = nblon(n)
            do 9040 ib = 1, 2
            nblk(ib,nblict) = nblg(nblk(ib,n)) + m
            lpoint = n
            if (m.gt.1) lpoint = nblict -1
            do 9030 ind = 1, 2
            isva(ib,ind,nblict) = isva(ib,ind,n)
 9030       continue
            do 9035 l = 1, 6
            itest1 = (limblk(ib,l,lpoint) - 1)/2 + 1
            itest2 = limblk(ib,l,lpoint)/2 + 1
            if (itest1 .ne. itest2 .and. itest1.ne.1) then
              write(iunit11,9850) ib 
              write(iunit11,*) itest1,itest2
              write(iunit11,*) limblk(ib,l,lpoint)
               call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
            end if
            if (itest1.eq.1) itest1 = limblk(ib,l,lpoint)
            limblk(ib,l,nblict) = itest1
 9035       continue
 9040       continue
c           write(15,9300) (nblk(ibl,nblict),
c    .                     (limblk(ibl,l,nblict),l=1,6),
c    .                     (isva(ibl,ind,nblict),ind=1,2),ibl=1,2)
 9050       continue
         end if
 9750    format(13h **** adding ,i1,26h block interfaces to allow,
     .          30h coarse grid 1:1 blocking ****)
 9800    format(39h blk1 ist jst kst ind jnd knd ind1 ind2,
     .          39h blk2 ist jst kst ind jnd knd ind1 ind2)
 9850    format(53h one of the points for blocking cannot be mapped to a
     .         ,23h coarser grid in block ,i1)
 9300    format(3x,i2,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,2x,i3,2x,i3,
     .          3x,i2,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,2x,i3,2x,i3)
c
         nblk(1,n) = nblg(nblk(1,n))
         nblk(2,n) = nblg(nblk(2,n))
 9080    continue
       	 if (nbli.gt.0) nbli =  nblict
         if (nbli.lt.0) nbli = -nblict
c
 5000    continue
      end if
c
c      write summary of 1:1 blocking indices
c
c     write(15,*) '0 Summary of 1:1 blocking indices below'
c     write(15,9800)
      do 9913 n=1,abs(nbli)
c
c    resetting limits to correspond to cell center locations
c
         do 9090 ib = 1, 2
         do 9090 m = 1, 3
         if (m + isva(ib,1,n) + isva(ib,2,n) .ne. 6) then
            if (limblk(ib,m,n) .lt. limblk(ib,m+3,n)) then
               limblk(ib,m+3,n) = limblk(ib,m+3,n) - 1
            else
               limblk(ib,m,n) = limblk(ib,m,n) - 1
            end if
         end if
 9090    continue
c
c        write(15,9300) (nblk(ibl,n),
c    .                  (limblk(ibl,l,n),l=1,6),
c    .                  (isva(ibl,ind,n),ind=1,2),ibl=1,2)
c
 9913 continue
c
c     determine more array size requirements
c
      if (icall.eq.0) then
         lmaxgr   = max(ngrid,1)
         lmaxbl   = max(nblock,1)
         lmxseg   = max(msegment,1)
         lnplts   = max(nplot3d,nprint,1)
         lmaxcs   = max(ncs2+1,1)
         lmxbli   = max(ntest,1)
         lncycm   = max(ncyctot+ntr,1)
         nintr    = 0
         lmptch   = 1
         lintmx   = 1
         lmxxe    = 1
         lmsub1   = 1
         lnmds    = mxmds
         lmaxaes  = mxaes
         lmxsegdg = mxdefseg
         lnmaster = 21
      end if
c
      return
      end
 
       subroutine blocking_skip(nskip1,nskip2,nskipt1,nskipt2)
c
c      This subroutine checks if the iskip,jskip,kskip
c      values across 1-1 block interfaces is the same.
c      If not, the skip values are set at a common divisor
c      of the two values.
c
       if(nskip1.gt.nskip2) then
        if(real(nskip1)/real(nskip2).eq.real(nskip1/nskip2)) then
         nskipt1 = nskip2
         nskipt2 = nskip2
        else
         do i = 1,8
           nskip2 = nskip2/2
           if(real(nskip1)/real(nskip2).eq.real(nskip1/nskip2)) then
            nskipt1 = nskip2
            nskipt2 = nskip2
            return
           end if
         enddo
        end if
       else
        if(real(nskip2)/real(nskip1).eq.real(nskip2/nskip1)) then
         nskipt1 = nskip1
         nskipt2 = nskip1
        else
         do i = 1,8
           nskip1 = nskip1/2
           if(real(nskip2)/real(nskip1).eq.real(nskip2/nskip1)) then
            nskipt1 = nskip1
            nskipt2 = nskip1
            return
           end if
         enddo
        end if
       end if
 
       return
       end
       subroutine skordr(nbl,nskp,i1,ijkskip,maxbl)
       dimension ijkskip(maxbl,500)
       do ii1 = nskp,1,-1
         if(ijkskip(nbl,ii1).eq.i1) return
       enddo
       nskp = nskp + 1 
       ijkskip(nbl,nskp) = ijkskip(nbl,nskp-1) 
       do ii1 = nskp-1,2,-1
         if(ijkskip(nbl,ii1-1).lt.i1) then
           ijkskip(nbl,ii1) = i1
           return 
         else
           ijkskip(nbl,ii1) = ijkskip(nbl,ii1-1) 
         end if
       enddo
       return
       end
       subroutine bndchk(n,nbl,iseg,n1,n2,n3,maxbl,maxsegdg,maxseg,
     .                   ijkskip,nskpijk1,nbcinfo,nbc,ncs) 
       dimension nbcinfo(maxbl,maxseg,7,2),ijkskip(maxbl,500)
       dimension nskpijk1(maxbl),nbc(maxbl)
       dimension ncs(maxbl,maxsegdg)
c
       do nseg=1,nbc(nbl)
        if(nbcinfo(nbl,nseg,1,n3).eq.2004.or. 
     .     nbcinfo(nbl,nseg,1,n3).eq.1005.or. 
     .     nbcinfo(nbl,nseg,1,n3).eq.1006.or.
     .     nbcinfo(nbl,nseg,1,n3).eq.2014.or.
     .     nbcinfo(nbl,nseg,1,n3).eq.2024.or.
     .     nbcinfo(nbl,nseg,1,n3).eq.2034.or.
     .     nbcinfo(nbl,nseg,1,n3).eq.2016) then
           if(nbcinfo(nbl,nseg,n1,n3).lt.ncs(nbl,iseg).and.      
     .        nbcinfo(nbl,nseg,n2,n3).ge.ncs(nbl,iseg)) 
     .        call skordr(nbl,nskpijk1(n),ncs(nbl,iseg)-1
     .                    ,ijkskip,maxbl) 
           if(nbcinfo(nbl,nseg,n1,n3).le.ncs(nbl,iseg).and.      
     .        nbcinfo(nbl,nseg,n2,n3).gt.ncs(nbl,iseg)) 
     .        call skordr(nbl,nskpijk1(n),ncs(nbl,iseg)+1
     .                    ,ijkskip,maxbl) 
        end if  
       enddo
       return
       end
       subroutine cpadd(nbl1,nbl2,nblk1,nblk2,lmblk11,lmblk12,
     .                  lmblk21,lmblk22,nskp1,ijkskip1,nskp2,
     .                  ijkskip2,maxbl,ijkdim1,ijkdim2,icrd1,icrd2) 
       dimension ijkskip1(500),ijkskip2(500)
       dimension ijkskpt1(maxbl,500),ijkskpt2(maxbl,500)
c
       lmskp11 = 0
       lmskp12 = 0
       lmskp21 = 0
       lmskp22 = 0
       ijkskpt1= 0
       ijkskpt2= 0
       do n = 1,nskp1
         ijkskpt1(nbl1,n) = ijkskip1(n)
         if(ijkskpt1(nbl1,n).eq.lmblk11) lmskp11 = n
         if(ijkskpt1(nbl1,n).eq.lmblk12) lmskp12 = n
       end do
       do n = 1,nskp2
         ijkskpt2(nbl2,n) = ijkskip2(n)
         if(ijkskpt2(nbl2,n).eq.lmblk21) lmskp21 = n
         if(ijkskpt2(nbl2,n).eq.lmblk22) lmskp22 = n
       end do
       if((nblk1.eq.nblk2).and.(icrd1.eq.icrd2)) then
         if(lmblk11.lt.lmblk12) then
           if(lmblk21.lt.lmblk22) then
             do n1 = lmskp11,lmskp12   
               itst1 = ijkskip1(n1) - lmblk11 + lmblk21
               do n2 = lmskp21,lmskp22   
                 if(ijkskip2(n2).eq.itst1) goto 100
               enddo
               call skordr(nbl2,nskp2,itst1,ijkskpt2,maxbl)
100            continue
             enddo
             do n2 = lmskp21,lmskp22   
               itst1 = ijkskip2(n2) - lmblk21 + lmblk11
               do n1 = lmskp11,lmskp12   
                 if(itst1.eq.ijkskip1(n1)) goto 200
               enddo
               call skordr(nbl2,nskp2,itst1,ijkskpt2,maxbl)
200            continue
             enddo
           else
             do n1 = lmskp11,lmskp12   
               itst1 = lmblk12 - ijkskip1(n1) + lmblk22
               do n2 = lmskp22,lmskp21   
                 if(ijkskip2(n2).eq.itst1) goto 300
               enddo
               call skordr(nbl2,nskp2,itst1,ijkskpt2,maxbl)
300            continue
             enddo
             do n2 = lmskp22,lmskp21   
               itst1 = lmblk21 - ijkskip2(n2) + lmblk11  
               do n1 = lmskp11,lmskp12   
                 if(itst1.eq.ijkskip1(n1)) goto 400
               enddo
               call skordr(nbl2,nskp2,itst1,ijkskpt2,maxbl)
400            continue
             enddo
           end if 
         else
           if(lmblk21.lt.lmblk22) then
             do n1 = lmskp12,lmskp11   
               itst1 = lmblk11 - ijkskip1(n1) + lmblk21
               do n2 = lmskp21,lmskp22   
                 if(ijkskip2(n2).eq.itst1) goto 500
               enddo
               call skordr(nbl2,nskp2,itst1,ijkskpt2,maxbl)
500            continue
             enddo
             do n2 = lmskp21,lmskp22   
               itst1 = lmblk22 - ijkskip2(n2) + lmblk12
               do n1 = lmskp12,lmskp11   
                 if(itst1.eq.ijkskip1(n1)) goto 600
               enddo
               call skordr(nbl2,nskp2,itst1,ijkskpt2,maxbl)
600            continue
             enddo 
           else
             do n1 = lmskp12,lmskp11   
               itst1 = ijkskip1(n1) - lmblk12 + lmblk22
               do n2 = lmskp22,lmskp21   
                 if(ijkskip2(n2).eq.itst1) goto 700
               enddo
               call skordr(nbl2,nskp2,itst1,ijkskpt2,maxbl)
700            continue
             enddo
             do n2 = lmskp22,lmskp21   
               itst1 = ijkskip2(n2) - lmblk22 + lmblk12
               do n1 = lmskp12,lmskp11   
                 if(itst1.eq.ijkskip1(n1)) goto 800
               enddo
               call skordr(nbl2,nskp2,itst1,ijkskpt2,maxbl)
800            continue
             enddo
           end if 
         end if
         nskp1 = nskp2
         do n = 1,nskp1
           ijkskip1(n) = ijkskpt2(nbl2,n)
         enddo 
         do n = 1,nskp2
           ijkskip2(n) = ijkskpt2(nbl2,n)
         enddo 
       else 
         if(lmblk11.lt.lmblk12) then
           if(lmblk21.lt.lmblk22) then
             do n1 = lmskp11,lmskp12   
               itst1 = ijkskip1(n1) - lmblk11 + lmblk21
               do n2 = lmskp21,lmskp22   
                 if(ijkskip2(n2).eq.itst1) goto 1100
               enddo
               call skordr(nbl2,nskp2,itst1,ijkskpt2,maxbl)
1100           continue
             enddo
             do n2 = lmskp21,lmskp22   
               itst1 = ijkskip2(n2) - lmblk21 + lmblk11
               do n1 = lmskp11,lmskp12   
                 if(itst1.eq.ijkskip1(n1)) goto 1200
               enddo
               call skordr(nbl1,nskp1,itst1,ijkskpt1,maxbl)
1200           continue
             enddo
           else
             do n1 = lmskp11,lmskp12   
               itst1 = lmblk12 - ijkskip1(n1) + lmblk22
               do n2 = lmskp22,lmskp21   
                 if(ijkskip2(n2).eq.itst1) goto 1300
               enddo
               call skordr(nbl2,nskp2,itst1,ijkskpt2,maxbl)
1300           continue
             enddo
             do n2 = lmskp22,lmskp21   
               itst1 = lmblk21 - ijkskip2(n2) + lmblk11  
               do n1 = lmskp11,lmskp12   
                 if(itst1.eq.ijkskip1(n1)) goto 1400
               enddo
               call skordr(nbl1,nskp1,itst1,ijkskpt1,maxbl)
1400           continue
             enddo
           end if 
         else
           if(lmblk21.lt.lmblk22) then
             do n1 = lmskp12,lmskp11   
               itst1 = lmblk11 - ijkskip1(n1) + lmblk21
               do n2 = lmskp21,lmskp22   
                 if(ijkskip2(n2).eq.itst1) goto 1500
               enddo
               call skordr(nbl2,nskp2,itst1,ijkskpt2,maxbl)
1500           continue
             enddo
             do n2 = lmskp21,lmskp22   
               itst1 = lmblk22 - ijkskip2(n2) + lmblk12
               do n1 = lmskp12,lmskp11   
                 if(itst1.eq.ijkskip1(n1)) goto 1600
               enddo
               call skordr(nbl1,nskp1,itst1,ijkskpt1,maxbl)
1600           continue
             enddo 
           else
             do n1 = lmskp12,lmskp11   
               itst1 = ijkskip1(n1) - lmblk12 + lmblk22
               do n2 = lmskp22,lmskp21   
                 if(ijkskip2(n2).eq.itst1) goto 1700
               enddo
               call skordr(nbl2,nskp2,itst1,ijkskpt2,maxbl)
1700           continue
             enddo
             do n2 = lmskp22,lmskp21   
               itst1 = ijkskip2(n2) - lmblk22 + lmblk12
               do n1 = lmskp12,lmskp11   
                 if(itst1.eq.ijkskip1(n1)) goto 1800
               enddo
               call skordr(nbl1,nskp1,itst1,ijkskpt1,maxbl)
1800           continue
             enddo
           end if 
         end if
         do n = 1,nskp1
           ijkskip1(n) = ijkskpt1(nbl1,n)
         enddo 
         do n = 1,nskp2
           ijkskip2(n) = ijkskpt2(nbl2,n)
         enddo 
       end if
       return
       end
