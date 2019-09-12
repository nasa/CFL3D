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
      subroutine trnsfr_vals(work,mwork,iwork,mworki,maxgr,maxbl,
     .                       maxseg,nbli,nblock,levelg,igridg,
     .                       iovrlp,nbci0,nbcj0,nbck0,
     .                       nbcidim,nbcjdim,nbckdim,iflimg,ifdsg,
     .                       iviscg,jdimg,kdimg,idimg,nblcg,idegg,
     .                       isg,jsg,ksg,ieg,jeg,keg,jlamlog,klamlog,
     .                       ilamlog,ilamhig,jlamhig,klamhig,mit,
     .                       iwfg,iadvance,iforce,lfgm,ihstry,ninter,
     .                       ngrid,ncgg,nblg,iemg,inewgg,idiagg,
     .                       nplot3d,nprint,ncs,ip3dsurf,mblk2nd,
     .                       utrans,vtrans,wtrans,omegax,omegay,
     .                       omegaz,xorig,yorig,zorig,dxmx,dymx,
     .                       dzmx,dthxmx,dthymx,dthzmx,thetax,thetay,
     .                       thetaz,rfreqt,rfreqr,xorig0,yorig0,zorig0,
     .                       time2,thetaxl,thetayl,thetazl,itrans,
     .                       irotat,idefrm,utrnsae,vtrnsae,wtrnsae,
     .                       omgxae,omgyae,omgzae,xorgae,yorgae,zorgae,
     .                       thtxae,thtyae,thtzae,rfrqtae,rfrqrae,
     .                       icsi,icsf,jcsi,jcsf,kcsi,kcsf,freq,gmass,
     .                       x0,gf0,damp,perturb,aesrfdat,bcvali,bcvalj,
     .                       bcvalk,ibcinfo,jbcinfo,kbcinfo,rkap0g,
     .                       bcfilei,bcfilej,bcfilek,mxbli,nblk,limblk,
     .                       nblon,isva,bcfiles,icsinfo,mxbcfil,maxcs,
     .                       nmds,maxaes,inpl3d,inpr,nplots,nou,bou,
     .                       nbuf,ibufdim,iskip,jskip,kskip,nsegdfrm,
     .                       idfrmseg,iaesurf,maxsegdg,xorgae0,yorgae0,
     .                       zorgae0,icouple,iprnsurf)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Load various data into work arrays for transfer from
c     host to nodes. For ease of modification, a special section is
c     reserved for keyword data (look for the word "keyword"); if 
c     additional keyword data is defined, make sure both the host
c     node sections are altered; the section provided is in the real
c     transfer section, so any integers keywords should first be 
c     converted to reals via "float" and then converted back to integer
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
#if defined DIST_MPI
c
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
#endif
c
      integer bcfilei,bcfilej,bcfilek
c
      character*120 bou(ibufdim,nbuf)
      character*80 bcfiles(mxbcfil)
      character*80 grid,plt3dg,plt3dq,output,residual,turbres,blomx,
     .             output2,printout,pplunge,ovrlap,patch,restrt,
     .             subres,subtur,grdmov,alphahist,errfile,preout,
     .             aeinp,aeout,sdhist,avgg,avgq
c
      dimension nou(nbuf)
      dimension work(mwork),iwork(mworki)
      dimension inpl3d(nplots,11),inpr(nplots,11)
      dimension levelg(maxbl),igridg(maxbl),mit(5,maxbl),
     .          iflimg(maxbl,3),ifdsg(maxbl,3),iviscg(maxbl,3),
     .          jdimg(maxbl),kdimg(maxbl),idimg(maxbl),idiagg(maxbl,3),
     .          nblcg(maxbl),idegg(maxbl,3),jsg(maxbl),
     .          ksg(maxbl),isg(maxbl),jeg(maxbl),keg(maxbl),ieg(maxbl),
     .          jlamlog(maxbl),klamlog(maxbl),ilamlog(maxbl),
     .          jlamhig(maxbl),klamhig(maxbl),ilamhig(maxbl),
     .          iwfg(maxbl,3),mblk2nd(maxbl),rkap0g(maxbl,3)
      dimension nbci0(maxbl),nbcidim(maxbl),nbcj0(maxbl),nbcjdim(maxbl),
     .          nbck0(maxbl),nbckdim(maxbl)
      dimension ncgg(maxgr),nblg(maxgr),iemg(maxgr),inewgg(maxgr)
      dimension iadvance(maxbl),iforce(maxbl),iovrlp(maxbl)
      dimension bcfilei(maxbl,maxseg,2),bcfilej(maxbl,maxseg,2),
     .          bcfilek(maxbl,maxseg,2),icsinfo(maxcs,9)
      dimension nblk(2,mxbli),limblk(2,6,mxbli),
     .          isva(2,2,mxbli),nblon(mxbli)
      dimension utrans(maxbl),vtrans(maxbl),wtrans(maxbl),
     .          omegax(maxbl),omegay(maxbl),omegaz(maxbl),xorig(maxbl),
     .          yorig(maxbl),zorig(maxbl),dxmx(maxbl),dymx(maxbl),
     .          dzmx(maxbl),dthxmx(maxbl),dthymx(maxbl),dthzmx(maxbl),
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
      dimension nsegdfrm(maxbl),idfrmseg(maxbl,maxsegdg),
     .          iaesurf(maxbl,maxsegdg)
      dimension freq(nmds,maxaes),gmass(nmds,maxaes),x0(2*nmds,maxaes),
     .          gf0(2*nmds,maxaes),damp(nmds,maxaes),
     .          perturb(nmds,maxaes,4)
      dimension aesrfdat(5,maxaes)
      dimension iskip(maxbl,500),jskip(maxbl,500),kskip(maxbl,500)
      dimension bcvali(maxbl,maxseg,12,2),
     .          bcvalj(maxbl,maxseg,12,2),bcvalk(maxbl,maxseg,12,2),
     .          ibcinfo(maxbl,maxseg,7,2),jbcinfo(maxbl,maxseg,7,2),
     .          kbcinfo(maxbl,maxseg,7,2)
c
      common /complx/ xmach_img,alpha_img,beta_img,reue_img,tinf_img,
     .                geom_img,surf_img,xrotrate_img,yrotrate_img,
     .                zrotrate_img
      common /cgns/ icgns,iccg,ibase,nzones,nsoluse,irind,jrind,krind
      common /des/ cdes,ides,cddes
      common /cfl/ dt0,dtold
      common /avgdata/ xnumavg,iteravg,xnumavg2,ipertavg,iclcd,isubit_r
      common /precond/ cprec,uref,avn
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /fluid2/ pr,prt,cbar
      common /konew/ ikoprod,isstdenom,pklimterm,ibeta8kzeta,i_bsl,
     .        keepambient,re_thetat0,i_wilcox06,i_wilcox06_chiw,
     .        i_turbprod_kterm,i_catris_kw,prod2d3dtrace,
     .        i_compress_correct,isstsf,i_wilcox98,i_wilcox98_chiw,
     .        isst2003
      common /axisym/ iaxi2plane,iaxi2planeturb,istrongturbdis,iforcev0
      common /singular/ atol
      common /deformz/ beta1,beta2,alpha1,alpha2,isktyp,negvol,meshdef,
     .                 nsprgit,ndgrd,ndwrt 
      common /mydist2/ nnodes,myhost,myid,mycomm
      common /alphait/ ialphit,cltarg,rlxalph,dalim,dalpha,icycupdt
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
      common /reyue/ reue,tinf,ivisc(3)
      common /sminn/ isminc,ismincforce
      common /twod/ i2d
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /moov/movie,nframes,icall1,lhdr,icoarsemovie,i2dmovie
      common /ginfo2/ lq2avg,iskip_blocks,inc_2d(3),inc_coarse(3) 
      common /wrbl/ nwrest
      common /wrestq/ irest,irest2
      common /igrdtyp/ ip3dgrd,ialph
      common /cpurate/ rate(5),ratesub(5),ncell(20)
      common /motionmc/ xmc0,ymc0,zmc0,utransmc,vtransmc,wtransmc,
     .                  omegaxmc,omegaymc,omegazmc,xorigmc,yorigmc,
     .                  zorigmc,xorig0mc,yorig0mc,zorig0mc,thetaxmc,
     .                  thetaymc,thetazmc,dxmxmc,dymxmc,dzmxmc,
     .                  dthxmxmc,dthymxmc,dthzmxmc,rfreqtmc,
     .                  rfreqrmc,itransmc,irotatmc,time2mc
      common /filenam/ grid,plt3dg,plt3dq,output,residual,turbres,blomx,
     .                 output2,printout,pplunge,ovrlap,patch,restrt,
     .                 subres,subtur,grdmov,alphahist,errfile,preout,
     .                 aeinp,aeout,sdhist,avgg,avgq
      common /elastic/ ndefrm,naesrf
      common /turbconv/ cflturb(7),edvislim,iturbprod,nsubturb,nfreeze,
     .                  iwarneddy,itime2read,itaturb,tur1cut,tur2cut,
     .                  iturbord,tur1cutlev,tur2cutlev
      common /entfix/ epsa_l,epsa_r
      common /elastic_ss/ idef_ss
      common /bin/ ibin,iblnk,iblnkfr,ip3dgrad
      common /memory/ memadd,memaddi
      common /is_blockbc/ is_blk(5),ie_blk(5),ivolint
      common /ghost/ irghost,iwghost
      common /noninertial/ xcentrot,ycentrot,zcentrot,xrotrate,
     .                     yrotrate,zrotrate,noninflag
      common /key/ nkey
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .        tur10(7)
      common /gridtrans/ roll_angle
      common /fullns/ ifullns
      common /curvat/ isarc2d,sarccr3,ieasmcc2d,isstrc,sstrc_crc,
     .        isar,crot,isarc3d
      common /rigidbody/ irigb,irbtrim
      common /rbstmt2/ tmass,yinert,uinfrb,qinfrb,greflrb,gaccel,crefrb,
     .                 xtmref,areat
      common /trim/ dmtrmn,dmtrmnm,dlcln,dlclnm,trtol,cmy,cnw,alf0,
     .              alf1,dzdt,thtd0,thtd1,zrg0,zrg1,dtrmsmx,dtrmsmn,
     .              dalfmx,ddtmx,ddtrm0,ddtrm1,itrmt,itrminc,fp(4,4),
     .              tp(4,4),zlfct,epstr,relax,ittrst
      common /lesinfo/ les_model,les_wallscale,cs_smagorinsky,cs_wale,
     .                 cs_vreman
      common /sourceterm/ xdir_only_source
      common /random_input/ randomize
      common /mms/ iexact_trunc,iexact_disc,iexact_ring
      common /initfac/ scal_ic
      common /plot3dtyp/ ifunct
      common /memry/ lowmem_ux
      common /constit/ i_nonlin,c_nonlin,snonlin_lim,i_tauijs,i_qcr2000,
     .                 i_qcr2013,i_qcr2013v
      common /easmlim/ cmulim
      common /easmv/ c10,c11,c2,c3,c4,c5,sigk1,cmuc1,ieasm_type
      common /is_patch/ is_pat(5),ie_pat(5),ipatch1st
      common /reystressmodel/ issglrrw2012,i_sas_rsm,i_yapterm
      common /writestuff/ ifort50write, j_ifort50write, i_ifort50write
      common /sa_options/ i_saneg,i_sanoft2,sa_cw2,sa_cw3,sa_cv1,sa_ct3,
     .                    sa_ct4,sa_cb1,sa_cb2,sa_sigma,sa_karman
      common /specialtop_kmax1001/ i_specialtop_kmax1001,
     .        a_specialtop_kmax1001,xc_specialtop_kmax1001,
     .        sig_specialtop_kmax1001,vtp_specialtop_kmax1001,
     .        wc_specialtop_kmax1001,fac_specialtop_kmax1001,
     .        cc_specialtop_kmax1001,xerf_specialtop_kmax1001,
     .        sigerf_specialtop_kmax1001
      common /iupdate/ iupdatemean
      common /wallvel_2034/ uub_2034,vvb_2034,wwb_2034,vvb_xaxisrot2034
#if defined DIST_MPI
c
c***********************
c     transfer real data
c***********************
c
      if (myid.eq.myhost) then
         do iii=1,maxbl
            work(iii)          = utrans(iii)
            work(iii+   maxbl) = vtrans(iii)
            work(iii+ 2*maxbl) = wtrans(iii)
            work(iii+ 3*maxbl) = omegax(iii)
            work(iii+ 4*maxbl) = omegay(iii)
            work(iii+ 5*maxbl) = omegaz(iii)
            work(iii+ 6*maxbl) = xorig(iii)
            work(iii+ 7*maxbl) = yorig(iii)
            work(iii+ 8*maxbl) = zorig(iii)
            work(iii+ 9*maxbl) = dxmx(iii)
            work(iii+10*maxbl) = dymx(iii)
            work(iii+11*maxbl) = dzmx(iii)
            work(iii+12*maxbl) = dthxmx(iii)
            work(iii+13*maxbl) = dthymx(iii)
            work(iii+14*maxbl) = dthzmx(iii)
            work(iii+15*maxbl) = thetax(iii)
            work(iii+16*maxbl) = thetay(iii)
            work(iii+17*maxbl) = thetaz(iii)
            work(iii+18*maxbl) = rfreqt(iii)
            work(iii+19*maxbl) = rfreqr(iii)
            work(iii+20*maxbl) = xorig0(iii)
            work(iii+21*maxbl) = yorig0(iii)
            work(iii+22*maxbl) = zorig0(iii)
            work(iii+23*maxbl) = time2(iii)
            work(iii+24*maxbl) = thetaxl(iii)
            work(iii+25*maxbl) = thetayl(iii)
            work(iii+26*maxbl) = thetazl(iii)
            work(iii+27*maxbl) = rkap0g(iii,1)
            work(iii+28*maxbl) = rkap0g(iii,2)
            work(iii+29*maxbl) = rkap0g(iii,3)
         end do
         nlast = 30*maxbl
         iii   = nlast 
         do jj=1,maxbl
            do kk=1,maxsegdg
               iii                         = iii + 1
               work(iii)                   = utrnsae(jj,kk)
               work(iii+   maxbl*maxsegdg) = vtrnsae(jj,kk)
               work(iii+ 2*maxbl*maxsegdg) = wtrnsae(jj,kk)
               work(iii+ 3*maxbl*maxsegdg) = omgxae(jj,kk)
               work(iii+ 4*maxbl*maxsegdg) = omgyae(jj,kk)
               work(iii+ 5*maxbl*maxsegdg) = omgzae(jj,kk)
               work(iii+ 6*maxbl*maxsegdg) = xorgae(jj,kk)
               work(iii+ 7*maxbl*maxsegdg) = yorgae(jj,kk)
               work(iii+ 8*maxbl*maxsegdg) = zorgae(jj,kk)
               work(iii+ 9*maxbl*maxsegdg) = thtxae(jj,kk)
               work(iii+10*maxbl*maxsegdg) = thtyae(jj,kk)
               work(iii+11*maxbl*maxsegdg) = thtzae(jj,kk)
               work(iii+12*maxbl*maxsegdg) = rfrqtae(jj,kk)
               work(iii+13*maxbl*maxsegdg) = rfrqrae(jj,kk)
               work(iii+14*maxbl*maxsegdg) = xorgae0(jj,kk)
               work(iii+15*maxbl*maxsegdg) = yorgae0(jj,kk)
               work(iii+16*maxbl*maxsegdg) = zorgae0(jj,kk)
            end do
         end do
         nlast = nlast + 17*maxbl*maxsegdg
         work(nlast+ 1) = xmc0
         work(nlast+ 2) = ymc0
         work(nlast+ 3) = zmc0
         work(nlast+ 4) = utransmc
         work(nlast+ 5) = vtransmc
         work(nlast+ 6) = wtransmc
         work(nlast+ 7) = omegaxmc
         work(nlast+ 8) = omegaymc
         work(nlast+ 9) = omegazmc
         work(nlast+10) = xorigmc
         work(nlast+11) = yorigmc
         work(nlast+12) = zorigmc
         work(nlast+13) = xorig0mc
         work(nlast+14) = yorig0mc
         work(nlast+15) = zorig0mc
         work(nlast+16) = thetaxmc
         work(nlast+17) = thetaymc
         work(nlast+18) = thetazmc
         work(nlast+19) = dxmxmc
         work(nlast+20) = dymxmc
         work(nlast+21) = dzmxmc
         work(nlast+22) = dthxmxmc
         work(nlast+23) = dthymxmc
         work(nlast+24) = dthzmxmc
         work(nlast+25) = rfreqtmc
         work(nlast+26) = rfreqrmc
         work(nlast+27) = time2mc
         work(nlast+28) = xmach
         work(nlast+29) = alpha
         work(nlast+30) = beta
         work(nlast+31) = reue
         work(nlast+32) = tinf
         work(nlast+33) = sref
         work(nlast+34) = cref
         work(nlast+35) = bref
         work(nlast+36) = xmc
         work(nlast+37) = ymc
         work(nlast+38) = zmc
         work(nlast+39) = dt
         work(nlast+40) = fmax
         work(nlast+41) = cfltau
         work(nlast+42) = dtold
         do j = 1,3
            work(nlast+42+j) = epsssc(j)
            work(nlast+45+j) = epsssr(j)
            work(nlast+48+j) = rkap0(j)
         end do
         work(nlast+52) = beta1
         work(nlast+53) = alpha1
         work(nlast+54) = beta2
         work(nlast+55) = alpha2
         nlast = nlast + 55
         iii   = nlast
         do kk=1,maxaes
            do jj=1,nmds
               iii = iii + 1
               work(iii)               = freq(jj,kk)
               work(iii+  nmds*maxaes) = gmass(jj,kk)
               work(iii+2*nmds*maxaes) = damp(jj,kk)
            end do
         end do
         nlast = nlast + 3*nmds*maxaes
         iii   = nlast
         do kk=1,maxaes
            do jj=1,nmds
               do ll=1,4
                  iii = iii + 1
                  work(iii) = perturb(jj,kk,ll)
               end do
            end do
         end do
         nlast = nlast + 4*nmds*maxaes
         iii   = nlast
         do kk=1,maxaes
            do jj=1,2*nmds
               iii = iii + 1
               work(iii)               = x0(jj,kk)
               work(iii+2*nmds*maxaes) = gf0(jj,kk)
            end do
         end do
         nlast = nlast + 4*nmds*maxaes
         iii   = nlast 
         do kk=1,maxaes
            do jj=1,5
               iii = iii + 1
               work(iii) = aesrfdat(jj,kk)
            end do
         end do
         nlast = nlast + 5*maxaes
         iii   = nlast
         do nbl=1,maxbl
            do nseg=1,maxseg
               do mm=1,12
                  do ll=1,2
                     iii = iii + 1
                     work(iii) = bcvali(nbl,nseg,mm,ll)
                     work(iii+maxbl*maxseg*24) = bcvalj(nbl,nseg,mm,ll)
                     work(iii+maxbl*maxseg*48) = bcvalk(nbl,nseg,mm,ll)
                  end do
               end do
            end do
         end do
         nlast = nlast + maxbl*maxseg*72
c
c        keyword data
c
         work(nlast+1)  = gamma
         work(nlast+2)  = pr
         work(nlast+3)  = prt
         work(nlast+4)  = cbar
         work(nlast+5)  = atol
         work(nlast+6)  = xmach_img
         work(nlast+7)  = alpha_img
         work(nlast+8)  = beta_img
         work(nlast+9)  = reue_img
         work(nlast+10) = tinf_img
         work(nlast+11) = geom_img
         work(nlast+12) = icgns
         work(nlast+13) = cprec
         work(nlast+14) = uref
         work(nlast+15) = avn
         work(nlast+16) = cltarg
         work(nlast+17) = rlxalph
         work(nlast+18) = nsubturb
         work(nlast+19) = cflturb(1)
         work(nlast+20) = cflturb(2)
         work(nlast+21) = cflturb(3)
         work(nlast+22) = cflturb(4)
         work(nlast+23) = cflturb(5)
         work(nlast+24) = cflturb(6)
         work(nlast+25) = cflturb(7)
         work(nlast+26) = epsa_r
         work(nlast+27) = nfreeze
         work(nlast+28) = ivolint
         work(nlast+29) = idef_ss
         work(nlast+30) = ibin
         work(nlast+31) = iblnk
         work(nlast+32) = surf_img
         work(nlast+33) = memadd
         work(nlast+34) = memaddi
         work(nlast+35) = negvol
         work(nlast+36) = iblnkfr
         work(nlast+37) = edvislim
         work(nlast+38) = iturbprod
         work(nlast+39) = ip3dgrad
         work(nlast+40) = irghost
         work(nlast+41) = iwghost
         work(nlast+42) = dalim
         work(nlast+43) = icycupdt
         work(nlast+44) = xcentrot
         work(nlast+45) = ycentrot
         work(nlast+46) = zcentrot
         work(nlast+47) = xrotrate
         work(nlast+48) = yrotrate
         work(nlast+49) = zrotrate
         work(nlast+50) = noninflag
         work(nlast+51) = xrotrate_img
         work(nlast+52) = yrotrate_img
         work(nlast+53) = zrotrate_img
         work(nlast+54) = itime2read
         work(nlast+55) = itaturb
         work(nlast+56) = ides
         work(nlast+57) = cdes
         work(nlast+58) = iteravg
         work(nlast+59) = tur10(1)
         work(nlast+60) = tur10(2)
         work(nlast+61) = tur10(3)
         work(nlast+62) = tur10(4)
         work(nlast+63) = tur10(5)
         work(nlast+64) = tur10(6)
         work(nlast+65) = tur10(7)
         work(nlast+66) = tur1cut
         work(nlast+67) = roll_angle
         work(nlast+68) = ikoprod
         work(nlast+69) = isstdenom
         work(nlast+70) = pklimterm
         work(nlast+71) = iaxi2plane
         work(nlast+72) = iturbord
         work(nlast+73) = tur2cut
         work(nlast+74) = tur1cutlev
         work(nlast+75) = tur2cutlev
         work(nlast+76) = meshdef
         work(nlast+77) = ifullns
         work(nlast+78) = ibeta8kzeta
         work(nlast+79) = isarc2d
         work(nlast+80) = sarccr3
         work(nlast+81) = ieasmcc2d
         work(nlast+82) = ipertavg
         work(nlast+83) = icoarsemovie
         work(nlast+84) = i2dmovie
         work(nlast+85) = iclcd
         if(iskip_blocks .eq. 0) iskip_blocks = 1
         work(nlast+86) = iskip_blocks
         work(nlast+87) = cfltauMax
         work(nlast+88) = cfltau0
         work(nlast+89) = irbtrim
         work(nlast+90) = irigb
         work(nlast+91) = greflrb
         work(nlast+92) = tmass
         work(nlast+93) = yinert
         work(nlast+94) = gaccel
         work(nlast+95) = relax
         work(nlast+96) = itrminc
         work(nlast+97) = tp(1,1)
         work(nlast+98) = tp(2,1)
         work(nlast+99) = tp(1,2)
         work(nlast+100) = tp(2,2)
         work(nlast+101) = ndgrd
         work(nlast+102) = ndwrt
         work(nlast+103) = i_bsl
         work(nlast+104) = keepambient
         work(nlast+105) = re_thetat0
         work(nlast+106)= cs_smagorinsky
         work(nlast+107)= xdir_only_source
         work(nlast+108)= randomize
         work(nlast+109)= iexact_trunc
         work(nlast+110)= iexact_disc
         work(nlast+111)= iexact_ring
         work(nlast+112)= i_wilcox06
         work(nlast+113)= i_wilcox06_chiw
         work(nlast+114)= i_turbprod_kterm
         work(nlast+115)= i_catris_kw
         work(nlast+116)= ismincforce
         work(nlast+117)= prod2d3dtrace
         work(nlast+118)= i_compress_correct
         work(nlast+119)= les_model
         work(nlast+120)= les_wallscale
         work(nlast+121)= cs_wale
         work(nlast+122)= cs_vreman
         work(nlast+123)= isstrc
         work(nlast+124)= sstrc_crc
         work(nlast+125)= isstsf
         work(nlast+126)= cddes
         work(nlast+127)= scal_ic
         work(nlast+128)= ifunct
         work(nlast+129)= lowmem_ux
         work(nlast+130)= isar
         work(nlast+131)= crot
         work(nlast+132)= i_nonlin
         work(nlast+133)= c_nonlin
         work(nlast+134)= snonlin_lim
         work(nlast+135)= isarc3d
         work(nlast+136)= isubit_r
         work(nlast+137)= i_wilcox98
         work(nlast+138)= i_wilcox98_chiw
         work(nlast+139)= cmulim
         work(nlast+140)= iaxi2planeturb
         work(nlast+141)= istrongturbdis
         work(nlast+142)= ieasm_type
         work(nlast+143)= ipatch1st
         work(nlast+144)= isst2003
         work(nlast+145)= issglrrw2012
         work(nlast+146)= ifort50write
         work(nlast+147)= j_ifort50write
         work(nlast+148)= i_ifort50write
         work(nlast+149)= i_sas_rsm
         work(nlast+150)= iforcev0
         work(nlast+151)= i_saneg
         work(nlast+152)= i_sanoft2
         work(nlast+153)= i_lam_forcezero
         work(nlast+154)= i_specialtop_kmax1001
         work(nlast+155)= a_specialtop_kmax1001
         work(nlast+156)= xc_specialtop_kmax1001
         work(nlast+157)= sig_specialtop_kmax1001
         work(nlast+158)= vtp_specialtop_kmax1001
         work(nlast+159)= wc_specialtop_kmax1001
         work(nlast+160)= fac_specialtop_kmax1001
         work(nlast+161)= cc_specialtop_kmax1001
         work(nlast+162)= xerf_specialtop_kmax1001
         work(nlast+163)= sigerf_specialtop_kmax1001
         work(nlast+164)= sa_cw2
         work(nlast+165)= sa_cw3
         work(nlast+166)= sa_cv1
         work(nlast+167)= sa_ct3
         work(nlast+168)= sa_ct4
         work(nlast+169)= sa_cb1
         work(nlast+170)= sa_cb2
         work(nlast+171)= sa_sigma
         work(nlast+172)= sa_karman
         work(nlast+173)= iupdatemean
         work(nlast+174)= uub_2034
         work(nlast+175)= vvb_2034
         work(nlast+176)= wwb_2034
         work(nlast+177)= vvb_xaxisrot2034
         work(nlast+178)= i_tauijs
         work(nlast+179)= i_qcr2000
         work(nlast+180)= i_qcr2013
         work(nlast+181)= i_qcr2013v
         work(nlast+182)= i_yapterm
         nlast = nlast+nkey
      end if
c
      nvals = 30*maxbl + 72*maxbl*maxseg   + 55   + 11*nmds*maxaes
     .      + 5*maxaes + 17*maxbl*maxsegdg + nkey
      if (myid.eq.myhost) then
         if (nvals.ne.nlast) then
            write(11,'('' error in subroutine trnsfr_vals, real'',
     .      '' transfer: nvals.ne.nlast'')')
            write(11,'('' nvals,nlast '',i8,'' '',i8)') nvals,nlast
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)         
         end if
      end if
      call MPI_Bcast (work, nvals, MY_MPI_REAL,
     .                myhost,mycomm, ierr)
c
      if (myid.ne.myhost) then
         do iii=1,maxbl
            utrans(iii)   = work(iii)         
            vtrans(iii)   = work(iii+   maxbl)
            wtrans(iii)   = work(iii+ 2*maxbl)
            omegax(iii)   = work(iii+ 3*maxbl)
            omegay(iii)   = work(iii+ 4*maxbl)
            omegaz(iii)   = work(iii+ 5*maxbl)
            xorig(iii)    = work(iii+ 6*maxbl)
            yorig(iii)    = work(iii+ 7*maxbl)
            zorig(iii)    = work(iii+ 8*maxbl)
            dxmx(iii)     = work(iii+ 9*maxbl)
            dymx(iii)     = work(iii+10*maxbl)
            dzmx(iii)     = work(iii+11*maxbl)
            dthxmx(iii)   = work(iii+12*maxbl)
            dthymx(iii)   = work(iii+13*maxbl)
            dthzmx(iii)   = work(iii+14*maxbl)
            thetax(iii)   = work(iii+15*maxbl)
            thetay(iii)   = work(iii+16*maxbl)
            thetaz(iii)   = work(iii+17*maxbl)
            rfreqt(iii)   = work(iii+18*maxbl)
            rfreqr(iii)   = work(iii+19*maxbl)
            xorig0(iii)   = work(iii+20*maxbl)
            yorig0(iii)   = work(iii+21*maxbl)
            zorig0(iii)   = work(iii+22*maxbl)
            time2(iii)    = work(iii+23*maxbl)
            thetaxl(iii)  = work(iii+24*maxbl)
            thetayl(iii)  = work(iii+25*maxbl)
            thetazl(iii)  = work(iii+26*maxbl)
            rkap0g(iii,1) = work(iii+27*maxbl)
            rkap0g(iii,2) = work(iii+28*maxbl)
            rkap0g(iii,3) = work(iii+29*maxbl) 
         end do
         nlast = 30*maxbl
         iii   = nlast
         do jj=1,maxbl
            do kk=1,maxsegdg
               iii            = iii + 1
               utrnsae(jj,kk) = work(iii)
               vtrnsae(jj,kk) = work(iii+   maxbl*maxsegdg)
               wtrnsae(jj,kk) = work(iii+ 2*maxbl*maxsegdg)
               omgxae(jj,kk)  = work(iii+ 3*maxbl*maxsegdg)
               omgyae(jj,kk)  = work(iii+ 4*maxbl*maxsegdg)
               omgzae(jj,kk)  = work(iii+ 5*maxbl*maxsegdg)
               xorgae(jj,kk)  = work(iii+ 6*maxbl*maxsegdg)
               yorgae(jj,kk)  = work(iii+ 7*maxbl*maxsegdg)
               zorgae(jj,kk)  = work(iii+ 8*maxbl*maxsegdg)
               thtxae(jj,kk)  = work(iii+ 9*maxbl*maxsegdg)
               thtyae(jj,kk)  = work(iii+10*maxbl*maxsegdg)
               thtzae(jj,kk)  = work(iii+11*maxbl*maxsegdg)
               rfrqtae(jj,kk) = work(iii+12*maxbl*maxsegdg)
               rfrqrae(jj,kk) = work(iii+13*maxbl*maxsegdg)
               xorgae0(jj,kk) = work(iii+14*maxbl*maxsegdg)
               yorgae0(jj,kk) = work(iii+15*maxbl*maxsegdg)
               zorgae0(jj,kk) = work(iii+16*maxbl*maxsegdg)
            end do
         end do
         nlast = nlast + 17*maxbl*maxsegdg
         xmc0     = work(nlast+ 1)
         ymc0     = work(nlast+ 2)
         zmc0     = work(nlast+ 3)
         utransmc = work(nlast+ 4)
         vtransmc = work(nlast+ 5)
         wtransmc = work(nlast+ 6)
         omegaxmc = work(nlast+ 7)
         omegaymc = work(nlast+ 8)
         omegazmc = work(nlast+ 9)
         xorigmc  = work(nlast+10)
         yorigmc  = work(nlast+11)
         zorigmc  = work(nlast+12)
         xorig0mc = work(nlast+13)
         yorig0mc = work(nlast+14)
         zorig0mc = work(nlast+15)
         thetaxmc = work(nlast+16)
         thetaymc = work(nlast+17)
         thetazmc = work(nlast+18)
         dxmxmc   = work(nlast+19)
         dymxmc   = work(nlast+20)
         dzmxmc   = work(nlast+21)
         dthxmxmc = work(nlast+22)
         dthymxmc = work(nlast+23)
         dthzmxmc = work(nlast+24)
         rfreqtmc = work(nlast+25)
         rfreqrmc = work(nlast+26)
         time2mc  = work(nlast+27)
         xmach    = work(nlast+28)
         alpha    = work(nlast+29)
         beta     = work(nlast+30)
         reue     = work(nlast+31)
         tinf     = work(nlast+32)
         sref     = work(nlast+33)
         cref     = work(nlast+34)
         bref     = work(nlast+35)
         xmc      = work(nlast+36)
         ymc      = work(nlast+37)
         zmc      = work(nlast+38)
         dt       = work(nlast+39)
         fmax     = work(nlast+40)
         cfltau   = work(nlast+41)
         dtold    = work(nlast+42)
         do j = 1,3
            epsssc(j) = work(nlast+42+j)
            epsssr(j) = work(nlast+45+j)
            rkap0(j)  = work(nlast+48+j)
         end do
         beta1    = work(nlast+52)
         alpha1   = work(nlast+53)
         beta2    = work(nlast+54)
         alpha2   = work(nlast+55)
         nlast = nlast + 55
         iii   = nlast
         do kk=1,maxaes
            do jj=1,nmds
               iii = iii + 1
               freq(jj,kk)  = work(iii)
               gmass(jj,kk) = work(iii+  nmds*maxaes)
               damp(jj,kk)  = work(iii+2*nmds*maxaes)
            end do
         end do
         nlast = nlast + 3*nmds*maxaes
         iii   = nlast
         do kk=1,maxaes
            do jj=1,nmds
               do ll=1,4
                  iii = iii + 1
                  perturb(jj,kk,ll) = work(iii)
               end do
            end do
         end do
         nlast = nlast + 4*nmds*maxaes
         iii   = nlast 
         do kk=1,maxaes
            do jj=1,2*nmds
               iii = iii + 1
               x0(jj,kk)  = work(iii)
               gf0(jj,kk) = work(iii+2*nmds*maxaes)
            end do
         end do
         nlast = nlast + 4*nmds*maxaes
         iii   = nlast
         do kk=1,maxaes
            do jj=1,5
               iii = iii + 1
               aesrfdat(jj,kk) = work(iii)
            end do
         end do
         nlast = nlast + 5*maxaes
         iii   = nlast
         do nbl=1,maxbl
            do nseg=1,maxseg
               do mm=1,12
                  do ll=1,2
                     iii = iii + 1
                     bcvali(nbl,nseg,mm,ll) = work(iii)
                     bcvalj(nbl,nseg,mm,ll) = work(iii+maxbl*maxseg*24)
                     bcvalk(nbl,nseg,mm,ll) = work(iii+maxbl*maxseg*48)
                  end do
               end do
            end do
         end do
         nlast = nlast + maxbl*maxseg*72
c
c        keyword data
c
         gamma = work(nlast+1)
         pr    = work(nlast+2)
         prt   = work(nlast+3)
         cbar  = work(nlast+4)
         atol  = work(nlast+5)
         xmach_img = work(nlast+6)
         alpha_img = work(nlast+7)
         beta_img  = work(nlast+8)
         reue_img  = work(nlast+9)
         tinf_img  = work(nlast+10)
         geom_img  = work(nlast+11)
         icgns     = int(work(nlast+12))
         cprec     = work(nlast+13)
         uref      = work(nlast+14)
         avn       = work(nlast+15)
         cltarg    = work(nlast+16)
         rlxalph   = work(nlast+17)
         nsubturb  = int(work(nlast+18))
         cflturb(1)= work(nlast+19)
         cflturb(2)= work(nlast+20)
         cflturb(3)= work(nlast+21)
         cflturb(4)= work(nlast+22)
         cflturb(5)= work(nlast+23)
         cflturb(6)= work(nlast+24)
         cflturb(7)= work(nlast+25)
         epsa_r    = work(nlast+26)
         nfreeze   = int(work(nlast+27))
         ivolint   = int(work(nlast+28))
         idef_ss   = int(work(nlast+29))
         ibin      = int(work(nlast+30))
         iblnk     = int(work(nlast+31))
         surf_img  = work(nlast+32)
         memadd    = int(work(nlast+33))
         memaddi   = int(work(nlast+34))
         negvol    = int(work(nlast+35))
         iblnkfr   = int(work(nlast+36))
         edvislim  = work(nlast+37)
         iturbprod = int(work(nlast+38))
         ip3dgrad  = int(work(nlast+39))
         irghost   = int(work(nlast+40))
         iwghost   = int(work(nlast+41))
         dalim     = work(nlast+42)
         icycupdt  = int(work(nlast+43))
         xcentrot  = work(nlast+44)
         ycentrot  = work(nlast+45)
         zcentrot  = work(nlast+46)
         xrotrate  = work(nlast+47)
         yrotrate  = work(nlast+48)
         zrotrate  = work(nlast+49)
         noninflag = int(work(nlast+50))
         xrotrate_img = work(nlast+51)
         yrotrate_img = work(nlast+52)
         zrotrate_img = work(nlast+53)
         itime2read=int(work(nlast+54))
         itaturb   =int(work(nlast+55))
         ides      = int(work(nlast+56))
         cdes      = work(nlast+57)
         iteravg   = int(work(nlast+58))
         tur10(1)  = work(nlast+59)
         tur10(2)  = work(nlast+60)
         tur10(3)  = work(nlast+61)
         tur10(4)  = work(nlast+62)
         tur10(5)  = work(nlast+63)
         tur10(6)  = work(nlast+64)
         tur10(7)  = work(nlast+65)
         tur1cut   = work(nlast+66)
         roll_angle= work(nlast+67)
         ikoprod   = int(work(nlast+68))
         isstdenom = int(work(nlast+69))
         pklimterm = work(nlast+70)
         iaxi2plane= int(work(nlast+71))
         iturbord  = int(work(nlast+72))
         tur2cut   = work(nlast+73)
         tur1cutlev= work(nlast+74)
         tur2cutlev= work(nlast+75)
         meshdef   = int(work(nlast+76))
         ifullns   = int(work(nlast+77))
         ibeta8kzeta=int(work(nlast+78))
         isarc2d   = int(work(nlast+79))
         sarccr3   = work(nlast+80)
         ieasmcc2d = int(work(nlast+81))
         ipertavg  = int(work(nlast+82))
         icoarsemovie = int(work(nlast+83))
         i2dmovie = int(work(nlast+84))
         iclcd     = int(work(nlast+85))
         iskip_blocks = int(work(nlast+86))
         cfltauMax = work(nlast+87)
         cfltau0 = work(nlast+88)
         irbtrim   = int(work(nlast+89))
         irigb     = int(work(nlast+90))
         greflrb   = work(nlast+91)
         tmass     = work(nlast+92)
         yinert    = work(nlast+93)
         gaccel    = work(nlast+94)
         relax     = work(nlast+95)
         itrminc   = int(work(nlast+96))
         tp(1,1)   = work(nlast+97)
         tp(2,1)   = work(nlast+98)
         tp(1,2)   = work(nlast+99)
         tp(2,2)   = work(nlast+100)
         ndgrd     = int(work(nlast+101))
         ndwrt     = int(work(nlast+102))
         i_bsl     = int(work(nlast+103))
         keepambient = int(work(nlast+104))
         re_thetat0  = work(nlast+105)
         cs_smagorinsky=work(nlast+106)
         xdir_only_source = work(nlast+107)
         randomize= work(nlast+108)
         iexact_trunc= int(work(nlast+109))
         iexact_disc= int(work(nlast+110))
         iexact_ring= int(work(nlast+111))
         i_wilcox06= int(work(nlast+112))
         i_wilcox06_chiw= int(work(nlast+113))
         i_turbprod_kterm= int(work(nlast+114))
         i_catris_kw= int(work(nlast+115))
         ismincforce= int(work(nlast+116))
         prod2d3dtrace= work(nlast+117)
         i_compress_correct = int(work(nlast+118))
         les_model = int(work(nlast+119))
         les_wallscale = int(work(nlast+120))
         cs_wale = work(nlast+121)
         cs_vreman = work(nlast+122)
         isstrc = int(work(nlast+123))
         sstrc_crc = work(nlast+124)
         isstsf = int(work(nlast+125))
         cddes = work(nlast+126)
         scal_ic = work(nlast+127)
         ifunct = int(work(nlast+128))
         lowmem_ux = int(work(nlast+129))
         isar = int(work(nlast+130))
         crot = work(nlast+131)
         i_nonlin = int(work(nlast+132))
         c_nonlin = work(nlast+133)
         snonlin_lim = work(nlast+134)
         isarc3d = int(work(nlast+135))
         isubit_r = int(work(nlast+136))
         i_wilcox98= int(work(nlast+137))
         i_wilcox98_chiw= int(work(nlast+138))
         cmulim= work(nlast+139)
         iaxi2planeturb= int(work(nlast+140))
         istrongturbdis= int(work(nlast+141))
         ieasm_type= int(work(nlast+142))
         ipatch1st= int(work(nlast+143))
         isst2003= int(work(nlast+144))
         issglrrw2012= int(work(nlast+145))
         ifort50write= int(work(nlast+146))
         j_ifort50write= int(work(nlast+147))
         i_ifort50write= int(work(nlast+148))
         i_sas_rsm= int(work(nlast+149))
         iforcev0 = int(work(nlast+150))
         i_saneg = int(work(nlast+151))
         i_sanoft2 = int(work(nlast+152))
         i_lam_forcezero = int(work(nlast+153))
         i_specialtop_kmax1001 = int(work(nlast+154))
         a_specialtop_kmax1001 = work(nlast+155)
         xc_specialtop_kmax1001 = work(nlast+156)
         sig_specialtop_kmax1001 = work(nlast+157)
         vtp_specialtop_kmax1001 = work(nlast+158)
         wc_specialtop_kmax1001 = work(nlast+159)
         fac_specialtop_kmax1001 = work(nlast+160)
         cc_specialtop_kmax1001 = work(nlast+161)
         xerf_specialtop_kmax1001 = work(nlast+162)
         sigerf_specialtop_kmax1001 = work(nlast+163)
         sa_cw2 = work(nlast+164)
         sa_cw3 = work(nlast+165)
         sa_cv1 = work(nlast+166)
         sa_ct3 = work(nlast+167)
         sa_ct4 = work(nlast+168)
         sa_cb1 = work(nlast+169)
         sa_cb2 = work(nlast+170)
         sa_sigma = work(nlast+171)
         sa_karman = work(nlast+172)
         iupdatemean = int(work(nlast+173))
         uub_2034 = work(nlast+174)
         vvb_2034 = work(nlast+175)
         wwb_2034 = work(nlast+176)
         vvb_xaxisrot2034 = work(nlast+177)
         i_tauijs = int(work(nlast+178))
         i_qcr2000 = int(work(nlast+179))
         i_qcr2013 = int(work(nlast+180))
         i_qcr2013v = int(work(nlast+181))
         i_yapterm = int(work(nlast+182))
      end if
c
c**************************
c     transfer integer data
c**************************
c
      if (myid.eq.myhost) then
         do iii=1,maxbl
            iwork(iii)         = itrans(iii)
            iwork(iii+maxbl)   = irotat(iii)
            iwork(iii+2*maxbl) = idefrm(iii)
         end do
         nlast = 3*maxbl
         iwork(nlast+1)  = itransmc
         iwork(nlast+2)  = irotatmc
         iwork(nlast+3)  = ndefrm
         iwork(nlast+4)  = naesrf
         iwork(nlast+5)  = ihstry
         iwork(nlast+6)  = irest
         iwork(nlast+7)  = iflagts
         iwork(nlast+8)  = iunst
         iwork(nlast+9)  = nwrest
         iwork(nlast+10) = ichk
         iwork(nlast+11) = i2d
         iwork(nlast+12) = ntstep
         iwork(nlast+13) = ita
         iwork(nlast+14) = mseq
         iwork(nlast+15) = issc
         iwork(nlast+16) = issr
         iwork(nlast+17) = ncyctot
         iwork(nlast+18) = lfgm
         iwork(nlast+19) = nplot3d
         iwork(nlast+20) = movie
         iwork(nlast+21) = nprint
         iwork(nlast+22) = nbli
         iwork(nlast+23) = nblock
         iwork(nlast+24) = ngrid
         iwork(nlast+25) = ninter
         iwork(nlast+26) = mgflag
         iwork(nlast+27) = iconsf
         iwork(nlast+28) = mtt
         iwork(nlast+29) = ngam
         iwork(nlast+30) = ideg(1)
         iwork(nlast+31) = ideg(2)
         iwork(nlast+32) = ideg(3)
         iwork(nlast+33) = ifds(1)
         iwork(nlast+34) = ifds(2)
         iwork(nlast+35) = ifds(3)
         iwork(nlast+36) = ilamlo
         iwork(nlast+37) = ilamhi
         iwork(nlast+38) = jlamlo
         iwork(nlast+39) = jlamhi
         iwork(nlast+40) = klamlo
         iwork(nlast+41) = klamhi
         iwork(nlast+42) = ivmx
         iwork(nlast+43) = isminc
         iwork(nlast+44) = ip3dgrd
         iwork(nlast+45) = ialph
         iwork(nlast+46) = iipv
         iwork(nlast+47) = irest2
         iwork(nlast+48) = nnodes
         iwork(nlast+49) = ncs
         iwork(nlast+50) = ialphit
         iwork(nlast+51) = ip3dsurf
         iwork(nlast+52) = isktyp
         iwork(nlast+53) = iprnsurf
         iwork(nlast+54) = nsprgit
         nlast = nlast + 54
         do iii=1,maxgr
            iwork(nlast+iii)         = ncgg(iii)
            iwork(nlast+iii+maxgr)   = iemg(iii)
            iwork(nlast+iii+2*maxgr) = nblg(iii)
            iwork(nlast+iii+3*maxgr) = inewgg(iii)
         end do
         nlast = nlast + 4*maxgr
         do iii=1,20
            iwork(iii+nlast) = ncell(iii)
         end do
         nlast = nlast + 20
         do iii=1,maxbl
            iwork(nlast+iii)       = iadvance(iii)
            iwork(nlast+iii+maxbl) = iforce(iii)
         end do
         nlast = nlast + 2*maxbl
         iii   = nlast
         do nbl=1,maxbl
            do i = 1,3
               iii = iii + 1
               iwork(iii)         = iviscg(nbl,i)
               iwork(iii+3*maxbl) = idegg(nbl,i)
               iwork(iii+6*maxbl) = iwfg(nbl,i)
            end do
         end do
         nlast = nlast + 9*maxbl
         do iii=1,maxbl
            iwork(nlast+iii)          = idimg(iii)
            iwork(nlast+iii+maxbl)    = jdimg(iii)
            iwork(nlast+iii+2*maxbl)  = kdimg(iii)
            iwork(nlast+iii+3*maxbl)  = ilamlog(iii)
            iwork(nlast+iii+4*maxbl)  = ilamhig(iii)
            iwork(nlast+iii+5*maxbl)  = jlamlog(iii)
            iwork(nlast+iii+6*maxbl)  = jlamhig(iii)
            iwork(nlast+iii+7*maxbl)  = klamlog(iii)
            iwork(nlast+iii+8*maxbl)  = klamhig(iii)
            iwork(nlast+iii+9*maxbl)  = jsg(iii)
            iwork(nlast+iii+10*maxbl) = ksg(iii)
            iwork(nlast+iii+11*maxbl) = isg(iii)
            iwork(nlast+iii+12*maxbl) = jeg(iii)
            iwork(nlast+iii+13*maxbl) = keg(iii)
            iwork(nlast+iii+14*maxbl) = ieg(iii)
            iwork(nlast+iii+15*maxbl) = nblcg(iii)
         end do
         nlast = nlast + 16*maxbl
         iii   = nlast
         do nbl=1,maxbl
            do i = 1,3
               iii = iii + 1
               iwork(iii)         = idiagg(nbl,i)
               iwork(iii+3*maxbl) = iflimg(nbl,i)
               iwork(iii+6*maxbl) = ifdsg(nbl,i)
            end do
         end do
         nlast = nlast + 9*maxbl
         do iii=1,maxbl
            iwork(nlast+iii)         = nbci0(iii)
            iwork(nlast+iii+maxbl)   = nbcidim(iii)
            iwork(nlast+iii+2*maxbl) = nbcj0(iii)
            iwork(nlast+iii+3*maxbl) = nbcjdim(iii)
            iwork(nlast+iii+4*maxbl) = nbck0(iii)
            iwork(nlast+iii+5*maxbl) = nbckdim(iii)
            iwork(nlast+iii+6*maxbl) = iovrlp(iii)
            iwork(nlast+iii+7*maxbl) = igridg(iii)
            iwork(nlast+iii+8*maxbl) = levelg(iii)
            iwork(nlast+iii+9*maxbl) = mblk2nd(iii)
         end do
         nlast = nlast + 10*maxbl
         do iii=1,5
            iwork(nlast+iii)    = ncyc1(iii)
            iwork(nlast+iii+5)  = nitfo1(iii)
            iwork(nlast+iii+10) = levelt(iii)
            iwork(nlast+iii+15) = levelb(iii)
         end do
         nlast = nlast + 20
         iii   = nlast
         do m = 1,5
            do i = 1, 5
               iii = iii + 1
               iwork(iii) = mit(i,m)
            end do
         end do
         nlast = nlast + 25
         iii   = nlast
         do kk=1,maxsegdg
            do jj=1,maxbl
               iii = iii + 1
               iwork(iii)                  = icsi(jj,kk)
               iwork(iii+  maxbl*maxsegdg) = icsf(jj,kk)
               iwork(iii+2*maxbl*maxsegdg) = jcsi(jj,kk)
               iwork(iii+3*maxbl*maxsegdg) = jcsf(jj,kk)
               iwork(iii+4*maxbl*maxsegdg) = kcsi(jj,kk)
               iwork(iii+5*maxbl*maxsegdg) = kcsf(jj,kk)
               iwork(iii+6*maxbl*maxsegdg) = idfrmseg(jj,kk)
               iwork(iii+7*maxbl*maxsegdg) = iaesurf(jj,kk)
               iwork(iii+8*maxbl*maxsegdg) = icouple(jj,kk)
            end do
         end do
         nlast = nlast + 9*maxbl*maxsegdg
         do iii = 1,maxbl
            iwork(nlast+iii) = nsegdfrm(iii)
         end do
         nlast = nlast + maxbl
         iii   = nlast
         do nbl=1,maxbl
            do nseg=1,maxseg
               do ll=1,2
                  iii = iii + 1
                  iwork(iii)                = bcfilei(nbl,nseg,ll)
                  iwork(iii+maxbl*maxseg*2) = bcfilej(nbl,nseg,ll) 
                  iwork(iii+maxbl*maxseg*4) = bcfilek(nbl,nseg,ll)
               end do
            end do
         end do
         nlast = nlast + maxbl*maxseg*6
         iii   = nlast
         do nbl=1,maxbl
            do nseg=1,maxseg
               do mm=1,7
                  do ll=1,2
                     iii = iii + 1
                     iwork(iii) =
     .                          ibcinfo(nbl,nseg,mm,ll)
                     iwork(iii+maxbl*maxseg*14) =
     .                          jbcinfo(nbl,nseg,mm,ll)
                     iwork(iii+maxbl*maxseg*28) =
     .                          kbcinfo(nbl,nseg,mm,ll)
                  end do
               end do
            end do
         end do
         nlast = nlast + maxbl*maxseg*42
         iii   = nlast
         do ics=1,maxcs
            do ll=1,9
               iii = iii + 1
               iwork(iii) = icsinfo(ics,ll)
            end do
         end do
         nlast = nlast + 9*maxcs
         iii   = nlast
         do ll=1,2
            do mm=1,mxbli
               iii = iii + 1
               iwork(iii) = nblk(ll,mm)
            end do
         end do
         nlast = nlast + 2*mxbli
         iii   = nlast
         do ll=1,2
            do nn = 1,2
               do mm=1,mxbli
                  iii = iii + 1
                  iwork(iii) = isva(ll,nn,mm)
               end do
            end do
         end do
         nlast = nlast + 4*mxbli
         iii   = nlast
         do ll=1,2
            do nn = 1,6
               do mm=1,mxbli
                  iii = iii + 1
                  iwork(iii) = limblk(ll,nn,mm)
               end do
            end do
         end do
         nlast = nlast + 12*mxbli
         do iii=1,mxbli
            iwork(nlast+iii) = nblon(iii)
         end do
         nlast = nlast + mxbli
         iii   = nlast
         do npl=1,nplots
            do ll=1,11
               iii = iii + 1
               iwork(iii)           = inpl3d(npl,ll)
               iwork(iii+11*nplots) = inpr(npl,ll)
            end do
         end do
         nlast = nlast + 22*nplots
         do j = 1,500
          do nbl=1,maxbl
            iwork(nlast+j+500*(nbl-1))           = iskip(nbl,j)
            iwork(nlast+j+500*maxbl+500*(nbl-1)) = jskip(nbl,j)
            iwork(nlast+j+1000*maxbl+500*(nbl-1))= kskip(nbl,j)
          enddo
         end do
         nlast = nlast + 1500*maxbl
c
      end if
c
      nvals =1550*maxbl + 119      + 9*maxbl*maxsegdg + 48*maxbl*maxseg
     .      + 19*mxbli + 22*nplots + 4*maxgr          + 9*maxcs
      if (myid.eq.myhost) then
         if (nvals.ne.nlast) then
            write(11,'('' error in subroutine trnsfr_vals, integer'',
     .      '' transfer: nvals.ne.nlast'')')
            write(11,'('' nvals,nlast '',i8,'' '',i8)') nvals,nlast
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)      
         end if
      end if
      call MPI_Bcast (iwork, nvals, MPI_INTEGER, myhost,
     .                mycomm, ierr)
c
      if (myid.ne.myhost) then
         do iii=1,maxbl
            itrans(iii) = iwork(iii)
            irotat(iii) = iwork(iii+maxbl)
            idefrm(iii) = iwork(iii+2*maxbl)
         end do
         nlast = 3*maxbl
         itransmc = iwork(nlast+1)
         irotatmc = iwork(nlast+2)
         ndefrm   = iwork(nlast+3)
         naesrf   = iwork(nlast+4)
         ihstry   = iwork(nlast+5)
         irest    = iwork(nlast+6)
         iflagts  = iwork(nlast+7)
         iunst    = iwork(nlast+8)
         nwrest   = iwork(nlast+9)
         ichk     = iwork(nlast+10)
         i2d      = iwork(nlast+11)
         ntstep   = iwork(nlast+12)
         ita      = iwork(nlast+13)
         mseq     = iwork(nlast+14)
         issc     = iwork(nlast+15)
         issr     = iwork(nlast+16)
         ncyctot  = iwork(nlast+17)
         lfgm     = iwork(nlast+18)
         nplot3d  = iwork(nlast+19)
         movie    = iwork(nlast+20)
         nprint   = iwork(nlast+21)
         nbli     = iwork(nlast+22)
         nblock   = iwork(nlast+23)
         ngrid    = iwork(nlast+24)
         ninter   = iwork(nlast+25)
         mgflag   = iwork(nlast+26)
         iconsf   = iwork(nlast+27)
         mtt      = iwork(nlast+28)
         ngam     = iwork(nlast+29)
         ideg(1)  = iwork(nlast+30)
         ideg(2)  = iwork(nlast+31)
         ideg(3)  = iwork(nlast+32)
         ifds(1)  = iwork(nlast+33)
         ifds(2)  = iwork(nlast+34)
         ifds(3)  = iwork(nlast+35)
         ilamlo   = iwork(nlast+36)
         ilamhi   = iwork(nlast+37)
         jlamlo   = iwork(nlast+38)
         jlamhi   = iwork(nlast+39)
         klamlo   = iwork(nlast+40)
         klamhi   = iwork(nlast+41)
         ivmx     = iwork(nlast+42)
         isminc   = iwork(nlast+43)
         ip3dgrd  = iwork(nlast+44)
         ialph    = iwork(nlast+45)
         iipv     = iwork(nlast+46)
         irest2   = iwork(nlast+47)
         nnodes   = iwork(nlast+48)
         ncs      = iwork(nlast+49)
         ialphit  = iwork(nlast+50)
         ip3dsurf = iwork(nlast+51)
         isktyp   = iwork(nlast+52)
         iprnsurf = iwork(nlast+53)
         nsprgit  = iwork(nlast+54)
         nlast = nlast + 54
         do iii=1,maxgr
            ncgg(iii)   = iwork(nlast+iii)
            iemg(iii)   = iwork(nlast+iii+maxgr)
            nblg(iii)   = iwork(nlast+iii+2*maxgr)
            inewgg(iii) = iwork(nlast+iii+3*maxgr)
         end do
          nlast = nlast + 4*maxgr
         do iii=1,20
            ncell(iii) = iwork(iii+nlast)
         end do
         nlast = nlast + 20
         do iii=1,maxbl
            iadvance(iii) = iwork(nlast+iii)
            iforce(iii)   = iwork(nlast+iii+maxbl)
         end do
         nlast = nlast + 2*maxbl
         iii   = nlast
         do nbl=1,maxbl
            do i = 1,3
               iii = iii + 1
               iviscg(nbl,i) = iwork(iii)
               idegg(nbl,i)  = iwork(iii+3*maxbl)
               iwfg(nbl,i)   = iwork(iii+6*maxbl)
            end do
         end do
         nlast = nlast + 9*maxbl
         do iii=1,maxbl
            idimg(iii)   = iwork(nlast+iii)
            jdimg(iii)   = iwork(nlast+iii+maxbl)
            kdimg(iii)   = iwork(nlast+iii+2*maxbl)
            ilamlog(iii) = iwork(nlast+iii+3*maxbl)
            ilamhig(iii) = iwork(nlast+iii+4*maxbl)
            jlamlog(iii) = iwork(nlast+iii+5*maxbl)
            jlamhig(iii) = iwork(nlast+iii+6*maxbl)
            klamlog(iii) = iwork(nlast+iii+7*maxbl)
            klamhig(iii) = iwork(nlast+iii+8*maxbl)
            jsg(iii)     = iwork(nlast+iii+9*maxbl)
            ksg(iii)     = iwork(nlast+iii+10*maxbl)
            isg(iii)     = iwork(nlast+iii+11*maxbl)
            jeg(iii)     = iwork(nlast+iii+12*maxbl)
            keg(iii)     = iwork(nlast+iii+13*maxbl)
            ieg(iii)     = iwork(nlast+iii+14*maxbl)
            nblcg(iii)   = iwork(nlast+iii+15*maxbl)
         end do
         nlast = nlast + 16*maxbl
         iii   = nlast
         do nbl=1,maxbl
            do i = 1,3
               iii = iii + 1
               idiagg(nbl,i) = iwork(iii)
               iflimg(nbl,i) = iwork(iii+3*maxbl)
               ifdsg(nbl,i)  = iwork(iii+6*maxbl)
            end do
         end do
         nlast = nlast + 9*maxbl
         do iii=1,maxbl
            nbci0(iii)   = iwork(nlast+iii)
            nbcidim(iii) = iwork(nlast+iii+maxbl)
            nbcj0(iii)   = iwork(nlast+iii+2*maxbl)
            nbcjdim(iii) = iwork(nlast+iii+3*maxbl)
            nbck0(iii)   = iwork(nlast+iii+4*maxbl)
            nbckdim(iii) = iwork(nlast+iii+5*maxbl)
            iovrlp(iii)  = iwork(nlast+iii+6*maxbl)
            igridg(iii)  = iwork(nlast+iii+7*maxbl)
            levelg(iii)  = iwork(nlast+iii+8*maxbl)
            mblk2nd(iii) = iwork(nlast+iii+9*maxbl)
         end do
         nlast = nlast + 10*maxbl
         do iii=1,5
            ncyc1(iii)  = iwork(nlast+iii)
            nitfo1(iii) = iwork(nlast+iii+5)
            levelt(iii) = iwork(nlast+iii+10)
            levelb(iii) = iwork(nlast+iii+15)
         end do
         nlast = nlast + 20
         iii   = nlast
         do m = 1,5
            do i = 1, 5
               iii = iii + 1
               mit(i,m) = iwork(iii)
            end do
         end do
         nlast = nlast + 25
         iii   = nlast
         do kk=1,maxsegdg
            do jj=1,maxbl
               iii             = iii + 1
               icsi(jj,kk)     = iwork(iii)
               icsf(jj,kk)     = iwork(iii+  maxbl*maxsegdg)
               jcsi(jj,kk)     = iwork(iii+2*maxbl*maxsegdg)
               jcsf(jj,kk)     = iwork(iii+3*maxbl*maxsegdg)
               kcsi(jj,kk)     = iwork(iii+4*maxbl*maxsegdg)
               kcsf(jj,kk)     = iwork(iii+5*maxbl*maxsegdg)
               idfrmseg(jj,kk) = iwork(iii+6*maxbl*maxsegdg)
               iaesurf(jj,kk)  = iwork(iii+7*maxbl*maxsegdg)
               icouple(jj,kk)  = iwork(iii+8*maxbl*maxsegdg)
            end do
         end do
         nlast = nlast + 9*maxbl*maxsegdg
         do iii = 1,maxbl
            nsegdfrm(iii) = iwork(nlast+iii)
         end do
         nlast = nlast + maxbl
         iii   = nlast
         do nbl=1,maxbl
            do nseg=1,maxseg
               do ll=1,2
                  iii = iii + 1
                  bcfilei(nbl,nseg,ll) = iwork(iii)
                  bcfilej(nbl,nseg,ll) = iwork(iii+maxbl*maxseg*2)
                  bcfilek(nbl,nseg,ll) = iwork(iii+maxbl*maxseg*4)
               end do
            end do
         end do
         nlast = nlast + maxbl*maxseg*6
         iii   = nlast
         do nbl=1,maxbl
            do nseg=1,maxseg
               do mm=1,7
                  do ll=1,2
                     iii = iii + 1
                     ibcinfo(nbl,nseg,mm,ll) =
     .                          iwork(iii)
                     jbcinfo(nbl,nseg,mm,ll) = 
     .                          iwork(iii+maxbl*maxseg*14)
                     kbcinfo(nbl,nseg,mm,ll) =
     .                          iwork(iii+maxbl*maxseg*28)
                  end do
               end do
            end do
         end do
         nlast = nlast + maxbl*maxseg*42
         iii   = nlast
         do ics=1,maxcs
            do ll=1,9
               iii = iii + 1
               icsinfo(ics,ll) = iwork(iii)
            end do
         end do
         nlast = nlast + 9*maxcs
         iii   = nlast
         do ll=1,2
            do mm=1,mxbli
               iii = iii + 1
               nblk(ll,mm) = iwork(iii)
            end do
         end do
         nlast = nlast + 2*mxbli
         iii   = nlast
         do ll=1,2
            do nn = 1,2
               do mm=1,mxbli
                  iii = iii + 1
                  isva(ll,nn,mm) = iwork(iii)
               end do
            end do
         end do
         nlast = nlast + 4*mxbli
         iii   = nlast
         do ll=1,2
            do nn = 1,6
               do mm=1,mxbli
                  iii = iii + 1
                  limblk(ll,nn,mm) = iwork(iii)
               end do
            end do
         end do
         nlast = nlast + 12*mxbli
         do iii=1,mxbli
            nblon(iii) = iwork(nlast+iii)
         end do
         nlast = nlast + mxbli
         iii   = nlast
         do npl=1,nplots
            do ll=1,11
               iii = iii + 1
               inpl3d(npl,ll) = iwork(iii)
               inpr(npl,ll)   = iwork(iii+11*nplots)
            end do
         end do
         nlast = nlast + 22*nplots
         do j = 1,500
          do nbl=1,maxbl
            iskip(nbl,j) =  iwork(nlast+j+500*(nbl-1))
            jskip(nbl,j) =  iwork(nlast+j+500*maxbl+500*(nbl-1))
            kskip(nbl,j) =  iwork(nlast+j+1000*maxbl+500*(nbl-1))
          end do
         end do
         nlast = nlast + 1500*maxbl
      end if
c
c****************************
c     transfer character data
c****************************
c
      call MPI_Bcast (bcfiles, mxbcfil*80, MPI_CHARACTER, myhost,
     .                mycomm, ierr)
      call MPI_Bcast (grdmov, 80, MPI_CHARACTER, myhost,
     .                mycomm, ierr)
#endif
c
      return
      end
