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
      subroutine mgbl(mwork,mworki,nplots,minnode,iitot,intmax,maxxe,
     .                mxbli,nsub1,lbcprd,lbcemb,lbcrad,maxbl,maxgr,
     .                maxseg,maxcs,ncycmax,intmx,mxxe,mptch,msub1,
     .                ibufdim,nbuf,mxbcfil,istat_size,
     .                work,iwork,lwdat,lw,lw2,nou,bou,bcfilei,
     .                bcfilej,bcfilek,istat2_bl,istat2_pa,istat2_em,
     .                istat2_pe,nblk,limblk,isva,nblon,resmx,imx,
     .                jmx,kmx,vormax,ivmax,jvmax,kvmax,lig,lbg,
     .                iovrlp,qb,ibpntsg,iipntsg,iibg,kkbg,jjbg,ibcg,
     .                dxintg,dyintg,dzintg,iiig,jjig,kkig,rkap0g,
     .                levelg,igridg,iflimg,ifdsg,iviscg,jdimg,
     .                kdimg,idimg,idiagg,nblcg,idegg,jsg,ksg,isg,
     .                jeg,keg,ieg,mit,ilamlog,ilamhig,jlamlog,
     .                jlamhig,klamlog,klamhig,iwfg,utrans,vtrans,
     .                wtrans,omegax,omegay,omegaz,xorig,yorig,
     .                zorig,dxmx,dymx,dzmx,dthxmx,dthymx,dthzmx,
     .                thetax,thetay,thetaz,rfreqt,rfreqr,xorig0,
     .                yorig0,zorig0,time2,thetaxl,thetayl,thetazl,
     .                itrans,irotat,idefrm,bcvali,bcvalj,bcvalk,nbci0,
     .                nbcidim,nbcj0,nbcjdim,nbck0,nbckdim,ibcinfo,
     .                jbcinfo,kbcinfo,ncgg,nblg,iemg,inewgg,
     .                inpl3d,inpr,iadvance,iforce,rms,clw,cdw,
     .                cdpw,cdvw,cxw,cyw,czw,cmxw,cmyw,cmzw,
     .                n_clcd,clcd,nblocks_clcd,blocks_clcd,chdw,
     .                swetw,fmdotw,cfttotw,cftmomw,cftpw,cftvw,
     .                swett,clt,cdt,cxt,
     .                cyt,czt,cmxt,cmyt,cmzt,cdpt,cdvt,sx,sy,sz,
     .                stot,pav,ptav,tav,ttav,xmav,fmdot,cfxp,cfyp,
     .                cfzp,cfdp,cflp,cftp,cfxv,cfyv,cfzv,cfdv,cflv,
     .                cftv,cfxmom,cfymom,cfzmom,cfdmom,cflmom,
     .                cftmom,cfxtot,cfytot,cfztot,cfdtot,cfltot,
     .                cfttot,icsinfo,windex,iindex,nblkpt,windx,
     .                iindx,llimit,iitmax,mmcxie,mmceta,ncheck,
     .                iifit,mblkpt,iic0,iiorph,iitoss,ifiner,
     .                dthetxx,dthetyy,dthetzz,dx,dy,dz,dthetx,
     .                dthety,dthetz,lout,xif1,xif2,etf1,etf2,
     .                jjmax1,kkmax1,iiint1,iiint2,jimage,kimage,
     .                jte,kte,jmm,kmm,nblk1,nblk2,xte,yte,zte,xmi,
     .                ymi,zmi,xmie,ymie,zmie,sxie,seta,sxie2,
     .                seta2,xie2s,eta2s,temp,x2,y2,z2,x1,y1,z1,
     .                factjlo,factjhi,factklo,factkhi,ifrom,
     .                geom_miss,period_miss,isav_blk,isav_prd,
     .                isav_pat,isav_pat_b,isav_dpat,isav_dpat_b,
     .                isav_emb,mblk2nd,ntr,utrnsae,vtrnsae,wtrnsae,
     .                omgxae,omgyae,omgzae,xorgae,yorgae,zorgae,thtxae,
     .                thtyae,thtzae,rfrqtae,rfrqrae,icsi,icsf,jcsi,jcsf,
     .                kcsi,kcsf,freq,gmass,damp,x0,gf0,nmds,maxaes,
     .                aesrfdat,perturb,memuse,bcfiles,islavept,nslave,
     .                iskip,jskip,kskip,bmat,stm,stmi,xs,xxn,gforcn,
     .                gforcnm,gforcs,nsegdfrm,idfrmseg,iaesurf,
     .                maxsegdg,nmaster,aehist,timekeep,xorgae0,yorgae0,
     .                zorgae0,icouple,nblelst,iskmax,jskmax,kskmax,ue,
     .                irdrea)
#   ifdef CMPLX
#   else
      use module_kwstm, only:kws_get_nummem
#   endif
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Main routine of cfl3d; drives mesh sequencing and
c     multigrid solutions to steady/unsteady Euler/Navier-Stokes 
c     (thin layer form) equations
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
      dimension istat(MPI_STATUS_SIZE)
#endif
c
      character*80 grid,plt3dg,plt3dq,output,residual,turbres,blomx,
     .             output2,printout,pplunge,ovrlap,patch,restrt,
     .             subres,subtur,grdmov,alphahist,errfile,preout,
     .             aeinp,aeout,sdhist,avgg,avgq
      character*50 string
      character*80 bcfiles(mxbcfil)
      character*120 bou(ibufdim,nbuf)
c
      real*4 tim(3,3),tm(3),walltime,totaltime
c
      integer lout,xif1,xif2,etf1,etf2
      integer bcfilei,bcfilej,bcfilek
      integer irdrea(maxgr),stats
c
      dimension bcfilei(maxbl,maxseg,2),bcfilej(maxbl,maxseg,2),
     .          bcfilek(maxbl,maxseg,2)
      dimension istat2_bl(istat_size,mxbli*5),
     .          istat2_pa(istat_size,intmax*nsub1*3),
     .          istat2_em(istat_size,lbcemb*3),
     .          istat2_pe(istat_size,lbcprd*5)
      dimension nou(nbuf)
      dimension work(mwork),iwork(mworki),lwdat(maxbl,maxseg,6)
      dimension lw(65,maxbl),lw2(43,maxbl)
      dimension nblk(2,mxbli),limblk(2,6,mxbli),
     .          isva(2,2,mxbli),nblon(mxbli)
      dimension resmx(maxbl),imx(maxbl),jmx(maxbl),kmx(maxbl)
      dimension vormax(maxbl),ivmax(maxbl),jvmax(maxbl),kvmax(maxbl)
      dimension lig(maxbl),lbg(maxbl),iovrlp(maxbl),qb(iitot,5,3)
      dimension ibpntsg(maxbl,4),iipntsg(maxbl)
      dimension iibg(iitot),kkbg(iitot),jjbg(iitot),ibcg(iitot)
      dimension dxintg(iitot),dyintg(iitot),dzintg(iitot),
     .          iiig(iitot),jjig(iitot),kkig(iitot)
      dimension rkap0g(maxbl,3),levelg(maxbl),igridg(maxbl),
     .          iflimg(maxbl,3),ifdsg(maxbl,3),iviscg(maxbl,3),
     .          jdimg(maxbl),kdimg(maxbl),idimg(maxbl),idiagg(maxbl,3),
     .          nblcg(maxbl),idegg(maxbl,3),
     .          jsg(maxbl),ksg(maxbl),isg(maxbl),jeg(maxbl),keg(maxbl),
     .          ieg(maxbl),mit(5,maxbl),ilamlog(maxbl),ilamhig(maxbl),
     .          jlamlog(maxbl),jlamhig(maxbl),klamlog(maxbl),
     .          klamhig(maxbl),iwfg(maxbl,3)
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
     .          perturb(nmds,maxaes,4),aehist(ncycmax,3,nmds,maxaes)
      dimension aesrfdat(5,maxaes),islavept(nslave,nmaster,5)
      dimension ue(3*nslave)
      dimension iskip(maxbl,500),jskip(maxbl,500),kskip(maxbl,500)
      dimension bmat(2*nmds,2*nmds,maxaes),gforcn(2*nmds,maxaes),
     .          gforcnm(2*nmds,maxaes),gforcs(2*nmds,maxaes),
     .          stm(2*nmds,2*nmds,maxaes),stmi(2*nmds,2*nmds,maxaes),
     .          xs(2*nmds,maxaes),xxn(2*nmds,maxaes)
      dimension bcvali(maxbl,maxseg,12,2),
     .          bcvalj(maxbl,maxseg,12,2),bcvalk(maxbl,maxseg,12,2),
     .          nbci0(maxbl),nbcidim(maxbl),nbcj0(maxbl),nbcjdim(maxbl),
     .          nbck0(maxbl),nbckdim(maxbl),ibcinfo(maxbl,maxseg,7,2),
     .          jbcinfo(maxbl,maxseg,7,2),kbcinfo(maxbl,maxseg,7,2)
      dimension ncgg(maxgr),nblg(maxgr),iemg(maxgr),inewgg(maxgr)
      dimension inpl3d(nplots,11),inpr(nplots,11)
      dimension iadvance(maxbl),iforce(maxbl)
      integer blocks_clcd
      dimension rms(ncycmax),clw(ncycmax),
     .          cdw(ncycmax),cdpw(ncycmax),cdvw(ncycmax),
     .          cxw(ncycmax),cyw(ncycmax),czw(ncycmax),
     .          cmxw(ncycmax),cmyw(ncycmax),cmzw(ncycmax),
     .          clcd(2,n_clcd,ncycmax), blocks_clcd(2,nblocks_clcd),
     .          chdw(ncycmax),swetw(ncycmax),
     .          fmdotw(ncycmax),cfttotw(ncycmax),
     .          cftmomw(ncycmax),cftpw(ncycmax),cftvw(ncycmax),
     .          timekeep(ncycmax)
      dimension swett(maxbl),clt(maxbl),cdt(maxbl),cxt(maxbl),
     .          cyt(maxbl),czt(maxbl),cmxt(maxbl),cmyt(maxbl),
     .          cmzt(maxbl),cdpt(maxbl),cdvt(maxbl)
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
      dimension dx(intmx,msub1),dy(intmx,msub1),dz(intmx,msub1),
     .          dthetx(intmx,msub1),dthety(intmx,msub1),
     .          dthetz(intmx,msub1)
      dimension lout(msub1),xif1(msub1),xif2(msub1),etf1(msub1),
     .          etf2(msub1),jjmax1(nsub1),kkmax1(nsub1),iiint1(nsub1),
     .          iiint2(nsub1),jimage(msub1,mptch+2,mptch+2),
     .          kimage(msub1,mptch+2,mptch+2),jte(msub1),kte(msub1),
     .          jmm(mptch+2),kmm(mptch+2),nblk1(mptch+2),nblk2(mptch+2),
     .          xte(mptch+2,mptch+2,msub1),yte(mptch+2,mptch+2,msub1),
     .          zte(mptch+2,mptch+2,msub1),xmi(mptch+2,mptch+2,msub1),
     .          ymi(mptch+2,mptch+2,msub1),zmi(mptch+2,mptch+2,msub1),
     .          xmie(mptch+2,mptch+2,msub1),ymie(mptch+2,mptch+2,msub1),
     .          zmie(mptch+2,mptch+2,msub1),sxie(mptch+2,mptch+2,msub1),
     .          seta(mptch+2,mptch+2,msub1),sxie2(mptch+2,mptch+2),
     .          seta2(mptch+2,mptch+2),xie2s(mptch+2,mptch+2),
     .          eta2s(mptch+2,mptch+2),temp((mptch+2)*(mptch+2)),
     .          x2(mptch+2,mptch+2),y2(mptch+2,mptch+2),
     .          z2(mptch+2,mptch+2),x1(mptch+2,mptch+2),
     .          y1(mptch+2,mptch+2),z1(mptch+2,mptch+2),
     .          factjlo(intmx,msub1),factjhi(intmx,msub1),
     .          factklo(intmx,msub1),factkhi(intmx,msub1),ifrom(msub1)
      dimension geom_miss(2*mxbli),period_miss(lbcprd)
      dimension isav_blk(2*mxbli,17)
      dimension isav_prd(lbcprd,12)
      dimension isav_pat(intmax,17),isav_pat_b(intmax,nsub1,6)
      dimension isav_dpat(intmx,17),isav_dpat_b(intmx,msub1,6)
      dimension isav_emb(lbcemb,12)
      dimension mblk2nd(maxbl)
      dimension nblelst(maxbl,2)
      dimension iskmax(maxbl)
      dimension jskmax(maxbl)
      dimension kskmax(maxbl)
c
      common /elastic/ ndefrm,naesrf
      common /alphait/ ialphit,cltarg,rlxalph,dalim,dalpha,icycupdt
      common /cpurate/ rate(5),ratesub(5),ncell(20)
      common /cfl/ dt0,dtold
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /maxiv/ ivmx
      common /sklton/ isklton
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /moov/movie,nframes,icall1,lhdr,icoarsemovie,i2dmovie
      common /filenam/ grid,plt3dg,plt3dq,output,residual,turbres,blomx,
     .                 output2,printout,pplunge,ovrlap,patch,restrt,
     .                 subres,subtur,grdmov,alphahist,errfile,preout,
     .                 aeinp,aeout,sdhist,avgg,avgq
      common /unit5/ iunit5
      common /conversion/ radtodeg
      common /mydist2/ nnodes,myhost,myid,mycomm
      common /proces/ numprocs
      common /time1/ tim,tm
      common /motionmc/ xmc0,ymc0,zmc0,utransmc,vtransmc,wtransmc,
     .                  omegaxmc,omegaymc,omegazmc,xorigmc,yorigmc,
     .                  zorigmc,xorig0mc,yorig0mc,zorig0mc,thetaxmc,
     .                  thetaymc,thetazmc,dxmxmc,dymxmc,dzmxmc,
     .                  dthxmxmc,dthymxmc,dthzmxmc,rfreqtmc,
     .                  rfreqrmc,itransmc,irotatmc,time2mc
      common /bin/ ibin,iblnk,iblnkfr,ip3dgrad
      common /deformz/ beta1,beta2,alpha1,alpha2,isktyp,negvol,meshdef,
     .                 nsprgit,ndgrd,ndwrt 

      allocatable :: rmstr(:,:)
      allocatable :: nneg(:,:)
c
c***********************************************************************
c     initialization
c***********************************************************************
c
c     initialize timing utility
c
      string = ' '
      call cputim(0,1,string,myhost,myid,mycomm,11)
c
c     initialize data for grids stored in array w
c
      do ll=1,mwork
         work(ll) = 0.0
      end do
c
      do ll=1,mworki
         iwork(ll) = 0
      end do
c
c     initialize internal write buffer
c
      do lou=1,nbuf
         do kou=1,ibufdim
            bou(kou,lou) = ' '
         end do
         nou(lou) = 0
      end do
c
c     initialize auxiliary file names for 2000 series bc's
c
      do nf=1,mxbcfil
         bcfiles(nf) = 'null'
      end do
c
c     initialize radtodeg
c
      pi       = 4.*atan(1.0)
      radtodeg = 180.e0/pi
c
      if (myid.eq.myhost) then
c
c***********************************************************************
c        write version information, etc. to main output file
c***********************************************************************
c
c        output banner
c
         write(11,83)
         write(11,83)
         write(11,87)
         write(11,9900)
 9900    format(2(2h *),43h   CFL3D :  CHARACTERISTIC-BASED SCHEME FOR,
     .   16h STEADY/UNSTEADY,3x,2(2h *)
     .   /2(2h *),12x,46hSOLUTIONS TO THE EULER/NAVIER-STOKES EQUATIONS,
     .   4x,2(2h *))
         write(11,87)
         write(11,9990)
 9990    format(2(2h *),43h   VERSION 6.7 :  Computational Fluids Lab,,
     .   15h Mail Stop 128,,4x,2(2h *),
     .   /2(2h *),18x,41hNASA Langley Research Center, Hampton, VA,
     .   3x,2(2h *),/2(2h *),18x,33hRelease Date:  February  1, 2017.,
     .   11x,2(2h *))
#ifdef CMPLX
         write(11,87)
         write(11,9991)
 9991    format(2(2h *),19x,24hCOMPLEX VARIABLE VERSION,19x,2(2h *))
#endif
         write(11,87)
         write(11,83)
         write(11,83)
   83    format(35(2h *))
   87    format(2(2h *),62x,2(2h *))
         write(11,'(''CFL3D is a structured-grid, cell-centered,'',
     +    '' upwind-biased, Reynolds-averaged'')')
         write(11,'(''Navier-Stokes (RANS) code. It can be run'',
     +    '' in parallel on multiple grid zones'')')
         write(11,'(''with point-matched, patched, overset, or'',
     +    '' embedded connectivities. Both'')')
         write(11,'(''multigrid and mesh sequencing are'',
     +    '' available in time-accurate or'')')
         write(11,'(''steady-state modes.'')')
         write(11,'('' '')')
         write(11,'(''-------------------------------------------'',
     +    ''----------------------------'')')
         write(11,'(''Copyright 2001 United States Government'',
     +    '' as represented by the Administrator'')')
         write(11,'(''of the National Aeronautics and Space'',
     +    '' Administration. All Rights Reserved.'')')
         write(11,'('' '')')
         write(11,'(''The CFL3D platform is licensed under the'',
     +    '' Apache License, Version 2.0'')')
         write(11,'(''(the "License"); you may not use this file'',
     +    '' except in compliance with the'')')
         write(11,'(''License. You may obtain a copy of the'',
     +    '' License at'')')
         write(11,'(''http://www.apache.org/licenses/LICENSE-2.0.'')')
         write(11,'('' '')')
         write(11,'(''Unless required by applicable law or agreed'',
     +    '' to in writing, software'')')
         write(11,'(''distributed under the License is distributed'',
     +    '' on an "AS IS" BASIS, WITHOUT'')')
         write(11,'(''WARRANTIES OR CONDITIONS OF ANY KIND, either'',
     +    '' express or implied. See the'')')
         write(11,'(''License for the specific language governing'',
     +    '' permissions and limitations'')')
         write(11,'(''under the License.'')')
         write(11,'(''-------------------------------------------'',
     +    ''----------------------------'')')
c
#ifdef CRAY
c        cray_time implies cray (always double precision)
         write(11,12) real(float(memuse))/1.e6
#else
#   ifdef DBLE_PRECSN
         write(11,12) real(float(memuse))/1.e6
#   else
         write(11,13) real(float(memuse))/1.e6
#   endif
#endif
   12    format(/,' memory allocation: ',f12.6,' Mbytes',
     .   ' double precision')
   13    format(/,' memory allocation: ',f12.6,' Mbytes',
     .   ' single precision')
c
         write(11,88)
   88    format(/19hinput/output files:)
c
         write(11,'(''  '',a60)')grid
         write(11,'(''  '',a60)')plt3dg
         write(11,'(''  '',a60)')plt3dq
         write(11,'(''  '',a60)')output
         write(11,'(''  '',a60)')residual
         write(11,'(''  '',a60)')turbres
         write(11,'(''  '',a60)')blomx
         write(11,'(''  '',a60)')output2
         write(11,'(''  '',a60)')printout
         write(11,'(''  '',a60)')pplunge
         write(11,'(''  '',a60)')ovrlap
         write(11,'(''  '',a60)')patch
         write(11,'(''  '',a60)')restrt
c
c***********************************************************************
c        read in global information
c***********************************************************************
c
         iwk1    = 1
         iwk2    = iwk1 + maxbl
         iwk3    = iwk2 + maxbl
         iwk4    = iwk3 + maxbl
         iwk5    = iwk4 + 11*nplots
         mworki1 = mworki - iwk5 + 1
         if (mworki1.le.0) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''stopping...not enough integer '',
     .           ''work space for subroutine global'')')
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
c
c        read input file up to title line
         nread = 14
         do n=1,nread
            read(iunit5,*)
         end do
c
         iunit   = 11
         icallgl = 1
         call global(myid,maxbl,maxgr,maxseg,maxcs,nplots,mxbli,
     .               bcvali,bcvalj,bcvalk,nbci0,nbcj0,nbck0,
     .               nbcidim,nbcjdim,nbckdim,ibcinfo,jbcinfo,
     .               kbcinfo,bcfilei,bcfilej,bcfilek,nblk,nbli,
     .               limblk,isva,nblon,rkap0g,nblock,levelg,
     .               igridg,iflimg,ifdsg,iviscg,jdimg,kdimg,
     .               idimg,idiagg,nblcg,idegg,jsg,ksg,isg,jeg,
     .               keg,ieg,mit,ilamlog,ilamhig,jlamlog,
     .               jlamhig,klamlog,klamhig,iwfg,utrans,vtrans,
     .               wtrans,omegax,omegay,omegaz,xorig,yorig,
     .               zorig,dxmx,dymx,dzmx,dthxmx,dthymx,
     .               dthzmx,thetax,thetay,thetaz,rfreqt,rfreqr,
     .               xorig0,yorig0,zorig0,itrans,irotat,idefrm,
     .               ngrid,ncgg,nblg,iemg,inewgg,iovrlp,ninter,
     .               nplot3d,inpl3d,ip3dsurf,nprint,inpr,
     .               iadvance,iforce,lfgm,ncs,icsinfo,ihstry,
     .               ncycmax,iwork(iwk1),time2,thetaxl,thetayl,thetazl,
     .               intmax,nsub1,iindex,lig,lbg,ibpntsg,
     .               iipntsg,icallgl,iunit,nou,bou,ibufdim,nbuf,
     .               iwork(iwk2),iwork(iwk3),iwork(iwk4),ntr,bcfiles,
     .               mxbcfil,utrnsae,vtrnsae,wtrnsae,omgxae,
     .               omgyae,omgzae,xorgae,yorgae,zorgae,thtxae,
     .               thtyae,thtzae,rfrqtae,rfrqrae,icsi,icsf,jcsi,jcsf,
     .               kcsi,kcsf,freq,gmass,damp,x0,gf0,nmds,maxaes,
     .               aesrfdat,perturb,iskip,jskip,kskip,nsegdfrm,
     .               idfrmseg,iaesurf,maxsegdg,xorgae0,yorgae0,zorgae0,
     .               icouple,iprnsurf)
         do iii = 1,iwk5
            iwork(iii) = 0
         end do
c
#if defined DIST_MPI
c
c        check that excess processors are not requested
c
         if (nnodes .gt. ngrid) then
            write(11,68) nnodes+1,ngrid+1,ngrid
   68       format(/,'stopping...requested ',i3,' processors',/,
     .               'you do not need more than ',i3,
     .               ' processors with ',i3,' grids!')
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
#endif
c
c***********************************************************************
c        set up block-to-node mapping array
c***********************************************************************
c
#if defined DIST_MPI
         ierrflg = -1
         call compg2n(nblock,ngrid,ncgg,nblg,idimg,jdimg,kdimg,
     .                nblcg,nnodes,iwork,myid,myhost,mblk2nd,
     .                mycomm,maxgr,maxbl,ierrflg,ibufdim,nbuf,bou,nou)
         do iii = 1,maxbl*3
            iwork(iii) = 0
         end do
#else
         do ibloc = 1,maxbl
            mblk2nd(ibloc)  = 0
         enddo
#endif
c
      end if
c
#if defined DIST_MPI
c***********************************************************************
c     transfer data from host to nodes
c***********************************************************************
c
      nvalsr = 30*maxbl + maxbl*maxseg*72 + 55 + 11*nmds*maxaes
     .       + 5*maxaes + 17*maxbl*maxsegdg + nkey
      nvalsi =1550*maxbl+ 119      + 9*maxbl*maxsegdg + 48*maxbl*maxseg
     .       + 19*mxbli + 22*nplots + 4*maxgr          + 9*maxcs
      if (mwork .lt. nvalsr) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),'(''stopping...not enough work '',
     .                  ''space for subroutine trnsfr_vals'')')
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
      if (mworki .lt. nvalsi) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),'(''stopping...not enough integer work '',
     .                  ''space for subroutine trnsfr_vals'')')
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
      call trnsfr_vals(work,mwork,iwork,mworki,maxgr,maxbl,
     .                 maxseg,nbli,nblock,levelg,igridg,
     .                 iovrlp,nbci0,nbcj0,nbck0,
     .                 nbcidim,nbcjdim,nbckdim,iflimg,ifdsg,
     .                 iviscg,jdimg,kdimg,idimg,nblcg,idegg,
     .                 isg,jsg,ksg,ieg,jeg,keg,jlamlog,klamlog,
     .                 ilamlog,ilamhig,jlamhig,klamhig,mit,
     .                 iwfg,iadvance,iforce,lfgm,ihstry,ninter,
     .                 ngrid,ncgg,nblg,iemg,inewgg,idiagg,
     .                 nplot3d,nprint,ncs,ip3dsurf,mblk2nd,
     .                 utrans,vtrans,wtrans,omegax,omegay,
     .                 omegaz,xorig,yorig,zorig,dxmx,dymx,
     .                 dzmx,dthxmx,dthymx,dthzmx,thetax,thetay,
     .                 thetaz,rfreqt,rfreqr,xorig0,yorig0,zorig0,
     .                 time2,thetaxl,thetayl,thetazl,itrans,
     .                 irotat,idefrm,utrnsae,vtrnsae,wtrnsae,
     .                 omgxae,omgyae,omgzae,xorgae,yorgae,zorgae,
     .                 thtxae,thtyae,thtzae,rfrqtae,rfrqrae,
     .                 icsi,icsf,jcsi,jcsf,kcsi,kcsf,freq,gmass,
     .                 x0,gf0,damp,perturb,aesrfdat,bcvali,bcvalj,
     .                 bcvalk,ibcinfo,jbcinfo,kbcinfo,rkap0g,bcfilei,
     .                 bcfilej,bcfilek,mxbli,nblk,limblk,nblon,
     .                 isva,bcfiles,icsinfo,mxbcfil,maxcs,nmds,
     .                 maxaes,inpl3d,inpr,nplots,nou,bou,nbuf,
     .                 ibufdim,iskip,jskip,kskip,nsegdfrm,idfrmseg,
     .                 iaesurf,maxsegdg,xorgae0,yorgae0,zorgae0,
     .                 icouple,iprnsurf)
c
      do iii = 1,nvalsr
         work(iii) = 0.0
      end do
      do iii = 1,nvalsi
         iwork(iii) = 0
      end do
#endif
c
c***********************************************************************
c     compute starting locations for arrays and print out block to
c     node mapping for mpi case
c***********************************************************************
c
      if (mworki .lt. 2*maxbl) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),'(''stopping...not enough integer work '',
     .                  ''space for subroutine pointers'')')
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
      icallpt = 1
#if defined DIST_MPI
      mpihost = 1
#else
      mpihost = 0
#endif
      call pointers (lw,lw2,maxl,lembed,nstart,nwork,mwork,maxbl,
     .               maxgr,maxseg,lwdat,levelg,igridg,iviscg,
     .               idimg,jdimg,kdimg,nblcg,itrans,irotat,idefrm,
     .               nbci0,nbcj0,nbck0,nbcidim,nbcjdim,nbckdim,
     .               ibcinfo,jbcinfo,kbcinfo,ngrid,ncgg,nblg,
     .               iemg,nblock,myhost,myid,mblk2nd,nou,bou,nbuf,
     .               ibufdim,iwork,ilamlog,jlamlog,
     .               klamlog,ilamhig,jlamhig,klamhig,idegg,iwfg,
     .               idiagg,iflimg,ifdsg,rkap0g,jsg,ksg,isg,jeg,
     .               keg,ieg,iwork(maxbl+1),icallpt,nmds,maxaes,mpihost)
c
      do iii = 1,mworki
         iwork(iii) = 0
      end do
c
      do m=1,nnodes
         do nbl=1,nblock
            if (mblk2nd(nbl).eq.m) then
#   ifdef FASTIO
               call writ_buffast(nbl,11,nou,bou,nbuf,ibufdim,myhost,
     .                       myid,mycomm,mblk2nd,maxbl,01)
#   else
               call writ_buf(nbl,11,nou,bou,nbuf,ibufdim,myhost,myid,
     .                       mycomm,mblk2nd,maxbl)
#   endif
               go to 777
            end if
         end do
  777    continue
      end do
c
#if defined DIST_MPI
      if (myid.eq.myhost) then
         write(11,200)
         write(11,201) nblock,nnodes
         do i=1,nblock
            write(11,202) i,mblk2nd(i)
         end do
  200    format(/,1x,'BLOCK TO NODE MAPPING')
  201    format(1x,'no. of blocks = ',i4,/,1x,'no. of  nodes = ',i4,/,
     .   1x,'block    node')
  202    format(i5,3x,i5)
      end if
#endif
c
c***********************************************************************
c      Read grid data and perform preliminary calculations.
c***********************************************************************
c
      isklton = 1
c
      iwk1    = 1
      iwk2    = iwk1 + maxgr
      iwk3    = iwk2 + maxgr
      iwk4    = iwk3 + maxgr
      iwk5    = iwk4 + maxbl*8
      mworki1 = mworki - iwk5 + 1
      if (mworki1.le.0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),'(''stopping...not enough integer '',
     .        ''work space for subroutine setup'')')
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
c
c     constants derived from gamma
c
      gm1  = gamma-1.0e0
      gp1  = gamma+1.0e0
      gm1g = gm1/gamma
      gp1g = gp1/gamma
      ggm1 = gamma*gm1
c     set nummem (usually 2 for 1 or 2 turbulence equations)
c     need to also set memory appropriately in pointers
      nummem=2
      if (ivmx .eq. 30) nummem=3
      if (ivmx .eq. 40) nummem=4
#   ifdef CMPLX
c     RSM not working in complex mode
#   else
      if (ivmx .eq. 72) nummem=kws_get_nummem()
#   endif
      allocate( rmstr(ncycmax,nummem),stat=stats )
      call umalloc(ncycmax*nummem,0,'rmstr',memuse,stats)
      allocate( nneg(ncycmax,nummem),stat=stats )
      call umalloc(ncycmax*nummem,1,'nneg',memuse,stats)
c
      call setup(lw,lw2,work,nstart,work(nstart+1),nwork,iwork(iwk5),
     .           mworki1,iwork(iwk1),iwork(iwk2),iwork(iwk3),
     .           maxbl,mxbli,maxgr,maxseg,nsub1,maxxe,intmax,
     .           iitot,ncycmax,lwdat,lig,lbg,iovrlp,
     .           qb,nblock,iviscg,jdimg,kdimg,idimg,utrans,
     .           vtrans,wtrans,omegax,omegay,omegaz,xorig,
     .           yorig,zorig,dxmx,dymx,dzmx,dthxmx,dthymx,
     .           dthzmx,thetax,thetay,thetaz,rfreqt,rfreqr,
     .           xorig0,yorig0,zorig0,time2,thetaxl,thetayl,
     .           thetazl,itrans,irotat,idefrm,bcvali,bcvalj,
     .           bcvalk,nbci0,nbcidim,nbcj0,nbcjdim,
     .           nbck0,nbckdim,ibcinfo,jbcinfo,kbcinfo,bcfilei,
     .           bcfilej,bcfilek,ngrid,ncgg,nblg,iemg,inewgg,
     .           rms,clw,cdw,cdpw,cdvw,cxw,cyw,czw,cmxw,cmyw,
     .           cmzw,n_clcd,clcd,nblocks_clcd,blocks_clcd,
     .           chdw,swetw,fmdotw,cfttotw,cftmomw,cftpw,
     .           cftvw,rmstr,nneg,ntr,windex,
     .           ninter,iindex,nblkpt,dthetxx,dthetyy,dthetzz,
     .           iibg,kkbg,jjbg,ibcg,dxintg,dyintg,dzintg,iiig,
     .           jjig,kkig,ibpntsg,iipntsg,mblk2nd,nou,bou,nbuf,
     .           ibufdim,iwork(iwk4),igridg,bcfiles,mxbcfil,
     .           utrnsae,vtrnsae,wtrnsae,omgxae,omgyae,omgzae,
     .           xorgae,yorgae,zorgae,thtxae,thtyae,thtzae,
     .           rfrqtae,rfrqrae,icsi,icsf,jcsi,jcsf,kcsi,kcsf,
     .           freq,gmass,damp,x0,gf0,nmds,maxaes,aesrfdat,perturb,
     .           islavept,nslave,iskip,jskip,kskip,bmat,stm,stmi,
     .           gforcn,gforcnm,xxn,nsegdfrm,idfrmseg,iaesurf,
     .           maxsegdg,nmaster,aehist,timekeep,inpl3d,nplots,
     .           nplot3d,levelg,iadvance,xs,gforcs,xorgae0,yorgae0,
     .           zorgae0,icouple,lfgm,nblk,limblk,isva,nblelst,iskmax,
     .           jskmax,kskmax,ue,irdrea,nbli,nummem)
c
      do iii = 1,mworki
         iwork(iii) = 0
      end do
      do iii = nstart+1,mwork
         work(iii) = 0.0
      end do
c
c     write header to cfl3d.alpha file if specified Cl option is used
c
      if (ialphit.ne.0 .and. myid.eq.myhost) then
            write(27,'(''       it    log(res)           cl'',
     .      ''            alpha           dalpha'')')
      end if
c
c     open plot3d/printout files
c
      if (myid.eq.myhost) then
c
      if (ibin.eq.1) then
      open(unit=3,file=plt3dg,form='unformatted',status='unknown')
      open(unit=4,file=plt3dq,form='unformatted',status='unknown')
      else
         open(unit=3,file=plt3dg,form='formatted',status='unknown')
         open(unit=4,file=plt3dq,form='formatted',status='unknown')
      end if
      open(unit=17,file=printout,form='formatted',status='unknown')
c
      end if
c
c***********************************************************************
c     mesh sequencing - cycle through mesh sequences, coarser to finer
c***********************************************************************
c
      dt0 = dt
      nframes = 0
      icall1 = 0
c
      do 8000 iseq=1,mseq
c
c     skip mesh sequence levels with ncyc = 0
c
      if (ncyc1(iseq).gt.0) then
c
c        interpolate from coarser solution
c
         if (iseq.gt.1 .or. iseq.eq.mseq) then
            isklton = 1
            call qinter(iseq,lembed,maxl,lw,lw2,work,nstart,
     .                  work(nstart+1),nwork,maxbl,maxgr,
     .                  levelg,igridg,idimg,jdimg,kdimg,ngrid,
     .                  nblg,iemg,inewgg,itrans,irotat,idefrm,xorig,
     .                  yorig,zorig,xorig0,yorig0,zorig0,thetax,thetay,
     .                  thetaz,time2,nou,bou,nbuf,ibufdim,myid,myhost,
     .                  mycomm,mblk2nd,nsegdfrm,idfrmseg,xorgae,
     .                  yorgae,zorgae,thtxae,thtyae,thtzae,maxsegdg,
     .                  nummem)
c
            do iii = nstart+1,mwork
               work(iii) = 0.0
            end do
c
         end if
c
c***********************************************************************
c        multigrid/time-step cycles
c***********************************************************************
c
         string = ' '
         call cputim(+1,nnodes,string,myhost,myid,mycomm,11)
c
         nwork1 = nwork - maxbl
         if (nwork1.le.0) then 
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''stopping...not enough '',
     .           ''work space for subroutine mgblk'')')
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if 
         iwk1 = 1
         iwk2 = iwk1 + maxbl*7*3
         iwk3 = iwk2 + maxbl*8
         mworki1 = mworki - iwk3 + 1
         if (mworki1.le.0) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''stopping...not enough integer '',
     .           ''work space for subroutine mgblk'')')
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         call mgblk(iseq,lw,lw2,work,nstart,work(nstart+1),
     .              work(nstart+1+maxbl),nwork1,
     .              iwork(iwk1),iwork(iwk2),iwork(iwk3),mworki1,
     .              maxbl,maxgr,maxseg,iitot,intmax,nsub1,
     .              maxxe,ncycmax,iovrlp,lig,lbg,ibpntsg,iipntsg,
     .              iibg,kkbg,jjbg,ibcg,dxintg,dyintg,dzintg,iiig,
     .              jjig,kkig,qb,lwdat,nblk,nbli,limblk,isva,nblon,
     .              nblock,levelg,igridg,iviscg,idimg,jdimg,kdimg,
     .              jsg,ksg,isg,jeg,keg,ieg,nblcg,mit,bcvali,bcvalj,
     .              bcvalk,nbci0,nbcidim,nbcj0,nbcjdim,nbck0,
     .              nbckdim,ibcinfo,jbcinfo,kbcinfo,bcfilei,bcfilej,
     .              bcfilek,utrans,vtrans,wtrans,omegax,omegay,
     .              omegaz,xorig,yorig,zorig,dxmx,dymx,dzmx,dthxmx,
     .              dthymx,dthzmx,thetax,thetay,thetaz,rfreqt,
     .              rfreqr,xorig0,yorig0,zorig0,time2,thetaxl,
     .              thetayl,thetazl,itrans,irotat,idefrm,
     .              rms,clw,cdw,cdpw,cdvw,cxw,cyw,czw,cmxw,
     .              cmyw,cmzw,n_clcd,clcd,nblocks_clcd,blocks_clcd,
     .              chdw,swetw,fmdotw,cfttotw,cftmomw,
     .              cftpw,cftvw,rmstr,nneg,ntr,
     .              ihstry,iadvance,iforce,lfgm,resmx,imx,jmx,kmx,
     .              vormax,ivmax,jvmax,kvmax,sx,sy,sz,stot,pav,
     .              ptav,tav,ttav,xmav,fmdot,cfxp,cfyp,cfzp,cfdp,
     .              cflp,cftp,cfxv,cfyv,cfzv,cfdv,cflv,cftv,cfxmom,
     .              cfymom,cfzmom,cfdmom,cflmom,cftmom,cfxtot,
     .              cfytot,cfztot,cfdtot,cfltot,cfttot,icsinfo,ncs,
     .              windex,ninter,iindex,nblkpt,windx,nintr,iindx,
     .              llimit,iitmax,mmcxie,mmceta,ncheck,iifit,mblkpt,
     .              iic0,iiorph,iitoss,ifiner,msub1,dthetxx,dthetyy,
     .              dthetzz,swett,clt,cdt,cxt,cyt,czt,cmxt,cmyt,
     .              cmzt,cdpt,cdvt,mblk2nd,geom_miss,epsc0,
     .              period_miss,epsrot,isav_blk,isav_prd,lbcprd,
     .              isav_pat,isav_pat_b,isav_emb,lbcemb,mxbli,
     .              maxcs,intmx,mxxe,mptch,ncgg,nblg,iemg,ngrid,
     .              dx,dy,dz,dthetx,dthety,dthetz,isav_dpat,
     .              isav_dpat_b,lout,ifrom,xif1,xif2,etf1,
     .              etf2,jjmax1,kkmax1,iiint1,iiint2,nblk1,
     .              nblk2,jimage,kimage,jte,kte,jmm,kmm,
     .              xte,yte,zte,xmi,ymi,zmi,xmie,ymie,
     .              zmie,sxie,seta,sxie2,seta2,xie2s,
     .              eta2s,temp,x2,y2,z2,x1,y1,z1,ip3dsurf,
     .              factjlo,factjhi,factklo,factkhi,nplots,
     .              nou,bou,nbuf,ibufdim,istat2_bl,
     .              istat2_pa,istat2_pe,istat2_em,istat_size,
     .              nplot3d,nprint,inpl3d,inpr,bcfiles,mxbcfil,
     .              utrnsae,vtrnsae,wtrnsae,omgxae,omgyae,omgzae,
     .              xorgae,yorgae,zorgae,thtxae,thtyae,thtzae,
     .              rfrqtae,rfrqrae,icsi,icsf,jcsi,jcsf,kcsi,kcsf,
     .              freq,gmass,damp,x0,gf0,nmds,maxaes,aesrfdat,perturb,
     .              islavept,nslave,iskip,jskip,kskip,bmat,stm,stmi,
     .              xs,xxn,gforcn,gforcnm,gforcs,nsegdfrm,idfrmseg,
     .              iaesurf,maxsegdg,nmaster,aehist,timekeep,xorgae0,
     .              yorgae0,zorgae0,icouple,iprnsurf,nblelst,iskmax,
     .              jskmax,kskmax,ue,nummem)
c
         do iii = 1,mworki
            iwork(iii) = 0
         end do
         do iii = nstart+1,mwork
            work(iii) = 0.0
         end do
c
         string = ' '
         call cputim(+1,nnodes,string,myhost,myid,mycomm,11)
c
c        determine computational rate (note: ignores start-up overhead)
c        tim(3,2) contains the wall time between calls to cputim
c
c        need a representative value for ncells for each sequence:
c        for non-embeded cases, include only those cells on the finest
c        global level; for embeded cases, include the cells on the finest
c        global level *and* the finest embeded level
c
         lglobal = lfgm-(mseq-iseq)
         if (levelt(iseq).eq.lglobal) then
            ncells = ncell(lglobal)
         else
            ncells = ncell(lglobal) + ncell(levelt(iseq))
         end if
c
c        best attempt at correct wall time for > 24  hour runs:
c        the problem seems to be that tim(3,2), supposedly the wall
c        time, gets reset after 24 hours. have to creative in the case
c        the total time - system+user - is not *quite* 24 hours and 
c        the true wall time is just over 24 hours, so that ndays=0
c        but tim(3,2) is small. in that case use the total time as
c        the wall time.
c
         walltime  = tim(3,1)
         totaltime = tim(1,1) + tim(2,1)
         ndays = int(real(totaltime))/86400
         if (ndays .gt. 0) then
            walltime = walltime + ndays*86400
         else
            walltime = max(walltime,totaltime)
         end if
         if (real(dt).lt.0.) then
            rate(iseq) = 1.e6*walltime/
     .                   float(ncells)/float(ncyc1(iseq))
         else
            rate(iseq) = 1.e6*walltime/
     .                   float(ncells)/float(ntstep)
            ratesub(iseq) = rate(iseq)/float(ncyc1(iseq))
         end if
c
c   meshdef = 1 if flow solution is to be bipassed, and mesh deformation only:
c
      if(meshdef.eq.1) goto 8100
c
c***********************************************************************
c        output for the current sequence level
c***********************************************************************
c
c        output solution (plot3d and printout files)
c
         if (real(dt).lt.0. .or. (real(dt).gt.0..and.movie.eq.0)) then
	    iflag = 0 
	    if (iseq.eq.mseq) iflag = 1
	    if (iseq.lt.mseq .and. ncyc1(iseq+1).eq.0) iflag = 1
	    if (iflag.gt.0) then
	       isklton = 0
	       lhdr    = 1
               iwk1    = 1
               iwk2    = iwk1   + 3*nplots
               iwk3    = iwk2   + maxbl
               mworki1 = mworki - iwk3 + 1
               if (mworki1.le.0) then
                  nou(1) = min(nou(1)+1,ibufdim)
                  write(bou(nou(1),1),'(''stopping...not enough '',
     .                 ''integer work space for subroutine qout'')')
                  call termn8(myid,-1,ibufdim,nbuf,bou,nou)
               end if
               call qout(iseq,lw,lw2,work,nstart,work(nstart+1),nwork,
     .                   nplots,iovrlp,iibg,kkbg,jjbg,ibcg,lbg,
     .                   ibpntsg,qb,lwdat,nbci0,nbcj0,nbck0,
     .                   nbcidim,nbcjdim,nbckdim,jbcinfo,kbcinfo,
     .                   ibcinfo,bcfilei,bcfilej,bcfilek,itrans,
     .                   irotat,idefrm,nblock,levelg,igridg,iviscg,
     .                   jdimg,kdimg,idimg,nblg,clw,ncycmax,nplot3d,
     .                   inpl3d,ip3dsurf,nprint,inpr,iadvance,mycomm,
     .                   myid,myhost,mblk2nd,nou,bou,nbuf,ibufdim,maxbl,
     .                   maxgr,maxseg,iitot,jsg,ksg,isg,jeg,keg,ieg,
     .                   ninter,windex,iindex,nblkpt,intmax,nsub1,maxxe,
     .                   nblk,nbli,limblk,isva,nblon,mxbli,thetay,
     .                   iwork(iwk1),iwork(iwk2),iwork(iwk3),mworki1,
     .                   xorig,yorig,zorig,period_miss,geom_miss,
     .                   epsc0,epsrot,isav_blk,isav_pat,isav_pat_b,
     .                   isav_emb,isav_prd,lbcprd,
     .                   lbcemb,dthetxx,dthetyy,dthetzz,nblcg,
     .                   lfgm,istat2_bl,istat2_pa,istat2_pe,istat2_em,
     .                   istat_size,vormax,ivmax,jvmax,kvmax,bcfiles,
     .                   mxbcfil,iprnsurf,nummem)
c
               do iii = 1,mworki
                  iwork(iii) = 0
               end do
               do iii = nstart+1,mwork
                  work(iii) = 0.0
               end do
c
	    end if
	 end if
c
c        print out control surface data
c

	 iflag = 0
	 if (iseq.eq.mseq) iflag = 1
	 if (iseq.lt.mseq .and. ncyc1(iseq+1).eq.0) iflag = 1
	 if (iflag.gt.0) then
            if (myid.eq.myhost) then
               call csout(iseq,maxbl,maxcs,igridg,levelg,ncs,sx,sy,sz,
     .                    stot,pav,ptav,tav,ttav,xmav,fmdot,cfxp,cfyp,
     .                    cfzp,cfdp,cflp,cftp,cfxv,cfyv,cfzv,cfdv,cflv,
     .                    cftv,cfxmom,cfymom,cfzmom,cfdmom,cflmom,
     .                    cftmom,cfxtot,cfytot,cfztot,cfdtot,cfltot,
     .                    cfttot,icsinfo)
            end if
	 end if
c
c        print out force and moment data
c
         call forceout(iseq,maxbl,maxgr,maxseg,nblock,iforce,
     .                 igridg,nbci0,nbcj0,nbck0,nbcidim,nbcjdim,
     .                 nbckdim,levelg,ibcinfo,jbcinfo,kbcinfo,
     .                 swett,clt,cdt,cxt,cyt,czt,cmxt,cmyt,cmzt,
     .                 cdpt,cdvt,swetw,clw,cdw,cxw,cyw,czw,cmxw,
     .                 cmyw,cmzw,cdpw,cdvw,ncycmax,myhost,myid,
     .                 mycomm,mblk2nd)
c
c        calculate and print out yplus statistics for turbulent flows
c
	 iflag = 0
	 if (iseq.eq.mseq) iflag = 1
	 if (iseq.lt.mseq .and. ncyc1(iseq+1).eq.0) iflag = 1
	 if (ivmx.lt.2) iflag = 0
	 if (iflag.gt.0) then
	    nttuse=max(ntt-1,1)
            iwk1    = 1
            iwk2    = iwk1 + 6*maxbl
            iwk3    = iwk2 + 6*maxbl
            mworki1 = mworki - iwk3
            if (mworki1.lt.0) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),'(''stopping...not enough integer '',
     .              ''work space for subroutine yplusout'')')
               call termn8(myid,-1,ibufdim,nbuf,bou,nou)
            end if
            call yplusout(iseq,lw,lw2,work,nstart,work(nstart+1),nwork,
     .                    clw(nttuse),maxbl,maxgr,maxseg,nblock,
     .                    lwdat,levelg,igridg,jdimg,kdimg,idimg,
     .                    nbci0,nbcj0,nbck0,nbcidim,nbcjdim,
     .                    nbckdim,bcfilei,bcfilej,bcfilek,
     .                    itrans,irotat,idefrm,nblg,ibcinfo,jbcinfo,
     .                    kbcinfo,iadvance,iovrlp,myid,myhost,
     .                    mycomm,mblk2nd,iwork(iwk1),iwork(iwk2),
     .                    iwork(iwk3),nou,bou,nbuf,ibufdim,bcfiles,
     .                    mxbcfil,nummem)
c
            do iii = 1,mworki
               iwork(iii) = 0
            end do
            do iii = nstart+1,mwork
               work(iii) = 0.0
            end do
c
	 end if
c
      end if
c
 8000 continue
 8100 continue
c
c     end of mesh sequencing
c
c***********************************************************************
c    output convergence history
c***********************************************************************
c
#if defined DIST_MPI
      if (myid.eq.myhost) then
#endif
c
c   meshdef = 1 if flow solution is to be bipassed, and mesh deformation only:
c
      if(meshdef.eq.1) goto 8110
c
      call histout(ihstry,rms,clw,cdw,cdpw,cdvw,cxw,cyw,czw,
     .             cmxw,cmyw,cmzw,n_clcd,clcd,nblocks_clcd,blocks_clcd,
     .             chdw,swetw,fmdotw,cfttotw,
     .             cftmomw,cftpw,cftvw,rmstr,nneg,
     .             ncycmax,aehist,aesrfdat,nmds,maxaes,
     .             timekeep,nummem)
#   ifdef CMPLX
c
c     derivative convergence history
c
      call histout_img(ihstry,rms,clw,cdw,cdpw,cdvw,cxw,cyw,czw,
     .                 cmxw,cmyw,cmzw,chdw,swetw,fmdotw,cfttotw,
     .                 cftmomw,cftpw,cftvw,rmstr,nneg,
     .                 ncycmax,aehist,aesrfdat,nmds,maxaes,
     .                 timekeep,nummem)
#   endif
c
      if (abs(movie).gt.0) then
         write(11,161) nframes
  161    format(7h output,i4,
     .   36h frames of plot3d data for animation)
      end if
8110  continue
c
c***********************************************************************
c     output computational rates based on *wall time*
c***********************************************************************
c
      write(11,9994)
      do iseq = 1,mseq
         if (real(dt).lt.0.) then
            write(11,9995) iseq,real(rate(iseq))
         else
            write(11,9996) iseq,real(rate(iseq))
            write(11,9997) real(ratesub(iseq))
         end if
      end do
c
 9994 format(/1x,36hcomputational rate by mesh sequence ,
     .           21h(based on wall time):)
 9995 format(1x,4hiseq,i2,1x,f7.2,28h microseconds/cell/iteration)
 9996 format(1x,4hiseq,i2,1x,f7.2,28h microseconds/cell/time step)
 9997 format(8x,f7.2,31h microseconds/cell/subiteration)
c
#if defined DIST_MPI
      end if
c
#endif
c***********************************************************************
c     output final timings for this run
c***********************************************************************
c
      string = '    timing for complete run - time in seconds     '
      call cputim(-1,nnodes,string,myhost,myid,mycomm,11)
c
      deallocate(rmstr)
      deallocate(nneg)
      return
      end
