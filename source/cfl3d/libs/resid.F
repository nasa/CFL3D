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
      subroutine resid(nbl,ntime,jdim,kdim,idim,q,qj0,qk0,qi0,sj,sk,si,
     .           vol,dtj,x,y,z,vist3d,snj0,snk0,sni0,res,wj0,wk0,wi0,
     .           vmuk,vmuj,vmui,wk,nwork,isf,iwfa,wfa,delt,blank,iover,
     .           nblendg,nblstat,nblstag,xib,sig,sqtq,g,
     .           tj0,tk0,ti0,xkb,blnum,vj0,vk0,vi0,bcj,bck,bci,
     .           nt,sumn,negn,ux,xib2,cmuv,volj0,
     .           volk0,voli0,nou,bou,nbuf,ibufdim,myid,mblk2nd,maxbl,
     .           maxseg,nbci0,nbcj0,nbck0,nbcidim,nbcjdim,nbckdim,
     .           ibcinfo,jbcinfo,kbcinfo,vormax,ivmax,jvmax,kvmax,
     .           idefrm,iadvance,qavg,nummem)
#   ifdef CMPLX
#   else
      use module_kwstm, only: kws_main,kws_dump_movie
#   endif
c
c     $Id$
c
c***********************************************************************
c     Purpose: Compute the residual contributions to the right-hand-side.
c     Note: viscous flux routines called even if block is not advanced
c     in order to compute and store the vmui/j/k arrays (though these
c     routines are then exited as soon as those arrays are computed)
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension wk(nwork)
      dimension wk0(idim*jdim*22),wi0(kdim*jdim*22),wj0(kdim*idim*22)
      dimension x(jdim*kdim*idim),y(jdim*kdim*idim),z(jdim*kdim*idim)
c     dimension vmu(jdim-1,idim-1),vmuj1(kdim-1,idim-1)
      dimension vmuk(jdim-1,idim-1,2),vmuj(kdim-1,idim-1,2),
     .          vmui(jdim-1,kdim-1,2)
      dimension snj0(jdim-1,kdim-1,idim-1),snk0(jdim-1,kdim-1,idim-1),
     .          sni0(jdim-1,kdim-1,idim-1),blank(jdim,kdim,idim)
      dimension q(jdim,kdim,idim,5),   qj0(kdim*(idim-1),5,4),
     .          qk0(jdim*(idim-1),5,4),qi0(jdim*kdim,5,4)
      dimension si(jdim*kdim*idim,5),sj(jdim*kdim*(idim-1),5),
     .          sk(jdim*kdim*(idim-1),5)
      dimension res(jdim*kdim*(idim-1),5),vol(jdim*kdim*(idim-1)),
     .          dtj(jdim*kdim*(idim-1)),vist3d(jdim,kdim,idim)
      dimension xib(jdim,kdim,idim,nummem),sig(jdim-1,idim-1),
     .          sqtq(jdim-1,idim-1),g(-1:jdim+1,-1:idim+1),
     .          tj0(kdim,idim-1,nummem,4),tk0(jdim,idim-1,nummem,4),
     .          ti0(jdim,kdim,nummem,4),xkb(jdim-1,kdim-1,idim-1),
     .          vj0(kdim,idim-1,1,4),vk0(jdim,idim-1,1,4),
     .          vi0(jdim,kdim,1,4),blnum(jdim-1,kdim-1,idim-1)
      dimension volj0(kdim,idim-1,4),volk0(jdim,idim-1,4),
     .          voli0(jdim,kdim,4)
      dimension wfa(1),iwfa(1)
      dimension jbctyp(2),kbctyp(2),ibctyp(2)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
      dimension ux(jdim-1,kdim-1,idim-1,9),xib2(jdim,kdim,idim,2*nummem)
      dimension cmuv(jdim-1,kdim-1,idim-1)
      dimension nbci0(maxbl),nbcidim(maxbl),nbcj0(maxbl),nbcjdim(maxbl),
     .          nbck0(maxbl),nbckdim(maxbl),ibcinfo(maxbl,maxseg,7,2),
     .          jbcinfo(maxbl,maxseg,7,2),kbcinfo(maxbl,maxseg,7,2)
      dimension vormax(maxbl),ivmax(maxbl),jvmax(maxbl),kvmax(maxbl)
      dimension mblk2nd(maxbl),idefrm(maxbl),iadvance(maxbl)
      dimension qavg(jdim,kdim,idim,5)
      dimension sumn(nummem),negn(nummem)
c
      common /avgdata/ xnumavg,iteravg,xnumavg2,ipertavg,iclcd,isubit_r
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /cfl/ dt0,dtold
      common /fvfds/ rkap0(3),ifds(3)
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /nfablk/ nfajki(3)
      common /reyue/ reue,tinf,ivisc(3)
      common /sklton/ isklton
      common /twod/ i2d
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /wrbl/ nwrest
      common /account/ iaccnt,ioutsub
      common /noninertial/ xcentrot,ycentrot,zcentrot,xrotrate,
     .                     yrotrate,zrotrate,noninflag
      common /konew/ ikoprod,isstdenom,pklimterm,ibeta8kzeta,i_bsl,
     .        keepambient,re_thetat0,i_wilcox06,i_wilcox06_chiw,
     .        i_turbprod_kterm,i_catris_kw,prod2d3dtrace,
     .        i_compress_correct,isstsf,i_wilcox98,i_wilcox98_chiw,
     .        isst2003
      common /maxiv/ ivmx
      common /turbconv/ cflturb(7),edvislim,iturbprod,nsubturb,nfreeze,
     .                  iwarneddy,itime2read,itaturb,tur1cut,tur2cut,
     .                  iturbord,tur1cutlev,tur2cutlev
      common /fullns/ ifullns
      common /curvat/ isarc2d,sarccr3,ieasmcc2d,isstrc,sstrc_crc,
     .        isar,crot,isarc3d
      common /sourceterm/ xdir_only_source
      common /mms/ iexact_trunc,iexact_disc,iexact_ring
      common /des/ cdes,ides,cddes
      common /memry/ lowmem_ux
      common /constit/ i_nonlin,c_nonlin,snonlin_lim,i_tauijs,i_qcr2000,
     .                 i_qcr2013,i_qcr2013v
      common /reystressmodel/ issglrrw2012,i_sas_rsm,i_yapterm
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
c     check storage
c
      maxmem = max(jdim*kdim*35,jdim*idim*35)
      if (nwork.lt.maxmem) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),'('' Stopping.  Insufficient memory in'',
     .    '' resid.  You must increase mwork.'')')
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
c
      if (icyc.ge.nitfo+1 .and. level.ge.lglobal) then
c
c     finest global grids or embedded  grids 
c
         rkap(1) = rkap0(1)
         rkap(2) = rkap0(2)
         rkap(3) = rkap0(3)
      else  if (level.ge.lglobal) then
         rkap(1) = -3.
         rkap(2) = -3.
         rkap(3) = -3.
      else
c
c      coarser grids multigrid
c
         rkap(1) = rkap0(1)
         rkap(2) = rkap0(2)
         rkap(3) = rkap0(3)
      end if
c
c     time step  kode.ne.-1
c
      if (kode.eq.-1 .or. ntime.gt.1) go to 130
c
      if (isklton.eq.1 .and.level.eq.levt) then
         if (real(dt).gt.0.0) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2050) real(dt),icyc
         end if
         if (real(dt).lt.0.0) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),2051) abs(real(dt)),icyc
         end if
 2050    format(25h initial time step, icyc=,e12.3,i5)
 2051    format(26h initial cfl number, icyc=,e12.3,i5)
      end if
c
c     local time step: ramp initial CFL number to a value fmax times larger.
c     The time step distribution is recalculated for each cycle.
c
      dttol = 1.e-6
      if (iflagts.gt.0 .and. real(dt).lt.0 .and. icyc.gt.1) then 
        if (nt.eq.1 .and. abs(real(dt)) .lt. 
     .      real(fmax)*abs(real(dt0))-real(dttol)) then
c          ramp CFL on first block on top level
           if (level.eq.levt .and. nbl.eq.nblstat) then
              t2   = 1.e0/float(iflagts)
              fact = fmax**t2
              dt   = fact*dt
c             nou(1) = min(nou(1)+1,ibufdim)
c             write(bou(nou(1),1),*)'cycle, cfl = ',icyc,real(dt)
           end if
        end if
      end if
c
c     global time step:ramp initial time step to a value fmax times larger.
c     The CFL number distribution is recalculated for each time step.
c
      dttol = 1.e-6
      if (iflagts.gt.0 .and. real(dt).gt.0 .and. nt.gt.1) then
        if (icyc.eq.1 .and. abs(real(dt)) .lt. 
     .      real(fmax)*abs(real(dt0))-real(dttol)) then
c          ramp time step on first block on top level
           if (level.eq.levt .and. nbl.eq.nblstat) then
              t2   = 1.e0/float(iflagts)
              fact = fmax**t2
              dt   = fact*dt
c             nou(1) = min(nou(1)+1,ibufdim)
c             write(bou(nou(1),1),*)'time step, dt = ',nt,real(dt)
           end if
        end if
      end if
c
c     calculate time step/CFL based on CFL/timestep
      iout = 0
      if (isklton.eq.1) iout=1
      iterm = max( ivisc(1) , ivisc(2) , ivisc(3) )
c
c     don't call ctime1 during time-accurate subiterations if
c     only physical time terms are used...actually really only
c     need to call ctime1 once in that case
c
      if (real(dt).gt.0 .and. ita.gt.0 .and. icyc.gt.1) go to 130
c
      call ctime1(nbl,jdim,kdim,idim,q,vol,sj,sk,si,dtj,wk,delt,
     .           vist3d,iterm,dtmin,iout,ntime,nou,bou,nbuf,ibufdim,
     .           idefrm(nbl))
c
c     zero out time step at hole points for chimera scheme
c
      if (iover.eq.1) then
         call dthole(jdim,kdim,idim,dtj,vol,blank,dtmin,
     .               nou,bou,nbuf,ibufdim)
      end if
c
      if (iflagts.gt.0 .and. abs(real(fmax)).gt.0. .and.
     .    real(dt).lt.0. .and. myid.le.1) then
      if (((icyc.eq.iflagts+1).or.(icyc.eq.ncyc.and.iflagts+1.gt.ncyc))
     .   .and.level.eq.levt.and.nbl.eq.nblstat.and.nt.eq.1)then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),*)
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),3051) abs(real(dt)),icyc
 3051    format(24h final cfl number, icyc=,e12.3,i5)
         if (icyc.ne.ncyc) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)
         end if
      end if
      end if
c
      if (iflagts.gt.0 .and. abs(real(fmax)).gt.0. .and.
     .    real(dt).gt.0. .and. myid.le.1) then
      if (((nt.eq.iflagts+1).or.(nt.eq.ntstep.and.iflagts+1.gt.ntstep))
     .   .and.level.eq.levt.and.nbl.eq.nblstat.and.icyc.eq.1)then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),*)
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),3050) real(dt),nt
 3050    format(21h final dt, time step=,e12.3,i5)
         if (nt.ne.ntstep) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),*)
         end if
      end if
      end if
c
  130 continue
c
c      turbulent viscosity - finest grids only
c
      if (iadvance(nbl) .ge. 0) then
c
      if (ivisc(3).gt.1 .or. ivisc(2).gt.1 .or.ivisc(1).gt.1) then
      icall = 0
      if (icyc.eq.(icyc/1)*1 .and. level.ge.lglobal .and. ntime.le.nit)
     .    icall = 1
      if (lowmem_ux .eq. 1) then
        if (ivmx.eq.8 .or. ivmx.eq.9 .or. ivmx.ge.11 .or. 
     .      ivmx.eq.16 .or. ivmx.eq.72 .or.
     .    ((ivmx.eq.6 .or. ivmx.eq.7 .or. ivmx.eq.10).and. 
     .    ikoprod.eq.1) .or. 
     .    (ivmx.eq.7 .and. isstdenom.eq.1) .or.
     .    (ivmx.eq.7 .and. isst2003.eq.1) .or.
     .    ((ivmx.eq.6 .or. ivmx.eq.7) .and. isstrc.gt.0) .or.
     .    ((ivmx.eq.6 .or. ivmx.eq.7) .and. isstsf.eq.1) .or.
     .    ((ivmx.eq.5) .and. isarc2d.eq.1) .or. 
     .    ((ivmx.eq.5) .and. isarc3d.eq.1) .or. 
     .    ((ivmx.eq.5) .and. isar.eq.1) .or. 
     .    (ivmx.eq.6 .and. i_wilcox06.eq.1) .or.
     .    (ivmx.eq.6 .and. i_wilcox98.eq.1) .or.
     .    ides.ge.2 .or. i_nonlin.ne.0 .or. i_tauijs.ne.0) then
          call delv(jdim,kdim,idim,q,sj,sk,si,vol,ux,wk,blank,iover,
     .          qj0,qk0,qi0,bcj,bck,bci,nbl,volj0,volk0,voli0,
     .          maxbl,vormax,ivmax,jvmax,kvmax)
        end if
      else
          call delv(jdim,kdim,idim,q,sj,sk,si,vol,ux,wk,blank,iover,
     .          qj0,qk0,qi0,bcj,bck,bci,nbl,volj0,volk0,voli0,
     .          maxbl,vormax,ivmax,jvmax,kvmax)
      end if
      if (icall.gt.0) then
c
c      vorticity magnitude
c
      if (isklton.eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),*)
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),49)
      end if
   49 format(1x,29hcomputing vorticity magnitude)
c
      ipw = 0
      if (real(dt).lt.0.0) then
         if (icyc/nwrest*nwrest.eq.icyc .or. icyc.eq.ncyc) then
            ipw = 1
         end if
      else
         if (nt/nwrest*nwrest.eq.nt .or. nt.eq.ntstep) then
            if (icyc.eq.ioutsub) ipw = 1
         end if
      end if
      call wmag(jdim,kdim,idim,q,sj,sk,si,vol,res,res(1,2),wk,ipw,
     .          blank,iover,qj0,qk0,qi0,bcj,bck,bci,nbl,volj0,
     .          volk0,voli0,vormax,ivmax,jvmax,kvmax,maxbl)
c
c      wall/wake turbulence model
c
      if (isklton.eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),*)
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),89)
      end if
   89 format(1x,27hevaluating turbulence model)
c
      if (ivmx .eq. 2) then
c   Baldwin-Lomax
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
      nroom=nwork-(iwk9+inmx)
      if(nroom .lt. 0.) then
       nou(1) = min(nou(1)+1,ibufdim)
       write(bou(nou(1),1),'('' not enough memory for B-L turb'',
     . '' model.'')')
       nou(1) = min(nou(1)+1,ibufdim)
       write(bou(nou(1),1),'('' nroom='',i12)') nroom
       call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
      nblt = nbl
      ipw = 0
      if (real(dt).lt.0.) then
         if (icyc.eq.ncyc) ipw = 1
         if (icyc.eq.ncyc .and. icyc.le.3) ipw = 2
      else 
         if (icyc.eq.ioutsub) then
            if (nt.eq.ntstep) ipw = 1
            if (nt.eq.ntstep .and. ntstep.le.3) ipw = 2
         end if
      end if
      call blomax(jdim,kdim,idim,q,qi0,qj0,qk0,res,snj0,snk0,sni0,
     .            blnum,xkb,blnum,
     .            vist3d,res(1,2),ipw,inmx,wk(iwk1),wk(iwk2),wk(iwk3),
     .            wk(iwk4),wk(iwk5),wk(iwk6),wk(iwk7),wk(iwk8),wk(iwk9),
     .            nblt,x,y,z,blank,iover,bci,bcj,bck,
     .            nou,bou,nbuf,ibufdim,myid,mblk2nd,maxbl)
c
c   Johnson-King no longer available due to unsuitability for 
c   multi-block grids
c
      else if (ivmx .eq. 4) then
c   Baldwin-Barth
      iwk1=1
      if (iturbord .eq. 1) then
        iwk5=iwk1+(jdim+1)*(kdim+1)*(idim+1-(2*i2d))
        iex=0
        iex2=-i2d
      else
        iwk5=iwk1+(jdim+3)*(kdim+3)*(idim+3-(4*i2d))
        iex=1
        iex2=1-(2*i2d)
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
      nroom=nwork-(iwk19+inmx)
      if(nroom .lt. 0.) then
       nou(1) = min(nou(1)+1,ibufdim)
       write(bou(nou(1),1),'('' not enough memory for B-B turb'',
     . '' model.'')')
       nou(1) = min(nou(1)+1,ibufdim)
       write(bou(nou(1),1),'('' nroom='',i12)') nroom
       call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
      call barth3d(jdim,kdim,idim,q,sj,sk,si,vol,dtj,x,y,z,vist3d,res,
     + snk0,snj0,xib,xkb,wk(iwk1),res(1,2),res(1,3),res(1,4),res(1,5),
     + wk(iwk5),
     + wk(iwk6),wk(iwk7),wk(iwk8),wk(iwk9),wk(iwk10),wk(iwk11),
     + wk(iwk12),wk(iwk13),wk(iwk14),wk(iwk15),wk(iwk16),wk(iwk17),
     + wk(iwk18),wk(iwk19),ntime,tj0,tk0,ti0,nbl,blnum,blank,iover,
     + sumn(1),sumn(2),negn(1),negn(2),xib2,volj0,volk0,voli0,nou,bou,
     + nbuf,ibufdim,iex,iex2,nummem)
      else if (ivmx .eq. 5) then
c   Spalart
      iwk1=1
      if (iturbord .eq. 1) then
        iwk3=iwk1+(jdim+1)*(kdim+1)*(idim+1-(2*i2d))
        iex=0
        iex2=-i2d
      else
        iwk3=iwk1+(jdim+3)*(kdim+3)*(idim+3-(4*i2d))
        iex=1
        iex2=1-(2*i2d)
      end if
      iwk4=iwk3+(jdim+1)*(kdim+1)*(idim+1-(2*i2d))
      iex3=-i2d
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
        isarcnum=4
        inmx=(jdim+1)*(kdim+1)*(idim-1)*isarcnum
      else if (isarc3d .eq. 1) then
        isarcnum=6
        inmx=(jdim+1)*(kdim+1)*(idim-1)*isarcnum
      else
        isarcnum=1
        inmx=0
      end if
      nroom=nwork-(iwk19+inmx)
      if(nroom .lt. 0.) then
       nou(1) = min(nou(1)+1,ibufdim)
       write(bou(nou(1),1),'('' not enough memory for Spalart'',
     . '' turb model.'')')
       nou(1) = min(nou(1)+1,ibufdim)
       write(bou(nou(1),1),'('' nroom='',i12)') nroom
       call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
      call spalart(jdim,kdim,idim,q,sj,sk,si,vol,dtj,x,y,z,vist3d,
     + res,snk0,xib,wk(iwk1),res(1,2),res(1,3),wk(iwk3),wk(iwk4),
     + wk(iwk5),wk(iwk6),wk(iwk7),wk(iwk8),wk(iwk9),wk(iwk10),wk(iwk11),
     + wk(iwk12),wk(iwk13),wk(iwk14),wk(iwk15),wk(iwk16),wk(iwk17),
     + wk(iwk18),ntime,tj0,tk0,ti0,nbl,qj0,qk0,qi0,blank,iover,sumn(1),
     + sumn(2),negn(1),negn(2),xib2,volj0,volk0,voli0,nou,bou,nbuf,
     + ibufdim,iex,iex2,iex3,bcj,bck,bci,ux,wk(iwk19),nummem,res(1,4),
     + isarcnum)
      else if (ivmx .ge. 6 .and. ivmx .le. 16) then
c   two-equation turb model
      iwk1=1
      if (iturbord .eq. 1) then
        iwk4=iwk1+(jdim+1)*(kdim+1)*(idim+1-(2*i2d))*2
        iex=0
        iex2=-i2d
      else
        iwk4=iwk1+(jdim+3)*(kdim+3)*(idim+3-(4*i2d))*2
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
      inmx=(jdim-1)*(kdim-1)*(idim-1)*2
      iwk33=iwk32+inmx
      inmx=(jdim+1)*(kdim+1)*(idim+1-(2*i2d))
c allocate for gradients of vorticity and turbulence - k-zeta only
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
c allocate for easmcc2d and certain isstrc
      if (ieasmcc2d .eq. 1) then
        iccnum=4
        inmx=(jdim+1)*(kdim+1)*(idim-1)*4
      else if (isstrc .eq. 2) then
        iccnum=6
        inmx=(jdim+1)*(kdim+1)*(idim-1)*6
      else if (ivmx .eq. 16) then
        iccnum=4
        inmx=(jdim+1)*(kdim+1)*(idim-1)*4
      else
        iccnum=1
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
      nroom=nwork-(iwk82+inmx)
      if(nroom .lt. 0.) then
       nou(1) = min(nou(1)+1,ibufdim)
       write(bou(nou(1),1),'('' not enough memory for 2-eqn turb'',
     . '' model.'')')
       nou(1) = min(nou(1)+1,ibufdim)
       write(bou(nou(1),1),'('' nroom='',i12)') nroom
       call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
      call twoeqn(jdim,kdim,idim,q,sj,sk,si,vol,dtj,x,y,z,vist3d,
     + res,snk0,xib,wk(iwk1),res(1,2),res(1,3),res(1,4),wk(iwk4),
     + wk(iwk5),
     + wk(iwk6),wk(iwk7),wk(iwk8),wk(iwk9),wk(iwk10),wk(iwk11),
     + wk(iwk12),wk(iwk13),wk(iwk14),wk(iwk15),wk(iwk16),wk(iwk17),
     + wk(iwk18),wk(iwk19),wk(iwk20),wk(iwk21),wk(iwk22),wk(iwk23),
     + wk(iwk24),wk(iwk25),wk(iwk26),wk(iwk27),wk(iwk28),wk(iwk29),
     + wk(iwk30),wk(iwk31),ntime,tj0,tk0,ti0,nbl,qj0,qk0,qi0,
     + vj0,vk0,vi0,blank,iover,sumn(1),sumn(2),negn(1),negn(2),ux,
     + wk(iwk32),xib2,wk(iwk33),cmuv,bcj,bck,bci,nbci0,nbcidim,
     + nbcj0,nbcjdim,nbck0,nbckdim,ibcinfo,jbcinfo,kbcinfo,
     + maxbl,maxseg,volj0,volk0,voli0,nou,bou,nbuf,ibufdim,iex,iex2,
     + iex3,
     + wk(iwk34),wk(iwk35),wk(iwk36),wk(iwk37),wk(iwk38),wk(iwk39),
     + wk(iwk40),wk(iwk41),wk(iwk42),wk(iwk43),wk(iwk44),wk(iwk45),
     + wk(iwk46),wk(iwk47),wk(iwk48),wk(iwk49),wk(iwk50),wk(iwk51),
     + wk(iwk52),wk(iwk53),wk(iwk54),
     + wk(iwk55),wk(iwk56),wk(iwk57),wk(iwk58),wk(iwk59),wk(iwk60),
     + wk(iwk61),wk(iwk62),wk(iwk63),wk(iwk64),wk(iwk65),wk(iwk66),
     + wk(iwk67),wk(iwk68),wk(iwk69),wk(iwk70),wk(iwk71),wk(iwk72),
     + wk(iwk73),wk(iwk74),wk(iwk75),
     + wk(iwk76),wk(iwk77),wk(iwk78),wk(iwk79),wk(iwk80),wk(iwk81),
     + wk(iwk82),nummem,iccnum)
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
      nroom=nwork-(iwk45+inmx)
      if(nroom .lt. 0.) then
       nou(1) = min(nou(1)+1,ibufdim)
       write(bou(nou(1),1),'('' not enough memory for 3-eqn turb'',
     . '' model.'')')
       nou(1) = min(nou(1)+1,ibufdim)
       write(bou(nou(1),1),'('' nroom='',i12)') nroom
       call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
      call threeeqn(jdim,kdim,idim,q,sj,sk,si,vol,dtj,x,y,z,vist3d,
     + res,snk0,xib,wk(iwk1),res(1,2),res(1,3),res(1,4),wk(iwk4),
     + wk(iwk5),
     + wk(iwk6),wk(iwk7),wk(iwk8),wk(iwk9),wk(iwk10),wk(iwk11),
     + wk(iwk12),wk(iwk13),wk(iwk14),wk(iwk15),wk(iwk16),wk(iwk17),
     + wk(iwk18),wk(iwk19),wk(iwk20),wk(iwk21),wk(iwk22),wk(iwk23),
     + wk(iwk24),wk(iwk25),wk(iwk26),wk(iwk27),wk(iwk28),wk(iwk29),
     + wk(iwk30),wk(iwk31),ntime,tj0,tk0,ti0,nbl,qj0,qk0,qi0,
     + vj0,vk0,vi0,blank,iover,sumn(1),sumn(2),sumn(3),negn(1),negn(2),
     + negn(3),ux,wk(iwk32),
     + xib2,wk(iwk33),bcj,bck,bci,nbci0,nbcidim,
     + nbcj0,nbcjdim,nbck0,nbckdim,ibcinfo,jbcinfo,kbcinfo,
     + maxbl,maxseg,volj0,volk0,voli0,nou,bou,nbuf,ibufdim,iex,iex2,
     + iex3,wk(iwk34),wk(iwk35),wk(iwk36),wk(iwk37),
     + wk(iwk38),wk(iwk39),wk(iwk40),wk(iwk41),
     + wk(iwk42),wk(iwk43),wk(iwk44),wk(iwk45),nummem)
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
      nroom=nwork-(iwk59+inmx)
      if(nroom .lt. 0.) then
       nou(1) = min(nou(1)+1,ibufdim)
       write(bou(nou(1),1),'('' not enough memory for 4-eqn turb'',
     . '' model.'')')
       nou(1) = min(nou(1)+1,ibufdim)
       write(bou(nou(1),1),'('' nroom='',i12)') nroom
       call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
      call foureqn(jdim,kdim,idim,q,sj,sk,si,vol,dtj,x,y,z,vist3d,
     + res,snk0,xib,wk(iwk1),res(1,2),res(1,3),res(1,4),res(1,5),
     + wk(iwk4),wk(iwk5),
     + wk(iwk6),wk(iwk7),wk(iwk8),wk(iwk9),wk(iwk10),wk(iwk11),
     + wk(iwk12),wk(iwk13),wk(iwk14),wk(iwk15),wk(iwk16),wk(iwk17),
     + wk(iwk18),wk(iwk19),wk(iwk20),wk(iwk21),wk(iwk22),wk(iwk23),
     + wk(iwk24),wk(iwk25),wk(iwk26),wk(iwk27),wk(iwk28),wk(iwk29),
     + wk(iwk30),wk(iwk31),ntime,tj0,tk0,ti0,nbl,qj0,qk0,qi0,
     + vj0,vk0,vi0,blank,iover,sumn(1),sumn(2),sumn(3),sumn(4),
     + negn(1),negn(2),negn(3),negn(4),ux,wk(iwk32),
     + xib2,wk(iwk33),bcj,bck,bci,nbci0,nbcidim,
     + nbcj0,nbcjdim,nbck0,nbckdim,ibcinfo,jbcinfo,kbcinfo,
     + maxbl,maxseg,volj0,volk0,voli0,nou,bou,nbuf,ibufdim,iex,iex2,
     + iex3,wk(iwk34),wk(iwk35),wk(iwk36),wk(iwk37),
     + wk(iwk38),wk(iwk39),wk(iwk40),wk(iwk41),
     + wk(iwk42),wk(iwk43),wk(iwk44),wk(iwk45),
     + wk(iwk46),wk(iwk47),wk(iwk48),wk(iwk49),
     + wk(iwk50),wk(iwk51),wk(iwk52),wk(iwk53),
     + wk(iwk54),wk(iwk55),wk(iwk56),wk(iwk57),
     + wk(iwk58),wk(iwk59),nummem)
      else if (ivmx .eq. 72) then ! full stress model, k-omega based stress model
#   ifdef CMPLX
c          RSM not working in complex mode
#   else
         if (isklton.eq.1) then
            if (issglrrw2012 .eq. 1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'('' Computing turbulence using '',
     +        ''Eisfeld SSG/LRR-RSM-w2012 7-eqn RSM, block='',i5)') nbl
            else if (issglrrw2012 .eq. 2) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'('' Computing turbulence using '',
     +        ''Eisfeld SSG/LRR-RSM-w2012 7-eqn RSM with F1=1'',
     +        '', block='',i5)') nbl
            else if (issglrrw2012 .eq. 3) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'('' Computing turbulence using '',
     +        ''Eisfeld SSG/LRR-RSM-w2012 7-eqn RSM with simple '',
     +        ''diffusion model, block='',i5)') nbl
            else if (issglrrw2012 .eq. 4) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'('' Computing turbulence using '',
     +        ''Eisfeld SSG/LRR-RSM-w2012 7-eqn RSM with F1=1 and '',
     +        ''simple diffusion model, block='',i5)') nbl
            else if (issglrrw2012 .eq. 5) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'('' Computing turbulence using '',
     +        ''Eisfeld SSG/LRR-RSM-w2012 7-eqn RSM with Wilcox '',
     +        ''simple diffusion (non-blended), block='',i5)') nbl
            else if (issglrrw2012 .eq. 6) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'('' Computing turbulence using '',
     +        ''Eisfeld SSG/LRR-RSM-w2012-g 7-eqn RSM,'',
     +        '' block='',i5)') nbl
            else
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'('' Computing turbulence using '',
     +        ''Wilcox Stress-Omega 7-eqn RSM, block='',i5)') nbl
            end if
            if (i_sas_rsm .eq. 1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''   ...adding SAS-like term '',
     +        ''(AIAA-2014-0586)'')')
            else if (i_sas_rsm .eq. -1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''   ...subtracting SAS-like term '',
     +        ''(AIAA-2014-0586)'')')
            end if
            if (i_yapterm .eq. 1) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''   ...adding Yap term'')')
            end if
         end if
         if (icyc .le. nfreeze) then
           if (isklton .gt. 0) then
              nss=min(ncyc,nfreeze)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),*)
              nou(1) = min(nou(1)+1,ibufdim)
              write(bou(nou(1),1),'('' RSM turbulence model is '',
     +        ''frozen for '',i5,'' iterations or subits'')') nss
           end if
           do n=1,nummem
             sumn(n) = 0.
             negn(n) = 0
           enddo
         else
           call kws_main(jdim,kdim,idim,myid,myhost,nummem, 
     $          tj0,tk0,ti0,xib,
     $          ux,
     $          sj,sk,si,vol,volj0,volk0,voli0,q,qj0,qk0,qi0,dtj,vist3d,
     $          level,icyc,sumn,negn,xib2,snk0,x,y,z,nbl,issglrrw2012,
     $          i_sas_rsm,i_yapterm)
         end if
#   endif
      else if (ivmx .eq. 25) then
c   LES-type model and diagnostics
        call lesdiag(myid,jdim,kdim,idim,q,ux,vist3d,vol,si,sj,sk,
     +    res,snk0,snj0,xib,xkb,blnum,nou,bou,nbuf,ibufdim,nbl,nummem,
     +    x,y,z)
c   *********
      end if
      end if
      end if
c
c      residuals  res = r(q)
c
      if (isklton.gt.0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),*)
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),7) nbl
         if (ifullns .ne. 0 .or. ivmx .ge. 70) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),777)
         end if
      end if
    7 format(1x,29hcomputing residuals for block,i6)
  777 format(1x,37hFull N-S terms are ON - meanflow only)
c
      end if  ! iadvance .ge. 0
c
      iwk4   = 1
      iwk5   = 1
      iwk6   = 1
c
c      zero residuals
c
      n    = jdim*kdim
      nplq = min(idim1,999000/n)
      npl  = nplq
      do 6 i=1,idim1,nplq
      if (i+npl-1.gt.idim1) npl = idim1-i+1
      nnpl = n*npl
      ist  = (i-1)*n+1
      do 6 l=1,5
cdir$ ivdep
      do 6 izz=1,nnpl
      res(izz+ist-1,l) = 0.0e0
    6 continue 
c
c     imult: multi-plane vectorization flag
c            = 0 single plane at a time
c            > 0 multiple planes at a time
c
      imult = 1
c
c     residuals   J direction
c
      if (iadvance(nbl) .ge. 0) then
      nv = 35
      if (idefrm(nbl) .gt. 0) nv = 41
      nvtq  = min(999000,nwork/nv)
      n     = jdim*kdim
      nplq  = min(idim1,nvtq/n)
      if (imult.eq.0) nplq = 1
      npl   = nplq
      niwfac = 1
      do 200 i=1,idim1,nplq
      if (i+npl-1.gt.idim1) npl = idim1-i+1
      nvtq  = npl*jdim*kdim
      jbctyp(1) = jbcinfo(nbl,1,1,1)
      jbctyp(2) = jbcinfo(nbl,1,1,2)
      call gfluxr(i,npl,rkap(2),jdim,kdim,idim,res,q,qj0,sj,wk,nvtq,
     .            nv,nfajki(1),wfa,iwfa(niwfac),jbctyp,isf,nbl,bcj,
     .            nou,bou,nbuf,ibufdim,myid,mblk2nd,maxbl,idefrm(nbl))
  200 continue
      end if
c
c     residuals   J direction - viscous terms
c
      if (ivisc(2).gt.0) then
         nvtq = min(999000,nwork/32)
         n    = jdim*kdim
         nplq = min(idim1,nvtq/n)
         if (imult.eq.0) nplq = 1
         npl = nplq
         do 201 i=1,idim1,nplq
         if (i+npl-1.gt.idim1) npl = idim1-i+1
         nvtq = npl*jdim*kdim
         call gfluxv(i,npl,jdim,kdim,idim,res,q,qj0,sj,vol,wk,
     .               nvtq,wj0,vist3d,vmuj,vj0,bcj,xib,tj0,cmuv,
     .               volj0,nou,bou,nbuf,ibufdim,iadvance(nbl),nummem,ux)
c
c-- full viscous terms
c
         if (ifullns .ne. 0) then
         call gfluxv1(i,npl,jdim,kdim,idim,res,q,qj0,qk0,qi0,
     .                sj,sk,si,vol,wk,nvtq,wj0,vist3d,vj0,
     .                bcj,bck,bci,volj0,
     .                nou,bou,nbuf,ibufdim,iadvance(nbl))
         end if
  201    continue
      end if
c
      if (icyc.eq.1 .and. iadvance(nbl).ge.0)
     .   call l2norm(nbl,0,resd,+1,jdim,kdim,idim,res,vol)
c
c     residuals   K direction
c
      if (iadvance(nbl) .ge. 0) then
      nv = 35
      if (idefrm(nbl) .gt. 0) nv = 41
      nvtq = min(999000,nwork/nv)
      n    = jdim*kdim
      nplq = min(idim1,nvtq/n)
      if (imult.eq.0) nplq = 1
      npl  = nplq
      if (nfajki(2).gt.0) then
         niwfac = nfajki(1)*7+1
      end if
      do 210 i=1,idim1,nplq
      if (i+npl-1.gt.idim1) npl = idim1-i+1
      nvtq = npl*jdim*kdim
      kbctyp(1) = kbcinfo(nbl,1,1,1)
      kbctyp(2) = kbcinfo(nbl,1,1,2)
      call hfluxr(i,npl,rkap(3),jdim,kdim,idim,res,q,qk0,sk,wk,nvtq, 
     .            nv,nfajki(2),wfa,iwfa(niwfac),kbctyp,isf,nbl,bck,
     .            nou,bou,nbuf,ibufdim,myid,mblk2nd,maxbl,idefrm(nbl))
  210 continue
      end if
c
c     residuals   K direction - viscous terms
c
      if (ivisc(3).gt.0) then
         nvtq = min(999000,nwork/32)
         n    = jdim*kdim
         nplq = min(idim1,nvtq/n)
         if (imult.eq.0) nplq = 1
         npl = nplq
         do 211 i=1,idim1,nplq
         if (i+npl-1.gt.idim1) npl = idim1-i+1
         nvtq = npl*jdim*kdim
         call hfluxv(i,npl,jdim,kdim,idim,0,wk(iwk4),wk(iwk5),wk(iwk6),
     .               res,q,qk0,sk,vol,wk,nvtq,wk0,vist3d,vmuk,vk0,bck,
     .               xib,tk0,cmuv,volk0,nou,bou,nbuf,ibufdim,
     .               iadvance(nbl),nummem,ux)
c
c-- full viscous terms
c
         if (ifullns .ne. 0) then
         call hfluxv1(i,npl,jdim,kdim,idim,res,q,qj0,qk0,qi0,
     .                sj,sk,si,vol,wk,nvtq,wk0,vist3d,vk0,
     .                bcj,bck,bci,volk0,
     .                nou,bou,nbuf,ibufdim,iadvance(nbl))
         end if
  211    continue
      end if
c
      if (icyc.eq.1 .and. iadvance(nbl).ge.0)
     .   call l2norm(nbl,0,resd,+1,jdim,kdim,idim,res,vol)
c
c     residuals in  I direction
c
      if (i2d.eq.0) then
c
         if (iadvance(nbl) .ge. 0) then
         n  = jdim*idim
         nv = 35
         if (idefrm(nbl) .gt. 0) nv = 41
         nvtq = min(999000,nwork/nv)
         nplq = min(kdim,nvtq/n)
         if (imult.eq.0) nplq = 1
         npl  = nplq
         if (nfajki(3).gt.0) then
            niwfac = (nfajki(1)+nfajki(2))*7+1
         end if
         do 300 k=1,kdim1,nplq
         if (k+npl-1.gt.kdim1) npl = kdim1-k+1
         nvtq = nplq*jdim*idim
         ibctyp(1) = ibcinfo(nbl,1,1,1)
         ibctyp(2) = ibcinfo(nbl,1,1,2)
         call ffluxr(k,npl,rkap(1),jdim,kdim,idim,res,q,qi0,si,wk,
     .               nvtq,nv,nfajki(3),wfa,iwfa(niwfac),ibctyp,isf,
     .               nbl,bci,nou,bou,nbuf,ibufdim,myid,mblk2nd,maxbl,
     .               idefrm(nbl))
  300    continue
         end if
c
c        residuals   I direction - viscous terms
c
         if (ivisc(1).gt.0) then
c
            nvtq = min(999000,nwork/32)
            n    = jdim*idim
            nplq = min(kdim1,nvtq/n)
            if (imult.eq.0) nplq = 1
            npl  = nplq
            do 277 k=1,kdim1,nplq
            if (k+npl-1.gt.kdim1) npl = kdim1-k+1
            nvtq = npl*jdim*idim
            call ffluxv(k,npl,jdim,kdim,idim,res,q,qi0,si,vol,wk,
     .                  nvtq,wi0,vist3d,vmui,vi0,bci,xib,ti0,cmuv,
     .                  voli0,nou,bou,nbuf,ibufdim,iadvance(nbl),nummem
     .                  ,ux)
c
c-- full viscous terms
c
            if (ifullns .ne. 0) then
            call ffluxv1(k,npl,jdim,kdim,idim,res,q,qj0,qk0,qi0,
     .                   sj,sk,si,vol,wk,nvtq,wi0,vist3d,vi0,
     .                   bcj,bck,bci,voli0,
     .                   nou,bou,nbuf,ibufdim,iadvance(nbl))
            end if
  277       continue
         end if
c
         if (icyc.eq.1 .and. iadvance(nbl).ge.0)
     .      call l2norm(nbl,0,resd,+1,jdim,kdim,idim,res,vol)
c
      end if
c
c     call resnonin to add rotating noninertital source term to res
c
      if (noninflag.gt.0) then
        call resnonin(nbl,jdim,kdim,idim,q,x,y,z,sj,sk,si,vol,res,
     .                nou,bou,nbuf,ibufdim)
      endif
c
c     add source term in x-dir for problems with inflow/outflow bdys
c     that are periodic with each other... this source term drives the problem
c     (there must also be viscous walls in the problem - the source term
c     adjusts the flow until the source balances the viscous losses)
c     grid must be cartesian
c
      if (real(ccabs(xdir_only_source)) .gt. 1.e-12) then
        if (isklton.eq.1) then
          nou(1) = min(nou(1)+1,ibufdim)
          write(bou(nou(1),1),'('' Adding x-dir source term to eqns'')')
        end if
        izz=0
        do i=1,idim-1
        do k=1,kdim
        do j=1,jdim
          izz=izz+1
          fsrce = -xdir_only_source*vol(izz)
          res(izz,2) = res(izz,2) + fsrce
          res(izz,5) = res(izz,5) + fsrce*q(j,k,i,2)
        enddo
        enddo
        enddo
      end if
c
      if (iadvance(nbl) .ge. 0) then
c
c        zero out residuals in edges of arrays (multi-plane vectorization)
c
         do 551 i=1,idim1
         jkv = (i-1)*jdim*kdim
         do 551 l=1,5
         jj    = 1-jdim
         kdww2 = kdim-2
         do 4588 ii=1,kdww2
         jj    = jj+jdim
 4588    res(jkv+jdim+jj-1,l) = 0. 
cdir$ ivdep
         do 551 izz=1,jdim+1
         res(izz+jkv+jdim*kdim1-1,l) = 0.e0
  551    continue
c
         if (iover.eq.1) then
            nvtq = min(999000,nwork/35)
            n    = jdim*kdim
            nplq = min(idim1,nvtq/n)
            if (imult.eq.0) nplq = 1
            npl  = nplq
c
            do 219 i=1,idim1,nplq
            if (i+npl-1.gt.idim1) npl = idim1-i+1
            if (isklton.gt.0 .and. i.eq.1) then
               nou(1) = min(nou(1)+1,ibufdim)
               write(bou(nou(1),1),1447) nbl
            end if
            call hole(i,npl,jdim,kdim,idim,res,blank)
  219       continue
         end if
c
      end if
c
 1447 format(50h blanking the rhs for the cells in the hole  block,i5)
c
c   exact solution
c
      if (iexact_trunc .ne. 0 .or. iexact_disc .ne. 0) then
        if (isklton.gt.0 .and. iexact_trunc .ne. 0) then
          nou(1) = min(nou(1)+1,ibufdim)
          if (iexact_trunc .eq. 1) then
            write(bou(nou(1),1),'(''   Using exact MMS source in'',
     +        '' resid routine to evaluate truncation error, MS1'')')
          else if (iexact_trunc .eq. 2) then
            write(bou(nou(1),1),'(''   Using exact MMS source in'',
     +        '' resid routine to evaluate truncation error, MS2'')')
          else if (iexact_trunc .eq. 4) then
            write(bou(nou(1),1),'(''   Using exact MMS source in'',
     +        '' resid routine to evaluate truncation error, MS4'')')
          end if
        else if (isklton.gt.0 .and. iexact_disc .ne. 0) then
          nou(1) = min(nou(1)+1,ibufdim)
          if (iexact_disc .eq. 1) then
            write(bou(nou(1),1),'(''   Using exact MMS source in'',
     +       '' resid routine to evaluate discretization error, MS1'')')
          else if (iexact_disc .eq. 2) then
            write(bou(nou(1),1),'(''   Using exact MMS source in'',
     +       '' resid routine to evaluate discretization error, MS2'')')
          else if (iexact_disc .eq. 4) then
            write(bou(nou(1),1),'(''   Using exact MMS source in'',
     +       '' resid routine to evaluate discretization error, MS4'')')
          end if
        end if
        call exact_flow_force(jdim,kdim,idim,x,y,z,vol,res,snk0,
     +      iexact_trunc,iexact_disc)
        if (iexact_ring .eq. 1) then
          if (isklton.gt.0) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''   Setting res=0 in'',
     +       '' outer ring'')')
          end if
          call zero_resid_ring(jdim,kdim,idim,res,jdim,kdim,idim-1,
     +      5,2,iexact_trunc,iexact_disc)
        end if
c
        if(icyc.eq.ncyc1(1) .or. icyc.eq.ncyc1(2) .or.
     +     icyc.eq.ncyc1(3) .or. icyc.eq.ncyc1(4) .or.
     +     icyc.eq.ncyc1(5)) then
          if(ntime .eq. 1 .and. level.ge.lglobal) then
            call exact_norm(jdim,kdim,idim,x,y,z,q,xib,vol,res,snk0,
     +           iexact_trunc,iexact_disc)
          end if
        end if
      end if
c
      if (ipertavg .eq. 1 .or. ipertavg .eq. 2) then
         return
      end if

      if (iteravg .eq. 1 .or. iteravg .eq. 2) then
      if (level.ge.lglobal .and. ntime.eq.nit) then
        if (real(dt) .lt. 0. .or. (real(dt) .gt. 0. .and. 
     +      icyc .eq. ncyc)) then
c   get iteration-averaged Q values and increment xnumavg
c   note: the qavg values are kept as conserved variables, whereas 
c   q values are primitive (this is done to avoid round-off differences
c   during restarts, because the PLOT3D output is in conserved form)
          do i=1,idim-1
            do j=1,jdim-1
              do k=1,kdim-1
                qavg(j,k,i,1)=(qavg(j,k,i,1)*(xnumavg-1.)+
     +            q(j,k,i,1))/xnumavg
                qavg(j,k,i,2)=(qavg(j,k,i,2)*(xnumavg-1.)+
     +            q(j,k,i,1)*q(j,k,i,2))/xnumavg
                qavg(j,k,i,3)=(qavg(j,k,i,3)*(xnumavg-1.)+
     +            q(j,k,i,1)*q(j,k,i,3))/xnumavg
                qavg(j,k,i,4)=(qavg(j,k,i,4)*(xnumavg-1.)+
     +            q(j,k,i,1)*q(j,k,i,4))/xnumavg
                qavg(j,k,i,5)=(qavg(j,k,i,5)*(xnumavg-1.)+
     +            (q(j,k,i,5)/gm1+0.5*(q(j,k,i,2)**2+
     +            q(j,k,i,3)**2+q(j,k,i,4)**2)/q(j,k,i,1))
     +            )/xnumavg
              enddo
            enddo
          enddo
        end if
      end if
      end if
      return 
      end
