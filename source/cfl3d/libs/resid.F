      subroutine resid(nbl,ntime,jdim,kdim,idim,q,qj0,qk0,qi0,sj,sk,si,
     .           vol,dtj,x,y,z,vist3d,snj0,snk0,sni0,res,wj0,wk0,wi0,
     .           vmuk,vmuj,vmui,wk,nwork,isf,iwfa,wfa,delt,blank,iover,
     .           nblendg,nblstat,nblstag,xib,sig,sqtq,g,
     .           tj0,tk0,ti0,xkb,blnum,vj0,vk0,vi0,bcj,bck,bci,
     .           nt,sumn1,sumn2,negn1,negn2,ux,xib2,cmuv,volj0,
     .           volk0,voli0,nou,bou,nbuf,ibufdim,myid,mblk2nd,maxbl,
     .           maxseg,nbci0,nbcj0,nbck0,nbcidim,nbcjdim,nbckdim,
     .           ibcinfo,jbcinfo,kbcinfo,vormax,ivmax,jvmax,kvmax,
     .           idefrm,iadvance)
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
      dimension q(jdim*kdim*idim,5),   qj0(kdim*(idim-1),5,4),
     .          qk0(jdim*(idim-1),5,4),qi0(jdim*kdim,5,4)
      dimension si(jdim*kdim*idim,5),sj(jdim*kdim*(idim-1),5),
     .          sk(jdim*kdim*(idim-1),5)
      dimension res(jdim*kdim*(idim-1),5),vol(jdim*kdim*(idim-1)),
     .          dtj(jdim*kdim*(idim-1)),vist3d(jdim,kdim,idim)
      dimension xib(jdim,kdim,idim,2),sig(jdim-1,idim-1),
     .          sqtq(jdim-1,idim-1),g(-1:jdim+1,-1:idim+1),
     .          tj0(kdim,idim-1,2,4),tk0(jdim,idim-1,2,4),
     .          ti0(jdim,kdim,2,4),xkb(jdim-1,kdim-1,idim-1),
     .          vj0(kdim,idim-1,1,4),vk0(jdim,idim-1,1,4),
     .          vi0(jdim,kdim,1,4),blnum(jdim-1,kdim-1,idim-1)
      dimension volj0(kdim,idim-1,4),volk0(jdim,idim-1,4),
     .          voli0(jdim,kdim,4)
      dimension wfa(1),iwfa(1)
      dimension jbctyp(2),kbctyp(2),ibctyp(2)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
      dimension ux(jdim-1,kdim-1,idim-1,9),xib2(jdim,kdim,idim,4)
      dimension cmuv(jdim-1,kdim-1,idim-1)
      dimension nbci0(maxbl),nbcidim(maxbl),nbcj0(maxbl),nbcjdim(maxbl),
     .          nbck0(maxbl),nbckdim(maxbl),ibcinfo(maxbl,maxseg,7,2),
     .          jbcinfo(maxbl,maxseg,7,2),kbcinfo(maxbl,maxseg,7,2)
      dimension vormax(maxbl),ivmax(maxbl),jvmax(maxbl),kvmax(maxbl)
      dimension mblk2nd(maxbl),idefrm(maxbl),iadvance(maxbl)
c
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
      common /unst/ time,cfltau,ntstep,ita,iunst
      common /wrbl/ nwrest
      common /account/ iaccnt,ioutsub
      common /noninertial/ xcentrot,ycentrot,zcentrot,xrotrate,
     .                     yrotrate,zrotrate,noninflag
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
      if (ivisc(3).eq.8 .or. ivisc(2).eq.8 .or. ivisc(1).eq.8 .or.
     .    ivisc(3).eq.9 .or. ivisc(2).eq.9 .or. ivisc(1).eq.9 .or.
     .    ivisc(3).ge.11.or. ivisc(2).ge.11.or. ivisc(1).ge.11) then
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
      if (ivisc(3).eq.2 .or. ivisc(2).eq.2 .or. ivisc(1).eq.2) then
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
      else if (ivisc(3).eq.4 .or. ivisc(2).eq.4 .or. 
     .         ivisc(1).eq.4) then
c   Baldwin-Barth
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
     + sumn1,sumn2,negn1,negn2,xib2,volj0,volk0,voli0,nou,bou,nbuf,
     + ibufdim)
      else if (ivisc(3).eq.5 .or. ivisc(2).eq.5 .or. 
     .         ivisc(1).eq.5) then
c   Spalart
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
      nroom=nwork-(iwk18+inmx)
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
     + wk(iwk5),
     + wk(iwk6),wk(iwk7),wk(iwk8),wk(iwk9),wk(iwk10),wk(iwk11),
     + wk(iwk12),wk(iwk13),wk(iwk14),wk(iwk15),wk(iwk16),wk(iwk17),
     + wk(iwk18),ntime,tj0,tk0,ti0,nbl,qj0,qk0,qi0,blank,iover,
     + sumn1,sumn2,negn1,negn2,xib2,volj0,volk0,voli0,nou,bou,nbuf,
     + ibufdim)
      else if (ivisc(3).ge.6 .or. ivisc(2).ge.6 .or.
     .         ivisc(1).ge.6) then
c   two-equation turb model
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
      nroom=nwork-(iwk33+inmx)
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
     + vj0,vk0,vi0,blank,iover,sumn1,sumn2,negn1,negn2,ux,wk(iwk32),
     + xib2,wk(iwk33),cmuv,bcj,bck,bci,nbci0,nbcidim,
     + nbcj0,nbcjdim,nbck0,nbckdim,ibcinfo,jbcinfo,kbcinfo,
     + maxbl,maxseg,volj0,volk0,voli0,nou,bou,nbuf,ibufdim)
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
      end if
    7 format(1x,29hcomputing residuals for block,i4)
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
     .               volj0,nou,bou,nbuf,ibufdim,iadvance(nbl))
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
     .               iadvance(nbl))
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
     .                  voli0,nou,bou,nbuf,ibufdim,iadvance(nbl))
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
      return 
      end
