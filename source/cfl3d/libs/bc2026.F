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
      subroutine bc2026(jdim,kdim,idim,q,qj0,qk0,qi0,sj,sk,si,bcj,bck,
     .                  bci,xtbj,xtbk,xtbi,atbj,atbk,atbi,ista,iend,
     .                  jsta,jend,ksta,kend,nface,tursav,tj0,tk0,
     .                  ti0,smin,vist3d,vj0,vk0,vi0,mdim,ndim,bcdata,
     .                  filname,iuns,irelv,snj0,snk0,sni0,ntime,
     .                  snjm,snkm,snim,nou,bou,nbuf,ibufdim,myid,
     .                  nummem)
c
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Set "sweeping jet" blowing BC
c               with 9 pieces of auxiliary data and up to 2 turbulence quantities:
c     This BC is NOT treated as a "wall" bc; i.e., forces are not integrated on it,
c     and the min distance does not see it as a wall
c     (so ismincforce=2 should be invoked if restarting from a soln w/o bc2026)
c
c           1) vmag = jet velocity magnitude/a_ref
c           2) rfreq = reduced frequency of unsteady jet sweeping =
c                      thetajet=sideangj*sin(twopi*rfreq*time)
c           3) sideangj = max angle (deg) that the jet "sweeps" in each direction
c                           (+-sideangj)
c           4) sxa = direction number for surface jet prior to sweeping (x-direction)
c           5) sya = direction number for surface jet prior to sweeping (y-direction)
c           6) sza = direction number for surface jet prior to sweeping (z-direction)
c           7) sxb = x-direction number for vector normal to (sxa,sya,sza), in
c                       plane perpendicular to body surface and pointing downstream
c           8) syb = y-direction number for vector normal to (sxa,sya,sza), in
c                       plane perpendicular to body surface and pointing downstream
c           9) szb = z-direction number for vector normal to (sxa,sya,sza), in
c                       plane perpendicular to body surface and pointing downstream
c              NOTES:
c              1.  input direction numbers do not have to be normalized
c              2.  rfreq = freq*lref/a_inf where freq is frequency in Hertz,
c                  lref is reference length, and a_inf is speed of sound
c
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
      character*80 filname
c
      dimension nou(nbuf)
      dimension q(jdim,kdim,idim,5), qi0(jdim,kdim,5,4),
     .          qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
      dimension sk(jdim,kdim,idim-1,5),si(jdim,kdim,idim,5),
     .          sj(jdim,kdim,idim-1,5)
      dimension xtbj(kdim,idim-1,3,2),xtbk(jdim,idim-1,3,2),
     .          xtbi(jdim,kdim,3,2),atbj(kdim,idim-1,3,2),
     .          atbk(jdim,idim-1,3,2),atbi(jdim,kdim,3,2)
      dimension bcdata(mdim,ndim,2,12)
      dimension tursav(jdim,kdim,idim,nummem),tj0(kdim,idim-1,nummem,4),
     .          tk0(jdim,idim-1,nummem,4),ti0(jdim,kdim,nummem,4),
     .          smin(jdim-1,kdim-1,idim-1),
     .          vj0(kdim,idim-1,1,4),vk0(jdim,idim-1,1,4),
     .          vi0(jdim,kdim,1,4),vist3d(jdim,kdim,idim)
      dimension snj0(jdim-1,kdim-1,idim-1),snk0(jdim-1,kdim-1,idim-1),
     .          sni0(jdim-1,kdim-1,idim-1)
      dimension snjm(jdim-1,kdim-1,idim-1),snkm(jdim-1,kdim-1,idim-1),
     .          snim(jdim-1,kdim-1,idim-1)
      dimension a0(7),a1(5),a2(3)
c
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /fluid2/ pr,prt,cbar
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /maxiv/ ivmx
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .        tur10(7)
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /reyue/ reue,tinf,ivisc(3)
      common /sklton/ isklton
      common /wallfun/ iwf(3)
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /conversion/ radtodeg
c  
c
      twopi = 2.0*acos(-1.0)
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
      jend1 = jend-1
      kend1 = kend-1
      iend1 = iend-1
c
c
c            * * * * * * * * * * * * * * * * * * * * * *
c            * standard boundary condition bctype=2026 *
c            * * * * * * * * * * * * * * * * * * * * * *
c
c******************************************************************************
c      j=1 boundary        sweeping jet                             bctype 2026
c******************************************************************************
      if (nface.eq.3) then
c
c     check to see if turbulence data is input (itrflg1 = 1) or
c     if freestream values are to be used (itrflg1 = 0); the check
c     assumes if the first point has been set, all points have been
c
      ipp     = 1
      itrflg1 = 0
      if (real(bcdata(1,1,ipp,10)) .gt. -1.e10) itrflg1 = 1
c
      do 400 ip=1,2
      do 400 i=ista,iend1
      ii = i-ista+1
c
      do 300 k=ksta,kend1
      kk = k-ksta+1
      vmag          = bcdata(kk,ii,ip,1)
      rfreq         = bcdata(kk,ii,ip,2)
      sideangj      = bcdata(kk,ii,ip,3)/radtodeg
      sxa           = bcdata(kk,ii,ip,4)
      sya           = bcdata(kk,ii,ip,5)
      sza           = bcdata(kk,ii,ip,6)
      sxb           = bcdata(kk,ii,ip,7)
      syb           = bcdata(kk,ii,ip,8)
      szb           = bcdata(kk,ii,ip,9)
c
c     surface velocities
      uub = 0.
      vvb = 0.
      wwb = 0.
c
c     for dynamic mesh, set velocity at wall to grid velocity at wall
c     if irelv > 0; otherwise, set to zero
c
      if (iuns.gt.0 .and. irelv.gt.0) then
      uub = xtbj(k,i,1,1)
      vvb = xtbj(k,i,2,1)
      wwb = xtbj(k,i,3,1)
      end if
c
c     find angle of sweeping jet relative to frame pointing up in z
      thetajet=sideangj*sin(twopi*rfreq*time)
      sjetx1=0.
      sjety1=sin(thetajet)
      sjetz1=cos(thetajet)
c
c     normalize vectors
      snorma   = sqrt(sxa*sxa+sya*sya+sza*sza)
      sxa = sxa/snorma
      sya = sya/snorma
      sza = sza/snorma
      snormb   = sqrt(sxb*sxb+syb*syb+szb*szb)
      sxb = sxb/snormb
      syb = syb/snormb
      szb = szb/snormb
c
c     check to make sure vectors are normal to each other
      dotp=sxa*sxb+sya*syb+sza*szb
      if (abs(dotp) .gt. 1.e-8) then
c       error!
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),*)' stopping in bc2026 ',
     .  '- dot product not zero: ',dotp
        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
c
c     rotate vector to required position and find u,v,w
      sjetx=sxb*sjetx1+(sya*szb-sza*syb)*sjety1+sxa*sjetz1
      sjety=syb*sjetx1+(sza*sxb-sxa*szb)*sjety1+sya*sjetz1
      sjetz=szb*sjetx1+(sxa*syb-sya*sxb)*sjety1+sza*sjetz1
      uub = uub + sjetx*vmag
      vvb = vvb + sjety*vmag
      wwb = wwb + sjetz*vmag
c
      qj0(k,i,1,ip) = q(1,k,i,1)
      qj0(k,i,2,ip) = uub
      qj0(k,i,3,ip) = vvb
      qj0(k,i,4,ip) = wwb
      qj0(k,i,5,ip) = q(1,k,i,5)

      bcj(k,i,1)   = 0.0
  300 continue
  400 continue
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 191 i=ista,iend1
        do 191 k=ksta,kend1
          vj0(k,i,1,1) = vist3d(1,k,i)
          vj0(k,i,1,2) = vist3d(1,k,i)
  191   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 101 i=ista,iend1
        ii = i-ista+1
        do 101 k=ksta,kend1
          kk = k-ksta+1
          ip  = 1
          t11 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(kk,ii,ip,9+l)
          ip  = 2
          t12 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(kk,ii,ip,9+l)
          tj0(k,i,l,1) = t11
          tj0(k,i,l,2) = t12
  101   continue
        enddo
      end if
      end if
c
      end if
c
c******************************************************************************
c      j=jdim boundary     sweeping jet                             bctype 2026
c******************************************************************************
      if (nface.eq.4) then
c
c     check to see if turbulence data is input (itrflg1 = 1) or
c     if freestream values are to be used (itrflg1 = 0); the check
c     assumes if the first point has been set, all points have been
c
      ipp     = 1
      itrflg1 = 0
      if (real(bcdata(1,1,ipp,10)) .gt. -1.e10) itrflg1 = 1
c
      do 800 ip=1,2
      do 800 i=ista,iend1
      ii = i-ista+1
c
      do 700 k=ksta,kend1
      kk = k-ksta+1
      vmag          = bcdata(kk,ii,ip,1)
      rfreq         = bcdata(kk,ii,ip,2)
      sideangj      = bcdata(kk,ii,ip,3)/radtodeg
      sxa           = bcdata(kk,ii,ip,4)
      sya           = bcdata(kk,ii,ip,5)
      sza           = bcdata(kk,ii,ip,6)
      sxb           = bcdata(kk,ii,ip,7)
      syb           = bcdata(kk,ii,ip,8)
      szb           = bcdata(kk,ii,ip,9)
c
c     surface velocities
      uub = 0.
      vvb = 0.
      wwb = 0.
c
c     for dynamic mesh, set velocity at wall to grid velocity at wall
c     if irelv > 0; otherwise, set to zero
c
      if (iuns.gt.0 .and. irelv.gt.0) then
      uub = xtbj(k,i,1,2)
      vvb = xtbj(k,i,2,2)
      wwb = xtbj(k,i,3,2)
      end if
c
c     find angle of sweeping jet relative to frame pointing up in z
      thetajet=sideangj*sin(twopi*rfreq*time)
      sjetx1=0.
      sjety1=sin(thetajet)
      sjetz1=cos(thetajet)
c
c     normalize vectors
      snorma   = sqrt(sxa*sxa+sya*sya+sza*sza)
      sxa = sxa/snorma
      sya = sya/snorma
      sza = sza/snorma
      snormb   = sqrt(sxb*sxb+syb*syb+szb*szb)
      sxb = sxb/snormb
      syb = syb/snormb
      szb = szb/snormb
c
c     check to make sure vectors are normal to each other
      dotp=sxa*sxb+sya*syb+sza*szb
      if (abs(dotp) .gt. 1.e-8) then
c       error!
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),*)' stopping in bc2026 ',
     .  '- dot product not zero: ',dotp
        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
c
c     rotate vector to required position and find u,v,w
      sjetx=sxb*sjetx1+(sya*szb-sza*syb)*sjety1+sxa*sjetz1
      sjety=syb*sjetx1+(sza*sxb-sxa*szb)*sjety1+sya*sjetz1
      sjetz=szb*sjetx1+(sxa*syb-sya*sxb)*sjety1+sza*sjetz1
      uub = uub + sjetx*vmag
      vvb = vvb + sjety*vmag
      wwb = wwb + sjetz*vmag
c
      qj0(k,i,1,ip+2) = q(jdim1,k,i,1)
      qj0(k,i,2,ip+2) = uub
      qj0(k,i,3,ip+2) = vvb
      qj0(k,i,4,ip+2) = wwb
      qj0(k,i,5,ip+2) = q(jdim1,k,i,5)

      bcj(k,i,2)   = 0.0
  700 continue
  800 continue
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 291 i=ista,iend1
        do 291 k=ksta,kend1
          vj0(k,i,1,3) = vist3d(jdim1,k,i)
          vj0(k,i,1,4) = vist3d(jdim1,k,i)
  291   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 201 i=ista,iend1
        ii = i-ista+1
        do 201 k=ksta,kend1
          kk = k-ksta+1
          ip  = 1
          t13 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(kk,ii,ip,9+l)
          ip  = 2
          t14 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(kk,ii,ip,9+l)
          tj0(k,i,l,3) = t13
          tj0(k,i,l,4) = t14
  201   continue
        enddo
      end if
      end if
c
      end if
c
c******************************************************************************
c      k=1 boundary        sweeping jet                             bctype 2026
c******************************************************************************
      if (nface.eq.5) then
c
c     check to see if turbulence data is input (itrflg1 = 1) or
c     if freestream values are to be used (itrflg1 = 0); the check
c     assumes if the first point has been set, all points have been
c
      ipp     = 1
      itrflg1 = 0
      if (real(bcdata(1,1,ipp,10)) .gt. -1.e10) itrflg1 = 1
c
      do 1200 ip=1,2
      do 1200 i=ista,iend1
      ii = i-ista+1
c
      do 1100 j=jsta,jend1
      jj = j-jsta+1
      vmag          = bcdata(jj,ii,ip,1) 
      rfreq         = bcdata(jj,ii,ip,2)
      sideangj      = bcdata(jj,ii,ip,3)/radtodeg
      sxa           = bcdata(jj,ii,ip,4)
      sya           = bcdata(jj,ii,ip,5)
      sza           = bcdata(jj,ii,ip,6)
      sxb           = bcdata(jj,ii,ip,7)
      syb           = bcdata(jj,ii,ip,8)
      szb           = bcdata(jj,ii,ip,9)
c
c     surface velocities
      uub = 0.
      vvb = 0.
      wwb = 0.
c
c     for dynamic mesh, set velocity at wall to grid velocity at wall
c     if irelv > 0; otherwise, set to zero
      if (iuns.gt.0 .and. irelv.gt.0) then
      uub = xtbk(j,i,1,1)
      vvb = xtbk(j,i,2,1)
      wwb = xtbk(j,i,3,1)
      end if
c
c     find angle of sweeping jet relative to frame pointing up in z
      thetajet=sideangj*sin(twopi*rfreq*time)
      sjetx1=0.
      sjety1=sin(thetajet)
      sjetz1=cos(thetajet)
c
c     normalize vectors
      snorma   = sqrt(sxa*sxa+sya*sya+sza*sza)
      sxa = sxa/snorma
      sya = sya/snorma
      sza = sza/snorma
      snormb   = sqrt(sxb*sxb+syb*syb+szb*szb)
      sxb = sxb/snormb
      syb = syb/snormb
      szb = szb/snormb
c
c     check to make sure vectors are normal to each other
      dotp=sxa*sxb+sya*syb+sza*szb
      if (abs(dotp) .gt. 1.e-8) then
c       error!
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),*)' stopping in bc2026 ',
     .  '- dot product not zero: ',dotp
        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
c
c     rotate vector to required position and find u,v,w
      sjetx=sxb*sjetx1+(sya*szb-sza*syb)*sjety1+sxa*sjetz1
      sjety=syb*sjetx1+(sza*sxb-sxa*szb)*sjety1+sya*sjetz1
      sjetz=szb*sjetx1+(sxa*syb-sya*sxb)*sjety1+sza*sjetz1
      uub = uub + sjetx*vmag
      vvb = vvb + sjety*vmag
      wwb = wwb + sjetz*vmag
c
      qk0(j,i,1,ip) = q(j,1,i,1)
      qk0(j,i,2,ip) = uub
      qk0(j,i,3,ip) = vvb
      qk0(j,i,4,ip) = wwb
      qk0(j,i,5,ip) = q(j,1,i,5)

      bck(j,i,1)   = 0.0
 1100 continue
 1200 continue
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 391 i=ista,iend1
        do 391 j=jsta,jend1
          vk0(j,i,1,1) = vist3d(j,1,i)
          vk0(j,i,1,2) = vist3d(j,1,i)
  391   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 301 i=ista,iend1
        ii = i-ista+1
        do 301 j=jsta,jend1
          jj = j-jsta+1
          ip  = 1
          t11 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(jj,ii,ip,9+l)
          ip  = 2
          t12 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(jj,ii,ip,9+l)
          tk0(j,i,l,1) = t11
          tk0(j,i,l,2) = t12
  301   continue
        enddo
      end if
      end if
c
      end if
c
c******************************************************************************
c      k=kdim boundary     sweeping jet                             bctype 2026
c******************************************************************************
      if (nface.eq.6) then
c
c     check to see if turbulence data is input (itrflg1 = 1) or
c     if freestream values are to be used (itrflg1 = 0); the check
c     assumes if the first point has been set, all points have been
c
      ipp     = 1
      itrflg1 = 0
      if (real(bcdata(1,1,ipp,10)) .gt. -1.e10) itrflg1 = 1
c
      do 1400 ip=1,2
      do 1400 i=ista,iend1
      ii = i-ista+1
c
      do 1300 j=jsta,jend1
      jj = j-jsta+1
      vmag          = bcdata(jj,ii,ip,1)
      rfreq         = bcdata(jj,ii,ip,2)
      sideangj      = bcdata(jj,ii,ip,3)/radtodeg
      sxa           = bcdata(jj,ii,ip,4)
      sya           = bcdata(jj,ii,ip,5)
      sza           = bcdata(jj,ii,ip,6)
      sxb           = bcdata(jj,ii,ip,7)
      syb           = bcdata(jj,ii,ip,8)
      szb           = bcdata(jj,ii,ip,9)
c
c     surface velocities
c
      uub = 0.
      vvb = 0.
      wwb = 0.
c
c     for dynamic mesh, set velocity at wall to grid velocity at wall
c     if irelv > 0; otherwise, set to zero
      if (iuns.gt.0 .and. irelv.gt.0) then
      uub = xtbk(j,i,1,2)
      vvb = xtbk(j,i,2,2)
      wwb = xtbk(j,i,3,2)
      end if
c
c     find angle of sweeping jet relative to frame pointing up in z
      thetajet=sideangj*sin(twopi*rfreq*time)
      sjetx1=0.
      sjety1=sin(thetajet)
      sjetz1=cos(thetajet)
c
c     normalize vectors
      snorma   = sqrt(sxa*sxa+sya*sya+sza*sza)
      sxa = sxa/snorma
      sya = sya/snorma
      sza = sza/snorma
      snormb   = sqrt(sxb*sxb+syb*syb+szb*szb)
      sxb = sxb/snormb
      syb = syb/snormb
      szb = szb/snormb
c
c     check to make sure vectors are normal to each other
      dotp=sxa*sxb+sya*syb+sza*szb
      if (abs(dotp) .gt. 1.e-8) then
c       error!
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),*)' stopping in bc2026 ',
     .  '- dot product not zero: ',dotp
        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
c
c     rotate vector to required position and find u,v,w
      sjetx=sxb*sjetx1+(sya*szb-sza*syb)*sjety1+sxa*sjetz1
      sjety=syb*sjetx1+(sza*sxb-sxa*szb)*sjety1+sya*sjetz1
      sjetz=szb*sjetx1+(sxa*syb-sya*sxb)*sjety1+sza*sjetz1
      uub = uub + sjetx*vmag
      vvb = vvb + sjety*vmag
      wwb = wwb + sjetz*vmag
c
      qk0(j,i,1,ip+2) = q(j,kdim1,i,1)
      qk0(j,i,2,ip+2) = uub
      qk0(j,i,3,ip+2) = vvb
      qk0(j,i,4,ip+2) = wwb
      qk0(j,i,5,ip+2) = q(j,kdim1,i,5)

      bck(j,i,2)   = 0.0
 1300 continue
 1400 continue
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 491 i=ista,iend1
        do 491 j=jsta,jend1
          vk0(j,i,1,3) = vist3d(j,kdim1,i)
          vk0(j,i,1,4) = vist3d(j,kdim1,i)
  491   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 401 i=ista,iend1
        ii = i-ista+1
        do 401 j=jsta,jend1
          jj = j-jsta+1
          ip  = 1
          t13 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(jj,ii,ip,9+l)
          ip  = 2
          t14 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(jj,ii,ip,9+l)
          tk0(j,i,l,3) = t13
          tk0(j,i,l,4) = t14
  401   continue
        enddo
      end if
      end if
c
      end if
c
c******************************************************************************
c      i=1 boundary        sweeping jet                             bctype 2026
c******************************************************************************
      if (nface.eq.1) then
c
c     check to see if turbulence data is input (itrflg1 = 1) or
c     if freestream values are to be used (itrflg1 = 0); the check
c     assumes if the first point has been set, all points have been
c
      ipp     = 1
      itrflg1 = 0
      if (real(bcdata(1,1,ipp,10)) .gt. -1.e10) itrflg1 = 1
c 
      do 2000 ip=1,2
      do 2000 k=ksta,kend1
      kk = k-ksta+1
c
      do 1900 j=jsta,jend1
      jj = j-jsta+1
      vmag          = bcdata(jj,kk,ip,1)
      rfreq         = bcdata(jj,kk,ip,2)
      sideangj      = bcdata(jj,kk,ip,3)/radtodeg
      sxa           = bcdata(jj,kk,ip,4)
      sya           = bcdata(jj,kk,ip,5)
      sza           = bcdata(jj,kk,ip,6)
      sxb           = bcdata(jj,kk,ip,7)
      syb           = bcdata(jj,kk,ip,8)
      szb           = bcdata(jj,kk,ip,9)
c
c     surface velocities
      uub = 0.
      vvb = 0.
      wwb = 0.
c
c     for dynamic mesh, set velocity at wall to grid velocity at wall
c     if irelv > 0; otherwise, set to zero
c
      if (iuns.gt.0 .and. irelv.gt.0) then
      uub = xtbi(j,k,1,1)
      vvb = xtbi(j,k,2,1)
      wwb = xtbi(j,k,3,1)
      end if
c
c     find angle of sweeping jet relative to frame pointing up in z
      thetajet=sideangj*sin(twopi*rfreq*time)
      sjetx1=0.
      sjety1=sin(thetajet)
      sjetz1=cos(thetajet)
c
c     normalize vectors
      snorma   = sqrt(sxa*sxa+sya*sya+sza*sza)
      sxa = sxa/snorma
      sya = sya/snorma
      sza = sza/snorma
      snormb   = sqrt(sxb*sxb+syb*syb+szb*szb)
      sxb = sxb/snormb
      syb = syb/snormb
      szb = szb/snormb
c
c     check to make sure vectors are normal to each other
      dotp=sxa*sxb+sya*syb+sza*szb
      if (abs(dotp) .gt. 1.e-8) then
c       error!
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),*)' stopping in bc2026 ',
     .  '- dot product not zero: ',dotp
        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
c
c     rotate vector to required position and find u,v,w
      sjetx=sxb*sjetx1+(sya*szb-sza*syb)*sjety1+sxa*sjetz1
      sjety=syb*sjetx1+(sza*sxb-sxa*szb)*sjety1+sya*sjetz1
      sjetz=szb*sjetx1+(sxa*syb-sya*sxb)*sjety1+sza*sjetz1
      uub = uub + sjetx*vmag
      vvb = vvb + sjety*vmag
      wwb = wwb + sjetz*vmag
c
      qi0(j,k,1,ip) = q(j,k,1,1)
      qi0(j,k,2,ip) = uub
      qi0(j,k,3,ip) = vvb
      qi0(j,k,4,ip) = wwb
      qi0(j,k,5,ip) = q(j,k,1,5)

      bci(j,k,1)   = 0.0
 1900 continue
 2000 continue
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 591 k=ksta,kend1
        do 591 j=jsta,jend1
          vi0(j,k,1,1) = vist3d(j,k,1)
          vi0(j,k,1,2) = vist3d(j,k,1)
  591   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 501 k=ksta,kend1
        kk = k-ksta+1
        do 501 j=jsta,jend1
          jj = j-jsta+1
          ip  = 1
          t11 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(jj,kk,ip,9+l)
          ip  = 2
          t12 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(jj,kk,ip,9+l)
          ti0(j,k,l,1) = t11
          ti0(j,k,l,2) = t12
  501   continue
        enddo
      end if
      end if
c
      end if
c
c******************************************************************************
c      i=idim boundary     sweeping jet                             bctype 2026
c******************************************************************************
      if (nface.eq.2) then
c
c     check to see if turbulence data is input (itrflg1 = 1) or
c     if freestream values are to be used (itrflg1 = 0); the check
c     assumes if the first point has been set, all points have been
c
      ipp     = 1
      itrflg1 = 0
      if (real(bcdata(1,1,ipp,10)) .gt. -1.e10) itrflg1 = 1
c
      do 2400 ip=1,2
      do 2400 k=ksta,kend1
      kk = k-ksta+1
c
      do 2300 j=jsta,jend1
      jj = j-jsta+1
      vmag          = bcdata(jj,kk,ip,1)
      rfreq         = bcdata(jj,kk,ip,2)
      sideangj      = bcdata(jj,kk,ip,3)/radtodeg
      sxa           = bcdata(jj,kk,ip,4)
      sya           = bcdata(jj,kk,ip,5)
      sza           = bcdata(jj,kk,ip,6)
      sxb           = bcdata(jj,kk,ip,7)
      syb           = bcdata(jj,kk,ip,8)
      szb           = bcdata(jj,kk,ip,9)
c
c     surface velocities
      uub = 0.
      vvb = 0.
      wwb = 0.
c
c     for dynamic mesh, set velocity at wall to grid velocity at wall
c     if irelv > 0; otherwise, set to zero
      if (iuns.gt.0 .and. irelv.gt.0) then
      uub = xtbi(j,k,1,2)
      vvb = xtbi(j,k,2,2)
      wwb = xtbi(j,k,3,2)
      end if
c
c     find angle of sweeping jet relative to frame pointing up in z
      thetajet=sideangj*sin(twopi*rfreq*time)
      sjetx1=0.
      sjety1=sin(thetajet)
      sjetz1=cos(thetajet)
c
c     normalize vectors
      snorma   = sqrt(sxa*sxa+sya*sya+sza*sza)
      sxa = sxa/snorma
      sya = sya/snorma
      sza = sza/snorma
      snormb   = sqrt(sxb*sxb+syb*syb+szb*szb)
      sxb = sxb/snormb
      syb = syb/snormb
      szb = szb/snormb
c
c     check to make sure vectors are normal to each other
      dotp=sxa*sxb+sya*syb+sza*szb
      if (abs(dotp) .gt. 1.e-8) then
c       error!
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),*)' stopping in bc2026 ',
     .  '- dot product not zero: ',dotp
        call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
c
c     rotate vector to required position and find u,v,w
      sjetx=sxb*sjetx1+(sya*szb-sza*syb)*sjety1+sxa*sjetz1
      sjety=syb*sjetx1+(sza*sxb-sxa*szb)*sjety1+sya*sjetz1
      sjetz=szb*sjetx1+(sxa*syb-sya*sxb)*sjety1+sza*sjetz1
      uub = uub + sjetx*vmag
      vvb = vvb + sjety*vmag
      wwb = wwb + sjetz*vmag
c
      qi0(j,k,1,ip+2) = q(j,k,idim1,1)
      qi0(j,k,2,ip+2) = uub
      qi0(j,k,3,ip+2) = vvb
      qi0(j,k,4,ip+2) = wwb
      qi0(j,k,5,ip+2) = q(j,k,idim1,5)

      bci(j,k,2)   = 0.0
 2300 continue
 2400 continue
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 691 k=ksta,kend1
        do 691 j=jsta,jend1
          vi0(j,k,1,3) = vist3d(j,k,idim1)
          vi0(j,k,1,4) = vist3d(j,k,idim1)
  691   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 601 k=ksta,kend1
        kk = k-ksta+1
        do 601 j=jsta,jend1
          jj = j-jsta+1
          ip  = 1
          t13 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(jj,kk,ip,9+l)
          ip  = 2
          t14 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(jj,kk,ip,9+l)
          ti0(j,k,l,3) = t13
          ti0(j,k,l,4) = t14
  601   continue
        enddo
      end if
      end if
      end if
c
      return
      end
