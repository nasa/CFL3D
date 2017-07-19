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
      subroutine bc2008(jdim,kdim,idim,q,qj0,qk0,qi0,sj,sk,si,bcj,bck,
     .                  bci,xtbj,xtbk,xtbi,atbj,atbk,atbi,ista,iend,
     .                  jsta,jend,ksta,kend,nface,tursav,tj0,tk0,
     .                  ti0,vist3d,vj0,vk0,vi0,mdim,ndim,bcdata,
     .                  filname,iuns,nou,bou,nbuf,ibufdim,myid,iflag,
     .                  nummem)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Set rho/rho_inf, u/a_inf, v/a_inf, w/a_inf, and up to
c               two turbulence quantities; the pressure is extrapolated 
c               from the interior
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
      dimension sk(jdim,kdim,idim-1,5),si(jdim,kdim,idim,5),
     .          sj(jdim,kdim,idim-1,5)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
      dimension xtbj(kdim,idim-1,3,2),xtbk(jdim,idim-1,3,2),
     .          xtbi(jdim,kdim,3,2),atbj(kdim,idim-1,3,2),
     .          atbk(jdim,idim-1,3,2),atbi(jdim,kdim,3,2)
      dimension bcdata(mdim,ndim,2,12)
      dimension tursav(jdim,kdim,idim,nummem),tj0(kdim,idim-1,nummem,4),
     .          tk0(jdim,idim-1,nummem,4),ti0(jdim,kdim,nummem,4),
     .          vj0(kdim,idim-1,1,4),vk0(jdim,idim-1,1,4),
     .          vi0(jdim,kdim,1,4),vist3d(jdim,kdim,idim)
      integer   stats
c
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /maxiv/ ivmx
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .        tur10(7)
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /reyue/ reue,tinf,ivisc(3)
      common /sklton/ isklton
      common /igrdtyp/ ip3dgrd,ialph
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
c
      real, allocatable :: harvest(:)
c
      twopi = 8.0*atan(1.0)
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
      jend1 = jend-1
      kend1 = kend-1
      iend1 = iend-1
c
      if (iflag .eq. 2038) then
c       generate array of random numbers
        if (nface .eq. 1 .or. nface .eq. 2) jsize=
     .     (jend1-jsta+1)*(kend1-ksta+1)
        if (nface .eq. 3 .or. nface .eq. 4) jsize=
     .     (iend1-ista+1)*(kend1-ksta+1)
        if (nface .eq. 5 .or. nface .eq. 6) jsize=
     .     (iend1-ista+1)*(jend1-jsta+1)
        allocate( harvest(jsize),stat=stats )
        call umalloc(jsize,0,'harvest',memuse,stats)
        call random_seed
        call random_number(harvest)
      end if
c
c     this bc makes use of two planes of data (for ndata>0 in the input
c     deck, the two planes are identical; for ndata<0 in the input deck,
c     the two planes are set by the data file and may differ to reflect
c     two distinct planes of ghost cells)
c
c            * * * * * * * * * * * * * * * * * * * * * *
c            * standard boundary condition bctype=2008 *
c            * * * * * * * * * * * * * * * * * * * * * *
c
c******************************************************************************
c      j=1 boundary      set velocity and turb quantities           bctype 2008
c******************************************************************************
c
      if (nface.eq.3) then
c
c     check to see if turbulence data is input (itrflg1 = 1) or
c     if freestream values are to be used (itrflg1 = 0); the check
c     assumes if the first point has been set, all points have been
c
      ipp     = 1
      itrflg1 = 0
      if (real(bcdata(1,1,ipp,5)) .gt. -1.e10) itrflg1 = 1
c
      if (ialph.eq.0) then
         do 100 ip=1,2
         irnd_count=0
         do 100 i=ista,iend1
         ii = i-ista+1
         do 100 k=ksta,kend1
         kk = k-ksta+1
         irnd_count=irnd_count+1
         if (iflag .eq. 2008) then
         qj0(k,i,1,ip) = bcdata(kk,ii,ip,1)
         qj0(k,i,2,ip) = bcdata(kk,ii,ip,2)
         qj0(k,i,3,ip) = bcdata(kk,ii,ip,3)
         qj0(k,i,4,ip) = bcdata(kk,ii,ip,4)
         qj0(k,i,5,ip) = q(1,k,i,5)
         else if (iflag .eq. 2018) then
         qj0(k,i,5,ip) = q(1,k,i,5)
         qj0(k,i,1,ip) = gamma*qj0(k,i,5,ip)/bcdata(kk,ii,ip,1)
         qj0(k,i,2,ip) = bcdata(kk,ii,ip,2)/qj0(k,i,1,ip)
         qj0(k,i,3,ip) = bcdata(kk,ii,ip,3)/qj0(k,i,1,ip)
         qj0(k,i,4,ip) = bcdata(kk,ii,ip,4)/qj0(k,i,1,ip)
         else if (iflag .eq. 2028) then
         qj0(k,i,1,ip) = q(1,k,i,1)
         qj0(k,i,2,ip) = bcdata(kk,ii,ip,2)/qj0(k,i,1,ip)*
     +     cos(twopi*bcdata(kk,ii,ip,1)*time)
         qj0(k,i,3,ip) = bcdata(kk,ii,ip,3)/qj0(k,i,1,ip)*
     +     cos(twopi*bcdata(kk,ii,ip,1)*time)
         qj0(k,i,4,ip) = bcdata(kk,ii,ip,4)/qj0(k,i,1,ip)*
     +     cos(twopi*bcdata(kk,ii,ip,1)*time)
         qj0(k,i,5,ip) = q(1,k,i,5)
         else if (iflag .eq. 2038) then
         rnd_nbr=harvest(irnd_count)-0.5
         if (real(bcdata(kk,ii,ip,2)) .gt. 0.85*real(xmach)) rnd_nbr=0.
         if (real(bcdata(kk,ii,ip,2)) .gt. 0.75*real(xmach) .and.
     +       real(bcdata(kk,ii,ip,2)) .le. 0.85*real(xmach)) then
           rnd_nbr = (-10./xmach*bcdata(kk,ii,ip,2)+8.5)*rnd_nbr
         end if
         qj0(k,i,1,ip) = bcdata(kk,ii,ip,1)*(1.0+rnd_nbr*0.10)
         qj0(k,i,2,ip) = bcdata(kk,ii,ip,2)*(1.0+rnd_nbr*0.10)
         qj0(k,i,3,ip) = bcdata(kk,ii,ip,3)*(1.0+rnd_nbr*0.10)
         qj0(k,i,4,ip) = bcdata(kk,ii,ip,4)*(1.0+rnd_nbr*0.10)
         qj0(k,i,5,ip) = q(1,k,i,5)
         end if
         bcj(k,i,1) = 0.0
  100    continue
      else
         do 1000 ip=1,2
         irnd_count=0
         do 1000 i=ista,iend1
         ii = i-ista+1
         do 1000 k=ksta,kend1
         kk = k-ksta+1
         irnd_count=irnd_count+1
         if (iflag .eq. 2008) then
         qj0(k,i,1,ip) =  bcdata(kk,ii,ip,1)
         qj0(k,i,2,ip) =  bcdata(kk,ii,ip,2)
         qj0(k,i,3,ip) = -bcdata(kk,ii,ip,4)
         qj0(k,i,4,ip) =  bcdata(kk,ii,ip,3)
         qj0(k,i,5,ip) =  q(1,k,i,5)
         else if (iflag .eq. 2018) then
         qj0(k,i,5,ip) =  q(1,k,i,5)
         qj0(k,i,1,ip) =  gamma*qj0(k,i,5,ip)/bcdata(kk,ii,ip,1)
         qj0(k,i,2,ip) =  bcdata(kk,ii,ip,2)/qj0(k,i,1,ip)
         qj0(k,i,3,ip) = -bcdata(kk,ii,ip,4)/qj0(k,i,1,ip)
         qj0(k,i,4,ip) =  bcdata(kk,ii,ip,3)/qj0(k,i,1,ip)
         else if (iflag .eq. 2028) then
         qj0(k,i,1,ip) =  q(1,k,i,1)
         qj0(k,i,2,ip) =  bcdata(kk,ii,ip,2)/qj0(k,i,1,ip)*
     +     cos(twopi*bcdata(kk,ii,ip,1)*time)
         qj0(k,i,3,ip) = -bcdata(kk,ii,ip,4)/qj0(k,i,1,ip)*
     +     cos(twopi*bcdata(kk,ii,ip,1)*time)
         qj0(k,i,4,ip) =  bcdata(kk,ii,ip,3)/qj0(k,i,1,ip)*
     +     cos(twopi*bcdata(kk,ii,ip,1)*time)
         qj0(k,i,5,ip) =  q(1,k,i,5)
         else if (iflag .eq. 2038) then
         rnd_nbr=harvest(irnd_count)-0.5
         if (real(bcdata(kk,ii,ip,2)) .gt. 0.85*real(xmach)) rnd_nbr=0.
         if (real(bcdata(kk,ii,ip,2)) .gt. 0.75*real(xmach) .and.
     +       real(bcdata(kk,ii,ip,2)) .le. 0.85*real(xmach)) then
           rnd_nbr = (-10./xmach*bcdata(kk,ii,ip,2)+8.5)*rnd_nbr
         end if
         qj0(k,i,1,ip) = bcdata(kk,ii,ip,1)*(1.0+rnd_nbr*0.10)
         qj0(k,i,2,ip) = bcdata(kk,ii,ip,2)*(1.0+rnd_nbr*0.10)
         qj0(k,i,3,ip) =-bcdata(kk,ii,ip,4)*(1.0+rnd_nbr*0.10)
         qj0(k,i,4,ip) = bcdata(kk,ii,ip,3)*(1.0+rnd_nbr*0.10)
         qj0(k,i,5,ip) = q(1,k,i,5)
         end if
         bcj(k,i,1) = 0.0
 1000    continue
      end if
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
          t11 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(kk,ii,ip,4+l)
          ip  = 2
          t12 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(kk,ii,ip,4+l)
          tj0(k,i,l,1) = t11
          tj0(k,i,l,2) = t12
  101   continue
        enddo
      end if
      end if
      end if
c
c******************************************************************************
c      j=jdim boundary   set velocity and turb quantities           bctype 2008
c******************************************************************************
c 
      if (nface.eq.4) then
c
c     check to see if turbulence data is input (itrflg1 = 1) or
c     if freestream values are to be used (itrflg1 = 0); the check
c     assumes if the first point has been set, all points have been
c
      ipp     = 1
      itrflg1 = 0
      if (real(bcdata(1,1,ipp,5)) .gt. -1.e10) itrflg1 = 1
c
      if (ialph.eq.0) then
         do 200 ip=1,2
         irnd_count=0
         do 200 i=ista,iend1
         ii = i-ista+1
         do 200 k=ksta,kend1
         kk = k-ksta+1
         irnd_count=irnd_count+1
         if (iflag .eq. 2008) then
         qj0(k,i,1,ip+2) = bcdata(kk,ii,ip,1)
         qj0(k,i,2,ip+2) = bcdata(kk,ii,ip,2)
         qj0(k,i,3,ip+2) = bcdata(kk,ii,ip,3)
         qj0(k,i,4,ip+2) = bcdata(kk,ii,ip,4)
         qj0(k,i,5,ip+2) = q(jdim1,k,i,5)
         else if (iflag .eq. 2018) then
         qj0(k,i,5,ip+2) = q(jdim1,k,i,5)
         qj0(k,i,1,ip+2) = gamma*qj0(k,i,5,ip+2)/bcdata(kk,ii,ip,1)
         qj0(k,i,2,ip+2) = bcdata(kk,ii,ip,2)/qj0(k,i,1,ip+2)
         qj0(k,i,3,ip+2) = bcdata(kk,ii,ip,3)/qj0(k,i,1,ip+2)
         qj0(k,i,4,ip+2) = bcdata(kk,ii,ip,4)/qj0(k,i,1,ip+2)
         else if (iflag .eq. 2028) then
         qj0(k,i,1,ip+2) = q(jdim1,k,i,1)
         qj0(k,i,2,ip+2) = bcdata(kk,ii,ip,2)/qj0(k,i,1,ip+2)*
     +     cos(twopi*bcdata(kk,ii,ip,1)*time)
         qj0(k,i,3,ip+2) = bcdata(kk,ii,ip,3)/qj0(k,i,1,ip+2)*
     +     cos(twopi*bcdata(kk,ii,ip,1)*time)
         qj0(k,i,4,ip+2) = bcdata(kk,ii,ip,4)/qj0(k,i,1,ip+2)*
     +     cos(twopi*bcdata(kk,ii,ip,1)*time)
         qj0(k,i,5,ip+2) = q(jdim1,k,i,5)
         else if (iflag .eq. 2038) then
         rnd_nbr=harvest(irnd_count)-0.5
         if (real(bcdata(kk,ii,ip,2)) .gt. 0.85*real(xmach)) rnd_nbr=0.
         if (real(bcdata(kk,ii,ip,2)) .gt. 0.75*real(xmach) .and.
     +       real(bcdata(kk,ii,ip,2)) .le. 0.85*real(xmach)) then
           rnd_nbr = (-10./xmach*bcdata(kk,ii,ip,2)+8.5)*rnd_nbr
         end if
         qj0(k,i,1,ip+2) = bcdata(kk,ii,ip,1)*(1.0+rnd_nbr*0.10)
         qj0(k,i,2,ip+2) = bcdata(kk,ii,ip,2)*(1.0+rnd_nbr*0.10)
         qj0(k,i,3,ip+2) = bcdata(kk,ii,ip,3)*(1.0+rnd_nbr*0.10)
         qj0(k,i,4,ip+2) = bcdata(kk,ii,ip,4)*(1.0+rnd_nbr*0.10)
         qj0(k,i,5,ip+2) = q(jdim1,k,i,5)
         end if
         bcj(k,i,2) = 0.0
  200    continue
      else
         do 2000 ip=1,2
         irnd_count=0
         do 2000 i=ista,iend1
         ii = i-ista+1
         do 2000 k=ksta,kend1
         kk = k-ksta+1
         irnd_count=irnd_count+1
         if (iflag .eq. 2008) then
         qj0(k,i,1,ip+2) =  bcdata(kk,ii,ip,1)
         qj0(k,i,2,ip+2) =  bcdata(kk,ii,ip,2)
         qj0(k,i,3,ip+2) = -bcdata(kk,ii,ip,4)
         qj0(k,i,4,ip+2) =  bcdata(kk,ii,ip,3)
         qj0(k,i,5,ip+2) =  q(jdim1,k,i,5)
         else if (iflag .eq. 2018) then
         qj0(k,i,5,ip+2) =  q(jdim1,k,i,5)
         qj0(k,i,1,ip+2) =  gamma*qj0(k,i,5,ip+2)/bcdata(kk,ii,ip,1)
         qj0(k,i,2,ip+2) =  bcdata(kk,ii,ip,2)/qj0(k,i,1,ip+2)
         qj0(k,i,3,ip+2) = -bcdata(kk,ii,ip,4)/qj0(k,i,1,ip+2)
         qj0(k,i,4,ip+2) =  bcdata(kk,ii,ip,3)/qj0(k,i,1,ip+2)
         else if (iflag .eq. 2028) then
         qj0(k,i,1,ip+2) =  q(jdim1,k,i,1)
         qj0(k,i,2,ip+2) =  bcdata(kk,ii,ip,2)/qj0(k,i,1,ip+2)*
     +     cos(twopi*bcdata(kk,ii,ip,1)*time)
         qj0(k,i,3,ip+2) = -bcdata(kk,ii,ip,4)/qj0(k,i,1,ip+2)*
     +     cos(twopi*bcdata(kk,ii,ip,1)*time)
         qj0(k,i,4,ip+2) =  bcdata(kk,ii,ip,3)/qj0(k,i,1,ip+2)*
     +     cos(twopi*bcdata(kk,ii,ip,1)*time)
         qj0(k,i,5,ip+2) =  q(jdim1,k,i,5)
         else if (iflag .eq. 2038) then
         rnd_nbr=harvest(irnd_count)-0.5
         if (real(bcdata(kk,ii,ip,2)) .gt. 0.85*real(xmach)) rnd_nbr=0.
         if (real(bcdata(kk,ii,ip,2)) .gt. 0.75*real(xmach) .and.
     +       real(bcdata(kk,ii,ip,2)) .le. 0.85*real(xmach)) then
           rnd_nbr = (-10./xmach*bcdata(kk,ii,ip,2)+8.5)*rnd_nbr
         end if
         qj0(k,i,1,ip+2) = bcdata(kk,ii,ip,1)*(1.0+rnd_nbr*0.10)
         qj0(k,i,2,ip+2) = bcdata(kk,ii,ip,2)*(1.0+rnd_nbr*0.10)
         qj0(k,i,3,ip+2) =-bcdata(kk,ii,ip,4)*(1.0+rnd_nbr*0.10)
         qj0(k,i,4,ip+2) = bcdata(kk,ii,ip,3)*(1.0+rnd_nbr*0.10)
         qj0(k,i,5,ip+2) = q(jdim1,k,i,5)
         end if
         bcj(k,i,2) = 0.0
 2000    continue
      end if
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
          t13 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(kk,ii,ip,4+l)
          ip  = 2
          t14 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(kk,ii,ip,4+l)
          tj0(k,i,l,3) = t13
          tj0(k,i,l,4) = t14
  201   continue
        enddo
      end if
      end if
      end if
c
c******************************************************************************
c      k=1 boundary      set velocity and turb quantities           bctype 2008
c******************************************************************************
c
      if (nface.eq.5) then
c
c     check to see if turbulence data is input (itrflg1 = 1) or
c     if freestream values are to be used (itrflg1 = 0); the check
c     assumes if the first point has been set, all points have been
c
      ipp     = 1
      itrflg1 = 0
      if (real(bcdata(1,1,ipp,5)) .gt. -1.e10) itrflg1 = 1
c
      if (ialph.eq.0) then
         do 300 ip=1,2
         irnd_count=0
         do 300 i=ista,iend1
         ii = i-ista+1
         do 300 j=jsta,jend1
         jj = j-jsta+1
         irnd_count=irnd_count+1
         if (iflag .eq. 2008) then
         qk0(j,i,1,ip) = bcdata(jj,ii,ip,1)
         qk0(j,i,2,ip) = bcdata(jj,ii,ip,2)
         qk0(j,i,3,ip) = bcdata(jj,ii,ip,3)
         qk0(j,i,4,ip) = bcdata(jj,ii,ip,4)
         qk0(j,i,5,ip) = q(j,1,i,5)
         else if (iflag .eq. 2018) then
         qk0(j,i,5,ip) = q(j,1,i,5)
         qk0(j,i,1,ip) = gamma*qk0(j,i,5,ip)/bcdata(jj,ii,ip,1)
         qk0(j,i,2,ip) = bcdata(jj,ii,ip,2)/qk0(j,i,1,ip)
         qk0(j,i,3,ip) = bcdata(jj,ii,ip,3)/qk0(j,i,1,ip)
         qk0(j,i,4,ip) = bcdata(jj,ii,ip,4)/qk0(j,i,1,ip)
         else if (iflag .eq. 2028) then
         qk0(j,i,1,ip) = q(j,1,i,1)
         qk0(j,i,2,ip) = bcdata(jj,ii,ip,2)/qk0(j,i,1,ip)*
     +     cos(twopi*bcdata(jj,ii,ip,1)*time)
         qk0(j,i,3,ip) = bcdata(jj,ii,ip,3)/qk0(j,i,1,ip)*
     +     cos(twopi*bcdata(jj,ii,ip,1)*time)
         qk0(j,i,4,ip) = bcdata(jj,ii,ip,4)/qk0(j,i,1,ip)*
     +     cos(twopi*bcdata(jj,ii,ip,1)*time)
         qk0(j,i,5,ip) = q(j,1,i,5)
         else if (iflag .eq. 2038) then
         rnd_nbr=harvest(irnd_count)-0.5
         if (real(bcdata(jj,ii,ip,2)) .gt. 0.85*real(xmach)) rnd_nbr=0.
         if (real(bcdata(jj,ii,ip,2)) .gt. 0.75*real(xmach) .and.
     +       real(bcdata(jj,ii,ip,2)) .le. 0.85*real(xmach)) then
           rnd_nbr = (-10./xmach*bcdata(jj,ii,ip,2)+8.5)*rnd_nbr
         end if
         qk0(j,i,1,ip) = bcdata(jj,ii,ip,1)*(1.0+rnd_nbr*0.10)
         qk0(j,i,2,ip) = bcdata(jj,ii,ip,2)*(1.0+rnd_nbr*0.10)
         qk0(j,i,3,ip) = bcdata(jj,ii,ip,3)*(1.0+rnd_nbr*0.10)
         qk0(j,i,4,ip) = bcdata(jj,ii,ip,4)*(1.0+rnd_nbr*0.10)
         qk0(j,i,5,ip) = q(j,1,i,5)
         end if
         bck(j,i,1) = 0.0
  300    continue
      else
         do 3000 ip=1,2
         irnd_count=0
         do 3000 i=ista,iend1
         ii = i-ista+1
         do 3000 j=jsta,jend1
         jj = j-jsta+1
         irnd_count=irnd_count+1
         if (iflag .eq. 2008) then
         qk0(j,i,1,ip) =  bcdata(jj,ii,ip,1)
         qk0(j,i,2,ip) =  bcdata(jj,ii,ip,2)
         qk0(j,i,3,ip) = -bcdata(jj,ii,ip,4)
         qk0(j,i,4,ip) =  bcdata(jj,ii,ip,3)
         qk0(j,i,5,ip) =  q(j,1,i,5)
         else if (iflag .eq. 2018) then
         qk0(j,i,5,ip) =  q(j,1,i,5)
         qk0(j,i,1,ip) =  gamma*qk0(j,i,5,ip)/bcdata(jj,ii,ip,1)
         qk0(j,i,2,ip) =  bcdata(jj,ii,ip,2)/qk0(j,i,1,ip)
         qk0(j,i,3,ip) = -bcdata(jj,ii,ip,4)/qk0(j,i,1,ip)
         qk0(j,i,4,ip) =  bcdata(jj,ii,ip,3)/qk0(j,i,1,ip)
         else if (iflag .eq. 2028) then
         qk0(j,i,1,ip) =  q(j,1,i,1)
         qk0(j,i,2,ip) =  bcdata(jj,ii,ip,2)/qk0(j,i,1,ip)*
     +     cos(twopi*bcdata(jj,ii,ip,1)*time)
         qk0(j,i,3,ip) = -bcdata(jj,ii,ip,4)/qk0(j,i,1,ip)*
     +     cos(twopi*bcdata(jj,ii,ip,1)*time)
         qk0(j,i,4,ip) =  bcdata(jj,ii,ip,3)/qk0(j,i,1,ip)*
     +     cos(twopi*bcdata(jj,ii,ip,1)*time)
         qk0(j,i,5,ip) =  q(j,1,i,5)
         else if (iflag .eq. 2038) then
         rnd_nbr=harvest(irnd_count)-0.5
         if (real(bcdata(jj,ii,ip,2)) .gt. 0.85*real(xmach)) rnd_nbr=0.
         if (real(bcdata(jj,ii,ip,2)) .gt. 0.75*real(xmach) .and.
     +       real(bcdata(jj,ii,ip,2)) .le. 0.85*real(xmach)) then
           rnd_nbr = (-10./xmach*bcdata(jj,ii,ip,2)+8.5)*rnd_nbr
         end if
         qk0(j,i,1,ip) = bcdata(jj,ii,ip,1)*(1.0+rnd_nbr*0.10)
         qk0(j,i,2,ip) = bcdata(jj,ii,ip,2)*(1.0+rnd_nbr*0.10)
         qk0(j,i,3,ip) =-bcdata(jj,ii,ip,4)*(1.0+rnd_nbr*0.10)
         qk0(j,i,4,ip) = bcdata(jj,ii,ip,3)*(1.0+rnd_nbr*0.10)
         qk0(j,i,5,ip) = q(j,1,i,5)
         end if
         bck(j,i,1) = 0.0
 3000    continue
      end if
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
          t11 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(jj,ii,ip,4+l)
          ip  = 2
          t12 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(jj,ii,ip,4+l)
          tk0(j,i,l,1) = t11
          tk0(j,i,l,2) = t12
  301   continue
        enddo
      end if
      end if
      end if
c
c******************************************************************************
c      k=kdim boundary   set velocity and turb quantities           bctype 2008
c******************************************************************************
c
      if (nface.eq.6) then
c
c     check to see if turbulence data is input (itrflg1 = 1) or
c     if freestream values are to be used (itrflg1 = 0); the check
c     assumes if the first point has been set, all points have been
c
      ipp     = 1
      itrflg1 = 0
      if (real(bcdata(1,1,ipp,5)) .gt. -1.e10) itrflg1 = 1
c
      if (ialph.eq.0) then
         do 400 ip=1,2
         irnd_count=0
         do 400 i=ista,iend1
         ii = i-ista+1
         do 400 j=jsta,jend1
         jj = j-jsta+1
         irnd_count=irnd_count+1
         if (iflag .eq. 2008) then
         qk0(j,i,1,ip+2) = bcdata(jj,ii,ip,1)
         qk0(j,i,2,ip+2) = bcdata(jj,ii,ip,2)
         qk0(j,i,3,ip+2) = bcdata(jj,ii,ip,3)
         qk0(j,i,4,ip+2) = bcdata(jj,ii,ip,4)
         qk0(j,i,5,ip+2) = q(j,kdim1,i,5)
         else if (iflag .eq. 2018) then
         qk0(j,i,5,ip+2) = q(j,kdim1,i,5)
         qk0(j,i,1,ip+2) = gamma*qk0(j,i,5,ip+2)/bcdata(jj,ii,ip,1)
         qk0(j,i,2,ip+2) = bcdata(jj,ii,ip,2)/qk0(j,i,1,ip+2)
         qk0(j,i,3,ip+2) = bcdata(jj,ii,ip,3)/qk0(j,i,1,ip+2)
         qk0(j,i,4,ip+2) = bcdata(jj,ii,ip,4)/qk0(j,i,1,ip+2)
         else if (iflag .eq. 2028) then
         qk0(j,i,1,ip+2) = q(j,kdim1,i,1)
         qk0(j,i,2,ip+2) = bcdata(jj,ii,ip,2)/qk0(j,i,1,ip+2)*
     +     cos(twopi*bcdata(jj,ii,ip,1)*time)
         qk0(j,i,3,ip+2) = bcdata(jj,ii,ip,3)/qk0(j,i,1,ip+2)*
     +     cos(twopi*bcdata(jj,ii,ip,1)*time)
         qk0(j,i,4,ip+2) = bcdata(jj,ii,ip,4)/qk0(j,i,1,ip+2)*
     +     cos(twopi*bcdata(jj,ii,ip,1)*time)
         qk0(j,i,5,ip+2) = q(j,kdim1,i,5)
         else if (iflag .eq. 2038) then
         rnd_nbr=harvest(irnd_count)-0.5
         if (real(bcdata(jj,ii,ip,2)) .gt. 0.85*real(xmach)) rnd_nbr=0.
         if (real(bcdata(jj,ii,ip,2)) .gt. 0.75*real(xmach) .and.
     +       real(bcdata(jj,ii,ip,2)) .le. 0.85*real(xmach)) then
           rnd_nbr = (-10./xmach*bcdata(jj,ii,ip,2)+8.5)*rnd_nbr
         end if
         qk0(j,i,1,ip+2) = bcdata(jj,ii,ip,1)*(1.0+rnd_nbr*0.10)
         qk0(j,i,2,ip+2) = bcdata(jj,ii,ip,2)*(1.0+rnd_nbr*0.10)
         qk0(j,i,3,ip+2) = bcdata(jj,ii,ip,3)*(1.0+rnd_nbr*0.10)
         qk0(j,i,4,ip+2) = bcdata(jj,ii,ip,4)*(1.0+rnd_nbr*0.10)
         qk0(j,i,5,ip+2) = q(j,kdim1,i,5)
         end if
         bck(j,i,2) = 0.0
  400    continue
      else
         do 4000 ip=1,2
         irnd_count=0
         do 4000 i=ista,iend1
         ii = i-ista+1
         do 4000 j=jsta,jend1
         jj = j-jsta+1
         irnd_count=irnd_count+1
         if (iflag .eq. 2008) then
         qk0(j,i,1,ip+2) =  bcdata(jj,ii,ip,1)
         qk0(j,i,2,ip+2) =  bcdata(jj,ii,ip,2)
         qk0(j,i,3,ip+2) = -bcdata(jj,ii,ip,4)
         qk0(j,i,4,ip+2) =  bcdata(jj,ii,ip,3)
         qk0(j,i,5,ip+2) =  q(j,kdim1,i,5)
         else if (iflag .eq. 2018) then
         qk0(j,i,5,ip+2) =  q(j,kdim1,i,5)
         qk0(j,i,1,ip+2) =  gamma*qk0(j,i,5,ip+2)/bcdata(jj,ii,ip,1)
         qk0(j,i,2,ip+2) =  bcdata(jj,ii,ip,2)/qk0(j,i,1,ip+2)
         qk0(j,i,3,ip+2) = -bcdata(jj,ii,ip,4)/qk0(j,i,1,ip+2)
         qk0(j,i,4,ip+2) =  bcdata(jj,ii,ip,3)/qk0(j,i,1,ip+2)
         else if (iflag .eq. 2028) then
         qk0(j,i,1,ip+2) =  q(j,kdim1,i,1)
         qk0(j,i,2,ip+2) =  bcdata(jj,ii,ip,2)/qk0(j,i,1,ip+2)*
     +     cos(twopi*bcdata(jj,ii,ip,1)*time)
         qk0(j,i,3,ip+2) = -bcdata(jj,ii,ip,4)/qk0(j,i,1,ip+2)*
     +     cos(twopi*bcdata(jj,ii,ip,1)*time)
         qk0(j,i,4,ip+2) =  bcdata(jj,ii,ip,3)/qk0(j,i,1,ip+2)*
     +     cos(twopi*bcdata(jj,ii,ip,1)*time)
         qk0(j,i,5,ip+2) =  q(j,kdim1,i,5)
         else if (iflag .eq. 2038) then
         rnd_nbr=harvest(irnd_count)-0.5
         if (real(bcdata(jj,ii,ip,2)) .gt. 0.85*real(xmach)) rnd_nbr=0.
         if (real(bcdata(jj,ii,ip,2)) .gt. 0.75*real(xmach) .and.
     +       real(bcdata(jj,ii,ip,2)) .le. 0.85*real(xmach)) then
           rnd_nbr = (-10./xmach*bcdata(jj,ii,ip,2)+8.5)*rnd_nbr
         end if
         qk0(j,i,1,ip+2) = bcdata(jj,ii,ip,1)*(1.0+rnd_nbr*0.10)
         qk0(j,i,2,ip+2) = bcdata(jj,ii,ip,2)*(1.0+rnd_nbr*0.10)
         qk0(j,i,3,ip+2) =-bcdata(jj,ii,ip,4)*(1.0+rnd_nbr*0.10)
         qk0(j,i,4,ip+2) = bcdata(jj,ii,ip,3)*(1.0+rnd_nbr*0.10)
         qk0(j,i,5,ip+2) = q(j,kdim1,i,5)
         end if
         bck(j,i,2) = 0.0
 4000    continue
      end if
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
          t13 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(jj,ii,ip,4+l)
          ip  = 2
          t14 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(jj,ii,ip,4+l)
          tk0(j,i,l,3) = t13
          tk0(j,i,l,4) = t14
  401   continue
        enddo
      end if
      end if
      end if
c
c******************************************************************************
c      i=1 boundary      set velocity and turb quantities           bctype 2008
c******************************************************************************
c
      if (nface.eq.1) then
c
c     check to see if turbulence data is input (itrflg1 = 1) or
c     if freestream values are to be used (itrflg1 = 0); the check
c     assumes if the first point has been set, all points have been
c
      ipp     = 1
      itrflg1 = 0
      if (real(bcdata(1,1,ipp,5)) .gt. -1.e10) itrflg1 = 1
c
      if (ialph.eq.0) then
         do 500 ip=1,2
         irnd_count=0
         do 500 k=ksta,kend1
         kk = k-ksta+1
         do 500 j=jsta,jend1
         jj = j-jsta+1
         irnd_count=irnd_count+1
         if (iflag .eq. 2008) then
         qi0(j,k,1,ip) = bcdata(jj,kk,ip,1)
         qi0(j,k,2,ip) = bcdata(jj,kk,ip,2)
         qi0(j,k,3,ip) = bcdata(jj,kk,ip,3)
         qi0(j,k,4,ip) = bcdata(jj,kk,ip,4)
         qi0(j,k,5,ip) = q(j,k,1,5)
         else if (iflag .eq. 2018) then
         qi0(j,k,5,ip) = q(j,k,1,5)
         qi0(j,k,1,ip) = gamma*qi0(j,k,5,ip)/bcdata(jj,kk,ip,1)
         qi0(j,k,2,ip) = bcdata(jj,kk,ip,2)/qi0(j,k,1,ip)
         qi0(j,k,3,ip) = bcdata(jj,kk,ip,3)/qi0(j,k,1,ip)
         qi0(j,k,4,ip) = bcdata(jj,kk,ip,4)/qi0(j,k,1,ip)
         else if (iflag .eq. 2028) then
         qi0(j,k,1,ip) = q(j,k,1,1)
         qi0(j,k,2,ip) = bcdata(jj,kk,ip,2)/qi0(j,k,1,ip)*
     +     cos(twopi*bcdata(jj,kk,ip,1)*time)
         qi0(j,k,3,ip) = bcdata(jj,kk,ip,3)/qi0(j,k,1,ip)*
     +     cos(twopi*bcdata(jj,kk,ip,1)*time)
         qi0(j,k,4,ip) = bcdata(jj,kk,ip,4)/qi0(j,k,1,ip)*
     +     cos(twopi*bcdata(jj,kk,ip,1)*time)
         qi0(j,k,5,ip) = q(j,k,1,5)
         else if (iflag .eq. 2038) then
         rnd_nbr=harvest(irnd_count)-0.5
         if (real(bcdata(jj,kk,ip,2)) .gt. 0.85*real(xmach)) rnd_nbr=0.
         if (real(bcdata(jj,kk,ip,2)) .gt. 0.75*real(xmach) .and.
     +       real(bcdata(jj,kk,ip,2)) .le. 0.85*real(xmach)) then
           rnd_nbr = (-10./xmach*bcdata(jj,kk,ip,2)+8.5)*rnd_nbr
         end if
         qi0(j,k,1,ip) = bcdata(jj,kk,ip,1)*(1.0+rnd_nbr*0.10)
         qi0(j,k,2,ip) = bcdata(jj,kk,ip,2)*(1.0+rnd_nbr*0.10)
         qi0(j,k,3,ip) = bcdata(jj,kk,ip,3)*(1.0+rnd_nbr*0.10)
         qi0(j,k,4,ip) = bcdata(jj,kk,ip,4)*(1.0+rnd_nbr*0.10)
         qi0(j,k,5,ip) = q(j,k,1,5)
         end if
         bci(j,k,1) = 0.0
  500    continue
      else
         do 5000 ip=1,2
         irnd_count=0
         do 5000 k=ksta,kend1
         kk = k-ksta+1
         do 5000 j=jsta,jend1
         jj = j-jsta+1
         irnd_count=irnd_count+1
         if (iflag .eq. 2008) then
         qi0(j,k,1,ip) =  bcdata(jj,kk,ip,1)
         qi0(j,k,2,ip) =  bcdata(jj,kk,ip,2)
         qi0(j,k,3,ip) = -bcdata(jj,kk,ip,4)
         qi0(j,k,4,ip) =  bcdata(jj,kk,ip,3)
         qi0(j,k,5,ip) =  q(j,k,1,5)
         else if (iflag .eq. 2018) then
         qi0(j,k,5,ip) =  q(j,k,1,5)
         qi0(j,k,1,ip) =  gamma*qi0(j,k,5,ip)/bcdata(jj,kk,ip,1)
         qi0(j,k,2,ip) =  bcdata(jj,kk,ip,2)/qi0(j,k,1,ip)
         qi0(j,k,3,ip) = -bcdata(jj,kk,ip,4)/qi0(j,k,1,ip)
         qi0(j,k,4,ip) =  bcdata(jj,kk,ip,3)/qi0(j,k,1,ip)
         else if (iflag .eq. 2028) then
         qi0(j,k,1,ip) =  q(j,k,1,1)
         qi0(j,k,2,ip) =  bcdata(jj,kk,ip,2)/qi0(j,k,1,ip)*
     +     cos(twopi*bcdata(jj,kk,ip,1)*time)
         qi0(j,k,3,ip) = -bcdata(jj,kk,ip,4)/qi0(j,k,1,ip)*
     +     cos(twopi*bcdata(jj,kk,ip,1)*time)
         qi0(j,k,4,ip) =  bcdata(jj,kk,ip,3)/qi0(j,k,1,ip)*
     +     cos(twopi*bcdata(jj,kk,ip,1)*time)
         qi0(j,k,5,ip) =  q(j,k,1,5)
         else if (iflag .eq. 2038) then
         rnd_nbr=harvest(irnd_count)-0.5
         if (real(bcdata(jj,kk,ip,2)) .gt. 0.85*real(xmach)) rnd_nbr=0.
         if (real(bcdata(jj,kk,ip,2)) .gt. 0.75*real(xmach) .and.
     +       real(bcdata(jj,kk,ip,2)) .le. 0.85*real(xmach)) then
           rnd_nbr = (-10./xmach*bcdata(jj,kk,ip,2)+8.5)*rnd_nbr
         end if
         qi0(j,k,1,ip) = bcdata(jj,kk,ip,1)*(1.0+rnd_nbr*0.10)
         qi0(j,k,2,ip) = bcdata(jj,kk,ip,2)*(1.0+rnd_nbr*0.10)
         qi0(j,k,3,ip) =-bcdata(jj,kk,ip,4)*(1.0+rnd_nbr*0.10)
         qi0(j,k,4,ip) = bcdata(jj,kk,ip,3)*(1.0+rnd_nbr*0.10)
         qi0(j,k,5,ip) = q(j,k,1,5)
         end if
         bci(j,k,1) = 0.0
 5000    continue
      end if
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
          t11 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(jj,kk,ip,4+l)
          ip  = 2
          t12 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(jj,kk,ip,4+l)
          ti0(j,k,l,1) = t11
          ti0(j,k,l,2) = t12
  501   continue
        enddo
      end if
      end if
      end if
c
c******************************************************************************
c      i=idim boundary   set velocity and turb quantities           bctype 2008
c******************************************************************************
c
      if (nface.eq.2) then
c
c     check to see if turbulence data is input (itrflg1 = 1) or
c     if freestream values are to be used (itrflg1 = 0); the check
c     assumes if the first point has been set, all points have been
c
      ipp     = 1
      itrflg1 = 0
      if (real(bcdata(1,1,ipp,5)) .gt. -1.e10) itrflg1 = 1
c
      if (ialph.eq.0) then
         do 600 ip=1,2
         irnd_count=0
         do 600 k=ksta,kend1
         kk = k-ksta+1
         do 600 j=jsta,jend1
         jj = j-jsta+1
         irnd_count=irnd_count+1
         if (iflag .eq. 2008) then
         qi0(j,k,1,ip+2) = bcdata(jj,kk,ip,1)
         qi0(j,k,2,ip+2) = bcdata(jj,kk,ip,2)
         qi0(j,k,3,ip+2) = bcdata(jj,kk,ip,3)
         qi0(j,k,4,ip+2) = bcdata(jj,kk,ip,4)
         qi0(j,k,5,ip+2) = q(j,k,idim1,5)
         else if (iflag .eq. 2018) then
         qi0(j,k,5,ip+2) = q(j,k,idim1,5)
         qi0(j,k,1,ip+2) = gamma*qi0(j,k,5,ip+2)/bcdata(jj,kk,ip,1)
         qi0(j,k,2,ip+2) = bcdata(jj,kk,ip,2)/qi0(j,k,1,ip+2)
         qi0(j,k,3,ip+2) = bcdata(jj,kk,ip,3)/qi0(j,k,1,ip+2)
         qi0(j,k,4,ip+2) = bcdata(jj,kk,ip,4)/qi0(j,k,1,ip+2)
         else if (iflag .eq. 2028) then
         qi0(j,k,1,ip+2) = q(j,k,idim1,1)
         qi0(j,k,2,ip+2) = bcdata(jj,kk,ip,2)/qi0(j,k,1,ip+2)*
     +     cos(twopi*bcdata(jj,kk,ip,1)*time)
         qi0(j,k,3,ip+2) = bcdata(jj,kk,ip,3)/qi0(j,k,1,ip+2)*
     +     cos(twopi*bcdata(jj,kk,ip,1)*time)
         qi0(j,k,4,ip+2) = bcdata(jj,kk,ip,4)/qi0(j,k,1,ip+2)*
     +     cos(twopi*bcdata(jj,kk,ip,1)*time)
         qi0(j,k,5,ip+2) = q(j,k,idim1,5)
         else if (iflag .eq. 2038) then
         rnd_nbr=harvest(irnd_count)-0.5
         if (real(bcdata(jj,kk,ip,2)) .gt. 0.85*real(xmach)) rnd_nbr=0.
         if (real(bcdata(jj,kk,ip,2)) .gt. 0.75*real(xmach) .and.
     +       real(bcdata(jj,kk,ip,2)) .le. 0.85*real(xmach)) then
           rnd_nbr = (-10./xmach*bcdata(jj,kk,ip,2)+8.5)*rnd_nbr
         end if
         qi0(j,k,1,ip+2) = bcdata(jj,kk,ip,1)*(1.0+rnd_nbr*0.10)
         qi0(j,k,2,ip+2) = bcdata(jj,kk,ip,2)*(1.0+rnd_nbr*0.10)
         qi0(j,k,3,ip+2) = bcdata(jj,kk,ip,3)*(1.0+rnd_nbr*0.10)
         qi0(j,k,4,ip+2) = bcdata(jj,kk,ip,4)*(1.0+rnd_nbr*0.10)
         qi0(j,k,5,ip+2) = q(j,k,idim1,5)
         end if
         bci(j,k,2) = 0.0
  600    continue
      else
         do 6000 ip=1,2
         irnd_count=0
         do 6000 k=ksta,kend1
         kk = k-ksta+1
         do 6000 j=jsta,jend1
         jj = j-jsta+1
         irnd_count=irnd_count+1
         if (iflag .eq. 2008) then
         qi0(j,k,1,ip+2) =  bcdata(jj,kk,ip,1)
         qi0(j,k,2,ip+2) =  bcdata(jj,kk,ip,2)
         qi0(j,k,3,ip+2) = -bcdata(jj,kk,ip,4)
         qi0(j,k,4,ip+2) =  bcdata(jj,kk,ip,3)
         qi0(j,k,5,ip+2) =  q(j,k,idim1,5)
         else if (iflag .eq. 2018) then
         qi0(j,k,5,ip+2) =  q(j,k,idim1,5)
         qi0(j,k,1,ip+2) =  gamma*qi0(j,k,5,ip+2)/bcdata(jj,kk,ip,1)
         qi0(j,k,2,ip+2) =  bcdata(jj,kk,ip,2)/qi0(j,k,1,ip+2)
         qi0(j,k,3,ip+2) = -bcdata(jj,kk,ip,4)/qi0(j,k,1,ip+2)
         qi0(j,k,4,ip+2) =  bcdata(jj,kk,ip,3)/qi0(j,k,1,ip+2)
         else if (iflag .eq. 2028) then
         qi0(j,k,1,ip+2) =  q(j,k,idim1,1)
         qi0(j,k,2,ip+2) =  bcdata(jj,kk,ip,2)/qi0(j,k,1,ip+2)*
     +     cos(twopi*bcdata(jj,kk,ip,1)*time)
         qi0(j,k,3,ip+2) = -bcdata(jj,kk,ip,4)/qi0(j,k,1,ip+2)*
     +     cos(twopi*bcdata(jj,kk,ip,1)*time)
         qi0(j,k,4,ip+2) =  bcdata(jj,kk,ip,3)/qi0(j,k,1,ip+2)*
     +     cos(twopi*bcdata(jj,kk,ip,1)*time)
         qi0(j,k,5,ip+2) =  q(j,k,idim1,5)
         else if (iflag .eq. 2038) then
         rnd_nbr=harvest(irnd_count)-0.5
         if (real(bcdata(jj,kk,ip,2)) .gt. 0.85*real(xmach)) rnd_nbr=0.
         if (real(bcdata(jj,kk,ip,2)) .gt. 0.75*real(xmach) .and.
     +       real(bcdata(jj,kk,ip,2)) .le. 0.85*real(xmach)) then
           rnd_nbr = (-10./xmach*bcdata(jj,kk,ip,2)+8.5)*rnd_nbr
         end if
         qi0(j,k,1,ip+2) = bcdata(jj,kk,ip,1)*(1.0+rnd_nbr*0.10)
         qi0(j,k,2,ip+2) = bcdata(jj,kk,ip,2)*(1.0+rnd_nbr*0.10)
         qi0(j,k,3,ip+2) =-bcdata(jj,kk,ip,4)*(1.0+rnd_nbr*0.10)
         qi0(j,k,4,ip+2) = bcdata(jj,kk,ip,3)*(1.0+rnd_nbr*0.10)
         qi0(j,k,5,ip+2) = q(j,k,idim1,5)
         end if
         bci(j,k,2) = 0.0
 6000    continue
      end if
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
          t13 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(jj,kk,ip,4+l)
          ip  = 2
          t14 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(jj,kk,ip,4+l)
          ti0(j,k,l,3) = t13
          ti0(j,k,l,4) = t14
  601   continue
        enddo
      end if
      end if
      end if
c
      if (iflag .eq. 2038) then
        deallocate(harvest)
      end if
      return 
      end
