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
      subroutine bc2003(jdim,kdim,idim,q,qj0,qk0,qi0,sj,sk,si,bcj,bck,
     .                  bci,ista,iend,jsta,jend,ksta,kend,t,jvdim,nface,
     .                  tursav,tj0,tk0,ti0,vist3d,vj0,vk0,vi0,
     .                  mdim,ndim,bcdata,filname,iuns,
     .                  nou,bou,nbuf,ibufdim,myid,nummem)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Set characteristic inlet boundary conditions at engine 
c               inlet, given (estimated) inlet mach, total pressure, 
c               ratio, total temperature ratio, and flow angle, and up
c               to two turbulence quantities.
c
c     xme     = estimated engine inflow Mach number (note: it may end
c               up slightly different for converged solution)
c     pte     = total pressure ratio at inlet (ptotal/pinf)
c     tte     = total temperature ratio at inlet (Ttotal/Tinf)
c     alpe    = alpha angle of inlet flow, deg
c     betae   = beta angle of inlet flow, deg
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
      dimension t(jvdim,25),bcdata(mdim,ndim,2,12)
      dimension tursav(jdim,kdim,idim,nummem),tj0(kdim,idim-1,nummem,4),
     .          tk0(jdim,idim-1,nummem,4),ti0(jdim,kdim,nummem,4),
     .          vj0(kdim,idim-1,1,4),vk0(jdim,idim-1,1,4),
     .          vi0(jdim,kdim,1,4),vist3d(jdim,kdim,idim)
c
      common /maxiv/ ivmx
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .        tur10(7)
      common /reyue/ reue,tinf,ivisc(3)
      common /sklton/ isklton
      common /conversion/ radtodeg
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
      jend1 = jend-1
      kend1 = kend-1
      iend1 = iend-1
c
c     this bc makes use of only one plane of data    
c
      ip    = 1
c
c            * * * * * * * * * * * * * * * * * * * * * *
c            * standard boundary condition bctype=2003 *
c            * * * * * * * * * * * * * * * * * * * * * *
c
c******************************************************************************
c      j=1 boundary             inflow for engine                   bctype 2003
c******************************************************************************
c
      if (nface.eq.3) then
c
c     check to see if turbulence data is input (itrflg1 = 1) or
c     if freestream values are to be used (itrflg1 = 0); the check
c     assumes if the first point has been set, all points have been
c
      itrflg1 = 0
      if (real(bcdata(1,1,ip,6)) .gt. -1.e10) itrflg1 = 1
c
      do 300 i=ista,iend1
      ii = i-ista+1
      js = (i-ista)*(kend-ksta)+1
      do 400 l=1,5
c
      do 500 k=ksta,kend1
  500 t(js+k-ksta,l) = q(1,k,i,l)
  400 continue
      do 600 l=1,3
c
      do 700 k=ksta,kend1
  700 t(js+k-ksta,5+l) = -sj(1,k,i,l)
  600 continue
c
      do 800 k=ksta,kend1
      kk = k-ksta+1
c
      xme   = bcdata(kk,ii,ip,1)
      pte   = bcdata(kk,ii,ip,2)
      tte   = bcdata(kk,ii,ip,3)
      alpe  = bcdata(kk,ii,ip,4)
      betae = bcdata(kk,ii,ip,5)
c
      alpe=alpe/radtodeg
      betae=betae/radtodeg
      te=tte/(1.+0.5*gm1*xme*xme)
      ae=sqrt(te)
      ue=ae*xme*cos(alpe)*cos(betae)
      ve=-ae*xme*sin(betae)
      we=ae*xme*sin(alpe)*cos(betae)
c
      t(js+k-ksta,21) = ue
      t(js+k-ksta,22) = ve
      t(js+k-ksta,23) = we
      t(js+k-ksta,24) = ae
      t(js+k-ksta,25) = pte
      t(js+k-ksta,20) = -sj(1,k,i,5)
c
  800 continue
  300 continue
c
      jv = (kend-ksta)*(iend-ista)
      call rie1de(jvdim,t,jv)
c
      do 900 i=ista,iend1
      js = (i-ista)*(kend-ksta)+1
      do 900 l=1,5
      do 900 k=ksta,kend1
      qj0(k,i,l,1) = t(k-ksta+js,l)
      qj0(k,i,l,2) = qj0(k,i,l,1)
      bcj(k,i,1)   = 0.0
 900  continue
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 191 i=ista,iend1
        do 191 k=ksta,kend1
          vj0(k,i,1,1) = vist3d(1,k,i)
          vj0(k,i,1,2) = 0.0
  191   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 101 i=ista,iend1
        ii = i-ista+1
        do 101 k=ksta,kend1
          kk=k-ksta+1
          ubar=-(qj0(k,i,2,1)*sj(1,k,i,1)+qj0(k,i,3,1)*sj(1,k,i,2)+
     +           qj0(k,i,4,1)*sj(1,k,i,3))
          t1 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(kk,ii,ip,5+l)
          if (real(ubar) .lt. 0.) then
             tj0(k,i,l,1) = t1
             tj0(k,i,l,2) = t1
          else
             tj0(k,i,l,1) = tursav(1,k,i,l)
             tj0(k,i,l,2) = tj0(k,i,l,1)
          end if
  101   continue
        enddo
      end if
      end if
c
      end if
c
c******************************************************************************
c      j=jdim boundary          inflow for engine                   bctype 2003
c******************************************************************************
c
      if (nface.eq.4) then
c
c     check to see if turbulence data is input (itrflg1 = 1) or
c     if freestream values are to be used (itrflg1 = 0); the check
c     assumes if the first point has been set, all points have been
c
      itrflg1 = 0
      if (real(bcdata(1,1,ip,6)) .gt. -1.e10) itrflg1 = 1
c
      do 310 i=ista,iend1
      ii = i-ista+1
      js = (i-ista)*(kend-ksta)+1
      do 410 l=1,5
c
      do 510 k=ksta,kend1
  510 t(js+k-ksta,l) = q(jdim1,k,i,l)
  410 continue
      do 610 l=1,3
c
      do 710 k=ksta,kend1
  710 t(js+k-ksta,5+l) = sj(jdim,k,i,l)
  610 continue
c
      do 810 k=ksta,kend1
      kk = k-ksta+1
c
      xme   = bcdata(kk,ii,ip,1)
      pte   = bcdata(kk,ii,ip,2)
      tte   = bcdata(kk,ii,ip,3)
      alpe  = bcdata(kk,ii,ip,4)
      betae = bcdata(kk,ii,ip,5)
c
      alpe=alpe/radtodeg
      betae=betae/radtodeg
      te=tte/(1.+0.5*gm1*xme*xme)
      ae=sqrt(te)
      ue=ae*xme*cos(alpe)*cos(betae)
      ve=-ae*xme*sin(betae)
      we=ae*xme*sin(alpe)*cos(betae)
c
      t(js+k-ksta,21) = ue
      t(js+k-ksta,22) = ve
      t(js+k-ksta,23) = we
      t(js+k-ksta,24) = ae
      t(js+k-ksta,25) = pte
      t(js+k-ksta,20) = sj(jdim,k,i,5)
c
  810 continue
  310 continue
c
      jv = (kend-ksta)*(iend-ista)
      call rie1de(jvdim,t,jv)
c
      do 910 i=ista,iend1
      js = (i-ista)*(kend-ksta)+1
      do 910 l=1,5
      do 910 k=ksta,kend1
      qj0(k,i,l,3) = t(k-ksta+js,l)
      qj0(k,i,l,4) = qj0(k,i,l,3)
      bcj(k,i,2)   = 0.0
  910 continue
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 291 i=ista,iend1
        do 291 k=ksta,kend1
          vj0(k,i,1,3) = vist3d(jdim-1,k,i)
          vj0(k,i,1,4) = 0.0
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
          ubar=qj0(k,i,2,3)*sj(jdim,k,i,1)+qj0(k,i,3,3)*sj(jdim,k,i,2)+
     +         qj0(k,i,4,3)*sj(jdim,k,i,3)
          t1 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(kk,ii,ip,5+l)
          if (real(ubar) .lt. 0.) then
             tj0(k,i,l,3) = t1
             tj0(k,i,l,4) = t1
          else
             tj0(k,i,l,3) = tursav(jdim-1,k,i,l)
             tj0(k,i,l,4) = tj0(k,i,l,3)
          end if
  201   continue
        enddo
      end if
      end if
c
      end if
c
c******************************************************************************
c      k=1 boundary             inflow for engine                   bctype 2003
c******************************************************************************
c
      if (nface.eq.5) then
c
c     check to see if turbulence data is input (itrflg1 = 1) or
c     if freestream values are to be used (itrflg1 = 0); the check
c     assumes if the first point has been set, all points have been
c
      itrflg1 = 0
      if (real(bcdata(1,1,ip,6)) .gt. -1.e10) itrflg1 = 1
c
      do 320 i=ista,iend1
      ii = i-ista+1
      js = (i-ista)*(jend-jsta)+1
      do 420 l=1,5
cdir$ ivdep
      do 520 j=jsta,jend1
  520 t(js+j-jsta,l) = q(j,1,i,l)
  420 continue
      do 620 l=1,3
cdir$ ivdep
      do 720 j=jsta,jend1
  720 t(js+j-jsta,5+l) = -sk(j,1,i,l)
  620 continue
cdir$ ivdep
      do 820 j=jsta,jend1
      jj = j-jsta+1
c
      xme   = bcdata(jj,ii,ip,1)
      pte   = bcdata(jj,ii,ip,2)
      tte   = bcdata(jj,ii,ip,3)
      alpe  = bcdata(jj,ii,ip,4)
      betae = bcdata(jj,ii,ip,5)
c
      alpe=alpe/radtodeg
      betae=betae/radtodeg
      te=tte/(1.+0.5*gm1*xme*xme)
      ae=sqrt(te)
      ue=ae*xme*cos(alpe)*cos(betae)
      ve=-ae*xme*sin(betae)
      we=ae*xme*sin(alpe)*cos(betae)
c
      t(js+j-jsta,21) = ue
      t(js+j-jsta,22) = ve
      t(js+j-jsta,23) = we
      t(js+j-jsta,24) = ae
      t(js+j-jsta,25) = pte
      t(js+j-jsta,20) = -sk(j,1,i,5)
c
  820 continue
  320 continue
c
      jv = (jend-jsta)*(iend-ista)
      call rie1de(jvdim,t,jv)
c
      do 920 i=ista,iend1
      js = (i-ista)*(jend-jsta)+1
      do 920 l=1,5
      do 920 j=jsta,jend1
      qk0(j,i,l,1) = t(j-jsta+js,l)
      qk0(j,i,l,2) = qk0(j,i,l,1)
      bck(j,i,1)   = 0.0
 920  continue
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 391 i=ista,iend1
        do 391 j=jsta,jend1
          vk0(j,i,1,1) = vist3d(j,1,i)
          vk0(j,i,1,2) = 0.0
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
          ubar=-(qk0(j,i,2,1)*sk(j,1,i,1)+qk0(j,i,3,1)*sk(j,1,i,2)+
     +           qk0(j,i,4,1)*sk(j,1,i,3))
          t1 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(jj,ii,ip,5+l)
          if (real(ubar) .lt. 0.) then
             tk0(j,i,l,1) = t1
             tk0(j,i,l,2) = t1
          else
             tk0(j,i,l,1) = tursav(j,1,i,l)
             tk0(j,i,l,2) = tk0(j,i,l,1)
          end if
  301   continue
        enddo
      end if
      end if
c
      end if
c
c******************************************************************************
c      k=kdim boundary          inflow for engine                   bctype 2003
c******************************************************************************
c
      if (nface.eq.6) then
c
c     check to see if turbulence data is input (itrflg1 = 1) or
c     if freestream values are to be used (itrflg1 = 0); the check
c     assumes if the first point has been set, all points have been
c
      itrflg1 = 0
      if (real(bcdata(1,1,ip,6)) .gt. -1.e10) itrflg1 = 1
c
      do 330 i=ista,iend1
      ii = i-ista+1
      js = (i-ista)*(jend-jsta)+1
      do 430 l=1,5
cdir$ ivdep
      do 530 j=jsta,jend1
  530 t(js+j-jsta,l) = q(j,kdim1,i,l)
  430 continue
      do 630 l=1,3
cdir$ ivdep
      do 730 j=jsta,jend1
  730 t(js+j-jsta,5+l) = sk(j,kdim,i,l)
  630 continue
cdir$ ivdep
      do 830 j=jsta,jend1
      jj = j-jsta+1
c
      xme   = bcdata(jj,ii,ip,1)
      pte   = bcdata(jj,ii,ip,2)
      tte   = bcdata(jj,ii,ip,3)
      alpe  = bcdata(jj,ii,ip,4)
      betae = bcdata(jj,ii,ip,5)
c
      alpe=alpe/radtodeg
      betae=betae/radtodeg
      te=tte/(1.+0.5*gm1*xme*xme)
      ae=sqrt(te)
      ue=ae*xme*cos(alpe)*cos(betae)
      ve=-ae*xme*sin(betae)
      we=ae*xme*sin(alpe)*cos(betae)
c
      t(js+j-jsta,21) = ue
      t(js+j-jsta,22) = ve
      t(js+j-jsta,23) = we
      t(js+j-jsta,24) = ae
      t(js+j-jsta,25) = pte
      t(js+j-jsta,20) = sk(j,kdim,i,5)
c
  830 continue
  330 continue
c
      jv = (jend-jsta)*(iend-ista)
      call rie1de(jvdim,t,jv)
c
      do 930 i=ista,iend1
      js = (i-ista)*(jend-jsta)+1
      do 930 l=1,5
      do 930 j=jsta,jend1
      qk0(j,i,l,3) = t(j-jsta+js,l)
      qk0(j,i,l,4) = qk0(j,i,l,3)
      bck(j,i,2)   = 0.0
  930 continue
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 491 i=ista,iend1
        do 491 j=jsta,jend1
          vk0(j,i,1,3) = vist3d(j,kdim-1,i)
          vk0(j,i,1,4) = 0.0
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
          ubar=qk0(j,i,2,3)*sk(j,kdim,i,1)+qk0(j,i,3,3)*sk(j,kdim,i,2)+
     +         qk0(j,i,4,3)*sk(j,kdim,i,3)
          t1 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(jj,ii,ip,5+l)
          if (real(ubar) .lt. 0.) then
             tk0(j,i,l,3) = t1
             tk0(j,i,l,4) = t1
          else
             tk0(j,i,l,3) = tursav(j,kdim-1,i,l)
             tk0(j,i,l,4) = tk0(j,i,l,3)
          end if
  401   continue
        enddo
      end if
      end if
c
      end if
c
c******************************************************************************
c      i=1 boundary             inflow for engine                   bctype 2003
c******************************************************************************
c
      if (nface.eq.1) then
c
c     check to see if turbulence data is input (itrflg1 = 1) or
c     if freestream values are to be used (itrflg1 = 0); the check
c     assumes if the first point has been set, all points have been
c
      itrflg1 = 0
      if (real(bcdata(1,1,ip,6)) .gt. -1.e10) itrflg1 = 1
c
      do 340 k=ksta,kend1
      kk = k-ksta+1
      js = (k-ksta)*(jend-jsta)+1
      do 440 l=1,5
cdir$ ivdep
      do 540 j=jsta,jend1
  540 t(js+j-jsta,l) = q(j,k,1,l)
  440 continue
      do 640 l=1,3
cdir$ ivdep
      do 740 j=jsta,jend1
  740 t(js+j-jsta,5+l) = -si(j,k,1,l)
  640 continue
c
      do 840 j=jsta,jend1
      jj = j-jsta+1
c
      xme   = bcdata(jj,kk,ip,1)
      pte   = bcdata(jj,kk,ip,2)
      tte   = bcdata(jj,kk,ip,3)
      alpe  = bcdata(jj,kk,ip,4)
      betae = bcdata(jj,kk,ip,5)
c
      alpe=alpe/radtodeg
      betae=betae/radtodeg
      te=tte/(1.+0.5*gm1*xme*xme)
      ae=sqrt(te)
      ue=ae*xme*cos(alpe)*cos(betae)
      ve=-ae*xme*sin(betae)
      we=ae*xme*sin(alpe)*cos(betae)
c
      t(js+j-jsta,21) = ue
      t(js+j-jsta,22) = ve
      t(js+j-jsta,23) = we
      t(js+j-jsta,24) = ae
      t(js+j-jsta,25) = pte
      t(js+j-jsta,20) = -si(j,k,1,5)
c
  840 continue
  340 continue
c
      jv = (jend-jsta)*(kend-ksta)
      call rie1de(jvdim,t,jv)
c
      do 940 k=ksta,kend1
      js = (k-ksta)*(jend-jsta)+1
      do 940 l=1,5
      do 940 j=jsta,jend1
      qi0(j,k,l,1) = t(j-jsta+js,l)
      qi0(j,k,l,2) = qi0(j,k,l,1)
      bci(j,k,1)   = 0.0
 940  continue
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 591 k=ksta,kend1
        do 591 j=jsta,jend1
          vi0(j,k,1,1) = vist3d(j,k,1)
          vi0(j,k,1,2) = 0.0
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
          ubar=-(qi0(j,k,2,1)*si(j,k,1,1)+qi0(j,k,3,1)*si(j,k,1,2)+
     +           qi0(j,k,4,1)*si(j,k,1,3))
          t1 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(jj,kk,ip,5+l)
          if (real(ubar) .lt. 0.) then
             ti0(j,k,l,1) = t1
             ti0(j,k,l,2) = t1
          else
             ti0(j,k,l,1) = tursav(j,k,1,l)
             ti0(j,k,l,2) = ti0(j,k,l,1)
          end if
  501   continue
        enddo
      end if
      end if
c
      end if
c
c******************************************************************************
c      i=idim boundary          inflow for engine                   bctype 2003
c******************************************************************************
c
      if (nface.eq.2) then
c
c     check to see if turbulence data is input (itrflg1 = 1) or
c     if freestream values are to be used (itrflg1 = 0); the check
c     assumes if the first point has been set, all points have been
c
      itrflg1 = 0
      if (real(bcdata(1,1,ip,6)) .gt. -1.e10) itrflg1 = 1
c
      do 350 k=ksta,kend1
      kk = k-ksta+1
      js = (k-ksta)*(jend-jsta)+1
      do 450 l=1,5
cdir$ ivdep
      do 550 j=jsta,jend1
  550 t(js+j-jsta,l) = q(j,k,idim1,l)
  450 continue
      do 650 l=1,3
cdir$ ivdep
      do 750 j=jsta,jend1
  750 t(js+j-jsta,5+l) = si(j,k,idim,l)
  650 continue
cdir$ ivdep
      do 850 j=jsta,jend1
      jj = j-jsta+1
c
      xme   = bcdata(jj,kk,ip,1)
      pte   = bcdata(jj,kk,ip,2)
      tte   = bcdata(jj,kk,ip,3)
      alpe  = bcdata(jj,kk,ip,4)
      betae = bcdata(jj,kk,ip,5)
c
      alpe=alpe/radtodeg
      betae=betae/radtodeg
      te=tte/(1.+0.5*gm1*xme*xme)
      ae=sqrt(te)
      ue=ae*xme*cos(alpe)*cos(betae)
      ve=-ae*xme*sin(betae)
      we=ae*xme*sin(alpe)*cos(betae)
c
      t(js+j-jsta,21) = ue
      t(js+j-jsta,22) = ve
      t(js+j-jsta,23) = we
      t(js+j-jsta,24) = ae
      t(js+j-jsta,25) = pte
      t(js+j-jsta,20) = si(j,k,idim,5)
c
  850 continue
  350 continue
c
      jv = (jend-jsta)*(kend-ksta)
      call rie1de(jvdim,t,jv)
c
      do 950 k=ksta,kend1
      js = (k-ksta)*(jend-jsta)+1
      do 950 l=1,5
      do 950 j=jsta,jend1
      qi0(j,k,l,3) = t(j-jsta+js,l)
      qi0(j,k,l,4) = qi0(j,k,l,3)
      bci(j,k,2)   = 0.0
  950 continue
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 691 k=ksta,kend1
        do 691 j=jsta,jend1
          vi0(j,k,1,3) = vist3d(j,k,idim-1)
          vi0(j,k,1,4) = 0.0
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
          ubar=qi0(j,k,2,3)*si(j,k,idim,1)+qi0(j,k,3,3)*si(j,k,idim,2)+
     +         qi0(j,k,4,3)*si(j,k,idim,3)
          t1 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(jj,kk,ip,5+l)
          if (real(ubar) .lt. 0.) then
             ti0(j,k,l,3) = t1
             ti0(j,k,l,4) = t1
          else
             ti0(j,k,l,3) = tursav(j,k,idim-1,l)
             ti0(j,k,l,4) = ti0(j,k,l,3)
          end if
  601   continue
        enddo
      end if
      end if
c
      end if
c
      return
      end
