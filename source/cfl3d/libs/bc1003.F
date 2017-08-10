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
      subroutine bc1003(jdim,kdim,idim,q,qj0,qk0,qi0,sj,sk,si,bcj,bck,
     .              bci,ista,iend,jsta,jend,ksta,kend,t,jvdim,nface,
     .              tursav,tj0,tk0,ti0,vist3d,vj0,vk0,vi0,iuns,x,z,cl,
     .              nou,bou,nbuf,ibufdim,myid,nummem)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Set characteristic inflow/outflow boundary conditions
c               note: storage locations 21,22,23 in the t array are 
c               reserved for use in the subroutine rie1d
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension q(jdim,kdim,idim,5), qi0(jdim,kdim,5,4),
     .          qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
      dimension sk(jdim,kdim,idim-1,5),si(jdim,kdim,idim,5),
     .          sj(jdim,kdim,idim-1,5)
      dimension t(jvdim,23)
      dimension tursav(jdim,kdim,idim,nummem),tj0(kdim,idim-1,nummem,4),
     .          tk0(jdim,idim-1,nummem,4),ti0(jdim,kdim,nummem,4),
     .          vj0(kdim,idim-1,1,4),vk0(jdim,idim-1,1,4),
     .          vi0(jdim,kdim,1,4),vist3d(jdim,kdim,idim)
      dimension x(jdim,kdim,idim),z(jdim,kdim,idim)
c
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .        tur10(7)
      common /reyue/ reue,tinf,ivisc(3)
      common /sklton/ isklton
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
      jend1 = jend-1
      kend1 = kend-1
      iend1 = iend-1
c
c            * * * * * * * * * * * * * * * * * * * * * *
c            * standard boundary condition bctype=1003 *
c            * * * * * * * * * * * * * * * * * * * * * *
c
c******************************************************************************
c      j=1 boundary                inflow/outflow                   bctype 1003
c******************************************************************************
c
      if (nface.eq.3) then
      do 300 i=ista,iend1
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
  800 t(js+k-ksta,20) = -sj(1,k,i,5)
      if(iipv .eq. 1) then
      do k=ksta,kend1
        xi_ = 0.25*(x(1,k,i)+x(1,k+1,i) + x(2,k,i)+x(2,k+1,i))
        zi_ = 0.25*(z(1,k,i)+z(1,k+1,i) + z(2,k,i)+z(2,k+1,i))
        xo_ = 0.5*(x(1,k,i)+x(1,k+1,i))
        zo_ = 0.5*(z(1,k,i)+z(1,k+1,i))
        t(js+k-ksta,18) = xo_ + (xo_-xi_)
        t(js+k-ksta,19) = zo_ + (zo_-zi_)
      enddo
      end if
  300 continue
c
      jv = (kend-ksta)*(iend-ista)
      call rie1d(jvdim,t,jv,cl)
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
        do 101 k=ksta,kend1
          ubar=-(qj0(k,i,2,1)*sj(1,k,i,1)+qj0(k,i,3,1)*sj(1,k,i,2)+
     +           qj0(k,i,4,1)*sj(1,k,i,3))
          if (real(ubar) .lt. 0.) then
             tj0(k,i,l,1) = tur10(l)
             tj0(k,i,l,2) = tur10(l)
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
c      j=jdim boundary             inflow/outflow                   bctype 1003
c******************************************************************************
c
      if (nface.eq.4) then
      do 310 i=ista,iend1
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
  810 t(js+k-ksta,20) = sj(jdim,k,i,5)
      if (iipv .eq. 1) then
      do k=ksta,kend1
        xi_ = 0.25*(x(jdim ,k,i)+x(jdim,k+1,i)
     .            + x(jdim1,k,i)+x(jdim1,k+1,i))
        zi_ = 0.25*(z(jdim ,k,i)+z(jdim ,k+1,i)
     .            + z(jdim1,k,i)+z(jdim1,k+1,i))
        xo_ = 0.5*(x(jdim,k,i)+x(jdim,k+1,i))
        zo_ = 0.5*(z(jdim,k,i)+z(jdim,k+1,i))
        t(js+k-ksta,18) = xo_ + (xo_-xi_)
        t(js+k-ksta,19) = zo_ + (zo_-zi_)
      enddo
      end if
  310 continue
c
      jv = (kend-ksta)*(iend-ista)
      call rie1d(jvdim,t,jv,cl)
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
        do 201 k=ksta,kend1
          ubar=qj0(k,i,2,3)*sj(jdim,k,i,1)+qj0(k,i,3,3)*sj(jdim,k,i,2)+
     +         qj0(k,i,4,3)*sj(jdim,k,i,3)
          if (real(ubar) .lt. 0.) then
             tj0(k,i,l,3) = tur10(l)
             tj0(k,i,l,4) = tur10(l)
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
c      k=1 boundary                inflow/outflow                   bctype 1003
c******************************************************************************
c
      if (nface.eq.5) then
      do 320 i=ista,iend1
      js = (i-ista)*(jend-jsta)+1
      do 420 l=1,5
cdir$ ivdep
      do 520 jj=jsta,jend1
  520 t(js+jj-jsta,l) = q(jj,1,i,l)
  420 continue
      do 620 l=1,3
cdir$ ivdep
      do 720 jj=jsta,jend1
  720 t(js+jj-jsta,5+l) = -sk(jj,1,i,l)
  620 continue
cdir$ ivdep
      do 820 jj=jsta,jend1
  820 t(js+jj-jsta,20) = -sk(jj,1,i,5)
      if (iipv .eq. 1) then
      do jj=jsta,jend1
        xi_ = 0.25*(x(jj,1,i)+x(jj+1,1,i) + x(jj,2,i)+x(jj+1,2,i))
        zi_ = 0.25*(z(jj,1,i)+z(jj+1,1,i) + z(jj,2,i)+z(jj+1,2,i))
        xo_ = 0.5*(x(jj,1,i)+x(jj+1,1,i))
        zo_ = 0.5*(z(jj,1,i)+z(jj+1,1,i))
        t(js+jj-jsta,18) = xo_ + (xo_-xi_)
        t(js+jj-jsta,19) = zo_ + (zo_-zi_)
      enddo
      end if
  320 continue
c
      jv = (jend-jsta)*(iend-ista)
      call rie1d(jvdim,t,jv,cl)
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
        do 301 j=jsta,jend1
          ubar=-(qk0(j,i,2,1)*sk(j,1,i,1)+qk0(j,i,3,1)*sk(j,1,i,2)+
     +           qk0(j,i,4,1)*sk(j,1,i,3))
          if (real(ubar) .lt. 0.) then
             tk0(j,i,l,1) = tur10(l)
             tk0(j,i,l,2) = tur10(l)
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
c      k=kdim boundary             inflow/outflow                   bctype 1003
c******************************************************************************
c
      if (nface.eq.6) then
      do 330 i=ista,iend1
      js = (i-ista)*(jend-jsta)+1
      do 430 l=1,5
cdir$ ivdep
      do 530 jj=jsta,jend1
  530 t(js+jj-jsta,l) = q(jj,kdim1,i,l)
  430 continue
      do 630 l=1,3
cdir$ ivdep
      do 730 jj=jsta,jend1
  730 t(js+jj-jsta,5+l) = sk(jj,kdim,i,l)
  630 continue
cdir$ ivdep
      do 830 jj=jsta,jend1
  830 t(js+jj-jsta,20) = sk(jj,kdim,i,5)
      if (iipv .eq. 1) then
      do jj=jsta,jend1
        xi_ = 0.25*(x(jj,kdim ,i)+x(jj+1,kdim ,i)
     .            + x(jj,kdim1,i)+x(jj+1,kdim1,i))
        zi_ = 0.25*(z(jj,kdim ,i)+z(jj+1,kdim ,i)
     .            + z(jj,kdim1,i)+z(jj+1,kdim1,i))
        xo_ = 0.5*(x(jj,kdim,i)+x(jj+1,kdim,i))
        zo_ = 0.5*(z(jj,kdim,i)+z(jj+1,kdim,i))
        t(js+jj-jsta,18) = xo_ + (xo_-xi_)
        t(js+jj-jsta,19) = zo_ + (zo_-zi_)
      enddo
      end if
  330 continue
c
      jv = (jend-jsta)*(iend-ista)
      call rie1d(jvdim,t,jv,cl)
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
        do 401 j=jsta,jend1
          ubar=qk0(j,i,2,3)*sk(j,kdim,i,1)+qk0(j,i,3,3)*sk(j,kdim,i,2)+
     +         qk0(j,i,4,3)*sk(j,kdim,i,3)
          if (real(ubar) .lt. 0.) then
             tk0(j,i,l,3) = tur10(l)
             tk0(j,i,l,4) = tur10(l)
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
c      i=1 boundary                inflow/outflow                   bctype 1003
c******************************************************************************
c
      if (nface.eq.1) then
      do 340 k=ksta,kend1
      js = (k-ksta)*(jend-jsta)+1
      do 440 l=1,5
cdir$ ivdep
      do 540 jj=jsta,jend1
  540 t(js+jj-jsta,l) = q(jj,k,1,l)
  440 continue
      do 640 l=1,3
cdir$ ivdep
      do 740 jj=jsta,jend1
  740 t(js+jj-jsta,5+l) = -si(jj,k,1,l)
  640 continue
c
      do 840 jj=jsta,jend1
  840 t(js+jj-jsta,20) = -si(jj,k,1,5)
  340 continue
c
      jv = (jend-jsta)*(kend-ksta)
      call rie1d(jvdim,t,jv,cl)
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
        do 501 j=jsta,jend1
          ubar=-(qi0(j,k,2,1)*si(j,k,1,1)+qi0(j,k,3,1)*si(j,k,1,2)+
     +           qi0(j,k,4,1)*si(j,k,1,3))
          if (real(ubar) .lt. 0.) then
             ti0(j,k,l,1) = tur10(l)
             ti0(j,k,l,2) = tur10(l)
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
c      i=idim boundary             inflow/outflow                   bctype 1003
c******************************************************************************
c
      if (nface.eq.2) then
      do 350 k=ksta,kend1
      js = (k-ksta)*(jend-jsta)+1
      do 450 l=1,5
cdir$ ivdep
      do 550 jj=jsta,jend1
  550 t(js+jj-jsta,l) = q(jj,k,idim1,l)
  450 continue
      do 650 l=1,3
cdir$ ivdep
      do 750 jj=jsta,jend1
  750 t(js+jj-jsta,5+l) = si(jj,k,idim,l)
  650 continue
cdir$ ivdep
      do 850 jj=jsta,jend1
  850 t(js+jj-jsta,20) = si(jj,k,idim,5)
  350 continue
c
      jv = (jend-jsta)*(kend-ksta)
      call rie1d(jvdim,t,jv,cl)
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
        do 601 j=jsta,jend1
          ubar=qi0(j,k,2,3)*si(j,k,idim,1)+qi0(j,k,3,3)*si(j,k,idim,2)+
     +         qi0(j,k,4,3)*si(j,k,idim,3)
          if (real(ubar) .lt. 0.) then
             ti0(j,k,l,3) = tur10(l)
             ti0(j,k,l,4) = tur10(l)
          else
             ti0(j,k,l,3) = tursav(j,k,idim-1,l)
             ti0(j,k,l,4) = ti0(j,k,l,3)
          end if
  601   continue
        enddo
      end if
      end if
      end if
c
      return
      end
