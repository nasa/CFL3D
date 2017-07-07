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
      subroutine add2x(q,qc,jdim,kdim,idim,jj2,kk2,ii2,q1,dq,wq,
     .                 wqj,wqjk,js,ks,is,je,ke,ie,ipass,nbl,nblc,
     .                 nou,bou,nbuf,ibufdim,ll,myid)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Interpolate the solution or the correction from a 
c     coarser mesh to a finer embedded mesh.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension q(jdim,kdim,idim,ll)
      dimension qc(jj2,kk2,ii2,ll),q1(jj2,kk2,ii2,ll)
      dimension dq(jdim,kdim,idim,ll), wq(jj2,kk2,ii2,ll),
     .          wqj(jdim,kk2,ii2),  wqjk(jdim,kdim,ii2)
c
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
c
c      interpolate solution from coarser mesh to 
c       finer embedded mesh (mode=0)
c
c      interpolate correction from coarser mesh to
c        finer embedded mesh (mode=1)
c     Note:  mode=1 must ONLY be used for ll=5 (primitive variables):
      if(mode .ne. 0 .and. ll .ne. 5) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),'(''mode must = 0 when ll .ne. 5 in'',
     .     '' addx'')')
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
c
c      jdim,kdim,idim  finer mesh
c      jj2,kk2,ii2     coarser mesh
c      js,ks,is        coarser mesh starting indices
c      je,ke,ie        coarser mesh ending indices
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
      jdim2 = jdim-2
      kdim2 = kdim-2
      idim2 = idim-2
      jjl   = jj2-1
      kkl   = kk2-1
      iil   = ii2-1
c
c      semi-coarsening / directional refinement
c
      nsi   = (idim-1)/(ie-is)
      ista  = is-1
      iend  = ie
      if (nsi.eq.1) then
         ista  = is
         iend  = ie-1
      end if
c
      nn = jdim*kdim
      do 5 i=1,idim1
      do 5 n=1,ll
cdir$ ivdep
      do 5 izz=1,nn
      dq(izz,1,i,n) = 0.e0
    5 continue
c
      if (ipass.eq.1) then
c
c      wq=qc     mode=0    coarser grid
c      wq=qc-q1  mode=1    coarser grid
c
      if (mode.eq.0) then
         nn = jj2*kk2
         do 10 i=1,ii2
         do 10 n=1,ll
cdir$ ivdep
         do 10 izz=1,nn
         wq(izz,1,i,n) = qc(izz,1,i,n)
   10    continue
      else
         nn = jj2*kk2
         do 20 i=1,ii2
         do 20 n=1,ll
cdir$ ivdep
         do 20 izz=1,nn
         wq(izz,1,i,n) = qc(izz,1,i,n)-q1(izz,1,i,n)
   20    continue
      end if
      end if
c
      do 80 n=1,ll
c
c      interpolate in j onto finer mesh
c
      ksta = max(ks-1,1)
      do 40 i=ista,iend
      do 40 k=ksta,ke
      jj = 0
      do 40 j=js,je-1
      jm = max(1,j-1)
      jp = min(jjl,j+1)
      jj = jj+1
      wqj(jj,k,i) = 0.75e0*wq(j,k,i,n)+0.25e0*wq(jm,k,i,n)
      jj = jj+1
      wqj(jj,k,i) = 0.75e0*wq(j,k,i,n)+0.25e0*wq(jp,k,i,n)
   40 continue
c
c      interpolate in k onto finer mesh
c
      do 50 i=ista,iend
      kk = 0
      do 50 k=ks,ke-1
      km = max(1,k-1)
      kp = min(kkl,k+1)
      kk = kk+1
      do 1006 j=1,jdim1
      wqjk(j,kk,i) = 0.75e0*wqj(j,k,i)+0.25e0*wqj(j,km,i)
 1006 continue
      kk = kk+1
      do 1007 j=1,jdim1
      wqjk(j,kk,i) = 0.75e0*wqj(j,k,i)+0.25e0*wqj(j,kp,i)
 1007 continue
   50 continue
c
c      interpolate in i
c
      np = jdim*kdim1-1
      if (nsi.eq.2) then
         ii = 0
         do 60 i=is,ie-1
         im = max(1,i-1)
        ip = min(iil,i+1)
         ii = ii+1
cdir$ ivdep
         do 1009 izz=1,np
         dq(izz,1,ii,n) = 0.75e0*wqjk(izz,1,i)+0.25e0*wqjk(izz,1,im)
 1009    continue
         ii = ii+1
cdir$ ivdep
         do 1010 izz=1,np
         dq(izz,1,ii,n) = 0.75e0*wqjk(izz,1,i)+0.25e0*wqjk(izz,1,ip)
 1010    continue
   60    continue
      else
         ii = 0
         do 65 i=is,ie-1
         ii = ii+1
cdir$ ivdep
         do 64 izz=1,np
         dq(izz,1,ii,n) = wqjk(izz,1,i)
   64    continue
   65    continue
      end if
c
   80 continue
c
c      q=dq    mode=0      fine grid    interpolation
c      q=q+dq  mode=1      fine grid    correction
c
      if (mode.eq.0) then
         do 105 i=1,idim1
         nn = jdim*kdim-jdim-1
         do 105 n=1,ll
cdir$ ivdep
         do 105 izz=1,nn
         q(izz,1,i,n) = dq(izz,1,i,n)
  105    continue
      else
         do 110 i=1,idim1
         nn = jdim*kdim-jdim-1
c
c         update density and pressure to ensure positivity
c
         alpq  = -.2
         phiq  = 1./0.5
         betq  = 1. + alpq*phiq
cdir$ ivdep
         do 7013 izz=1,nn
         t1            = dq(izz,1,i,1)/q(izz,1,i,1)
         t2            = dq(izz,1,i,1)/( betq + ccabs(t1)*phiq )
         dq(izz,1,i,1) =ccvmgt(t2,dq(izz,1,i,1),
     .                  (real(t1).lt.real(alpq)))
         t1            = dq(izz,1,i,5)/q(izz,1,i,5)
         t2            = dq(izz,1,i,5)/( betq + ccabs(t1)*phiq )
         dq(izz,1,i,5) =ccvmgt(t2,dq(izz,1,i,5),
     .                  (real(t1).lt.real(alpq)))
 7013    continue
c
c         update primitive variables
c
cdir$ ivdep
         do 1013 izz=1,nn
         q(izz,1,i,1) = q(izz,1,i,1)+dq(izz,1,i,1)
         q(izz,1,i,2) = q(izz,1,i,2)+dq(izz,1,i,2)
         q(izz,1,i,3) = q(izz,1,i,3)+dq(izz,1,i,3)
         q(izz,1,i,4) = q(izz,1,i,4)+dq(izz,1,i,4)
         q(izz,1,i,5) = q(izz,1,i,5)+dq(izz,1,i,5)
 1013    continue
  110    continue
      end if
      call fill(jdim,kdim,idim,q,ll)
      return
      end
