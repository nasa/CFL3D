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
      subroutine addx(q,qq,jdim,kdim,idim,jj2,kk2,ii2,q1,dq,wq,wqj,nbl,
     .                blank,nou,bou,nbuf,ibufdim,ll,myid)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Interpolate the solution or the correction from a 
c     coarser mesh to a finer mesh.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension q(jdim,kdim,idim,ll),qq(jj2,kk2,ii2,ll),
     .          q1(jj2,kk2,ii2,ll)
      dimension dq(jdim,kdim,idim,ll),wq(jj2,kk2,ii2,ll),
     .          wqj(jdim,kk2,ii2,ll),blank(jdim,kdim,idim)
c
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /mgv/ epsssc(3),epsssr(3),issc,issr
c
c      interpolate solution to finer mesh (imode=0)
c
c      interpolate correction to finer mesh (imode=1)
c
c     Note:  mode=1 must ONLY be used for ll=5 (primitive variables):
      if(mode .ne. 0 .and. ll .ne. 5) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),'(''mode must = 0 when ll .ne. 5 in'',
     .     '' addx'')')
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
c      jdim,kdim,idim  finer mesh indices
c      jj2,kk2,ii2   coarser mesh indices
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
      nplq  = min(idim1,999000/(jdim*kdim))
      npl   = nplq
      do 5 i=1,idim1,nplq
      if (i+npl-1.gt.idim1) npl = idim1+1-i
      nnpl  = npl*jdim*kdim
      do 5 n=1,ll
cdir$ ivdep
      do 1000 izz=1,nnpl
      dq(izz,1,i,n) = 0.e0
 1000 continue
    5 continue
c
      nbl1 = nbl-1
      if (mode.eq.0) then
      nplq = min(iil,999000/(jj2*kk2))
      npl  = nplq
      do 10 i=1,iil,nplq
      if (i+npl-1.gt.iil) npl = iil+1-i
      nnpl = npl*jj2*kk2-jj2-1
      do 10 n=1,ll
cdir$ ivdep
      do 1001 izz=1,nnpl
      wq(izz,1,i,n) = qq(izz,1,i,n)
 1001 continue
   10 continue
c
      else
c
      nplq = min(iil,999000/(jj2*kk2))
      npl  = nplq
      do 20 i=1,iil,nplq
      if (i+npl-1.gt.iil) npl = iil+1-i
      nnpl = npl*jj2*kk2-jj2-1
      do 20 n=1,ll
cdir$ ivdep
      do 1002 izz=1,nnpl
      wq(izz,1,i,n) = qq(izz,1,i,n)-q1(izz,1,i,n)
 1002 continue
   20 continue
c
c      smooth corrections
c
      if (issc.gt.0) call rsmooth(epsssc,ii2,jj2,kk2,ii2,
     .                            wq,wqj,nou,bou,nbuf,ibufdim)
      end if
c
c      interpolate in j
c
      jjl1 = jjl-1
      kkl1 = kkl-1
      do 80 n=1,ll
      do 40 i=1,iil
      do 40 k=1,kkl
      wqj(1,k,i,n)     = wq(1,k,i,n)
      wqj(jdim1,k,i,n) = wq(jjl,k,i,n)
      jj = 1
      do 40 j=1,jjl1
      jj = jj+1
      wqj(jj,k,i,n) = 0.75e0*wq(j,k,i,n)+0.25e0*wq(j+1,k,i,n)
      jj = jj+1
      wqj(jj,k,i,n) = 0.25e0*wq(j,k,i,n)+0.75e0*wq(j+1,k,i,n)
   40 continue
c
c      interpolate in k
c
      ii = idim1-iil
      do 50 i=1,iil
      ii = ii+1
      do 1073 j=1,jdim1
      dq(j,1,ii,n) = wqj(j,1,i,n)
 1073 continue
c
cdir$ ivdep
      do 1005 izz=1,jdim1
      dq(izz,kdim1,ii,n) = wqj(izz,kkl,i,n)
 1005 continue
      kk = 1
      do 50 k=1,kkl1
      kk = kk+1
cdir$ ivdep
      do 1006 izz=1,jdim1
      dq(izz,kk,ii,n) = 0.75e0*wqj(izz,k,i,n)+0.25e0*wqj(izz,k+1,i,n)
 1006 continue
      kk = kk+1
cdir$ ivdep
      do 1007 izz=1,jdim1
      dq(izz,kk,ii,n) = 0.25e0*wqj(izz,k,i,n)+0.75e0*wqj(izz,k+1,i,n)
 1007 continue
   50 continue
c
c      interpolate in i
c
      if (idim1.gt.1) then
      is = idim1-iil+1
      np = jdim*kdim1-1
c
c     i=1
c
cdir$ ivdep
      do 1008 izz=1,np
      dq(izz,1,1,n) = dq(izz,1,is,n)
 1008 continue
      ii = 1
c
      do 60 i=is,idim2
      ii = ii+1
cdir$ ivdep
      do 1009 izz=1,np
      dq(izz,1,ii,n) = 0.75e0*dq(izz,1,i,n)+0.25e0*dq(izz,1,i+1,n)
 1009 continue
      ii = ii+1
cdir$ ivdep
      do 1010 izz=1,np
      dq(izz,1,ii,n) = 0.25e0*dq(izz,1,i,n)+0.75e0*dq(izz,1,i+1,n)
 1010 continue
   60 continue
c
      end if
   80 continue
c
      if (mode.eq.0) then
         nplq = min(idim1,999000/(jdim*kdim))
         npl  = nplq
         do 105 i=1,idim1,nplq
         if (i+npl-1.gt.idim1) npl = idim1+1-i
         nnpl = npl*jdim*kdim-jdim-1
         do 105 n=1,ll
cdir$ ivdep
         do 1012 izz=1,nnpl
         q(izz,1,i,n)=dq(izz,1,i,n)
 1012    continue
  105    continue
      else
         nplq = min(idim1,999000/(jdim*kdim))
         npl = nplq
         do 110 i=1,idim1,nplq
         if (i+npl-1.gt.idim1) npl = idim1+1-i
         nnpl = npl*jdim*kdim-jdim-1
c
c         update density and pressure to ensure positivity
c          - "cut-off" point is determined by alpq
c          - minimum value of density is equal to (1/phiq)
c
         alpq  = -.2
         phiq  = 1./0.5
         betq  = 1. + alpq*phiq
cdir$ ivdep
         do 7013 izz=1,nnpl
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
         do 1013 izz=1,nnpl
         q(izz,1,i,1) = q(izz,1,i,1)+blank(izz,1,i)*dq(izz,1,i,1)
         q(izz,1,i,2) = q(izz,1,i,2)+blank(izz,1,i)*dq(izz,1,i,2)
         q(izz,1,i,3) = q(izz,1,i,3)+blank(izz,1,i)*dq(izz,1,i,3)
         q(izz,1,i,4) = q(izz,1,i,4)+blank(izz,1,i)*dq(izz,1,i,4)
         q(izz,1,i,5) = q(izz,1,i,5)+blank(izz,1,i)*dq(izz,1,i,5)
 1013    continue
  110 continue
      end if
      call fill(jdim,kdim,idim,q,ll)
      return
      end
