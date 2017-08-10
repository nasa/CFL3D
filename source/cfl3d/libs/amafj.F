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
      subroutine amafj(i,npl,jdim,kdim,idim,q,aj,bj,cj,dtj,t,nvt,
     .                 dgp,dgm)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Formulate the implicit matrices in the J-direction for 
c     the 3-factor algorithm.
c     Modified for Weiss-Smith preconditioning by J.R. Edwards, NCSU
c       cprec = 0 ---> original code used
c             > 0 ---> modified code used
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension q(jdim,kdim,idim,5)
      dimension dgp(jdim,npl*(kdim-1),5,5),dgm(jdim,npl*(kdim-1),5,5)
      dimension t(nvt,20),dtj(jdim,kdim,idim-1)
      dimension aj(npl*(kdim-1),jdim,5,5),bj(npl*(kdim-1),jdim,5,5),
     .          cj(npl*(kdim-1),jdim,5,5)
c
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /precond/ cprec,uref,avn
c
c     assemble matrix equation - interior points
c
      kdim1 = kdim-1
      jdim1 = jdim-1
      kv    = npl*kdim1
      n     = jdim1*kv
      if (abs(ita).eq.1) then
        tfacp1=1.e0
      else
        tfacp1=1.5e0
      end if
c
      do 1232 k=1,5
      do 1232 l=1,5
      do 1230 j=1,jdim
      n0 = (j-1)*kv+1
c
      jj = 1-jdim
      do 8893 kk=1,kv
      jj = jj+jdim
      t(n0+kk-1,1) = dgp(j+jj-1,1,k,l)
 8893 t(n0+kk-1,2) = dgm(j+jj-1,1,k,l)
c      call q8vgathp(kv,dgp(j,1,k,l),jdim,kv,kv,t(n0,1))
c      call q8vgathp(kv,dgm(j,1,k,l),jdim,kv,kv,t(n0,2))
 1230 continue
cdir$ ivdep
      do 1000 izz=1,n
      bj(izz,1,k,l) = (t(izz+kv,1)-t(izz,2))
      aj(izz,1,k,l) = -t(izz,1)
      cj(izz,1,k,l) =  t(izz+kv,2)
 1000 continue
 1232 continue
c
c      assemble matrix equation - time terms
c
      if (real(cprec) .eq. 0.) then
         do 1216 ipl=1,npl
         ii  = i+ipl-1
         jkv = (ipl-1)*kdim1
         do 1216 j=1,jdim1
         n0  = (j-1)*kv + jkv + 1
c
         jj  = 1-jdim
cdir$ ivdep
         do 7886 kk=1,kdim1
         jj  = jj+jdim
         t(n0+kk-1,1) = q(j+jj-1,1,ii,1)
         t(n0+kk-1,2) = q(j+jj-1,1,ii,2)
         t(n0+kk-1,3) = q(j+jj-1,1,ii,3)
         t(n0+kk-1,4) = q(j+jj-1,1,ii,4)
         t(n0+kk-1,6) = tfacp1*dtj(j+jj-1,1,ii)
 7886    continue
c        call q8vgathp(kdim1,dtj(j,1,ii),jdim,kdim1,kdim1,t(n0,1))
 1216    continue
      else
         do 12161 ipl=1,npl
         ii  = i+ipl-1
         jkv = (ipl-1)*kdim1
         do 12161 j=1,jdim1
         n0  = (j-1)*kv + jkv + 1
c
         jj  = 1-jdim
cdir$ ivdep
         do 78861 kk=1,kdim1
         jj  = jj+jdim
         t(n0+kk-1,1) = q(j+jj-1,1,ii,1)
         t(n0+kk-1,2) = q(j+jj-1,1,ii,2)
         t(n0+kk-1,3) = q(j+jj-1,1,ii,3)
         t(n0+kk-1,4) = q(j+jj-1,1,ii,4)
         t(n0+kk-1,5) = q(j+jj-1,1,ii,5)
         t(n0+kk-1,6) = tfacp1*dtj(j+jj-1,1,ii)
78861    continue
c        call q8vgathp(kdim1,dtj(j,1,ii),jdim,kdim1,kdim1,t(n0,1))
12161    continue
      end if
c
      if (real(cprec) .eq. 0.) then
cdir$ ivdep
         do 1001 izz=1,n
         temp          = t(izz,6)*t(izz,1)
         bj(izz,1,1,1) = bj(izz,1,1,1) + t(izz,6)
         bj(izz,1,2,1) = bj(izz,1,2,1) + t(izz,6)*t(izz,2)
         bj(izz,1,2,2) = bj(izz,1,2,2) + temp
         bj(izz,1,3,1) = bj(izz,1,3,1) + t(izz,6)*t(izz,3)
         bj(izz,1,3,3) = bj(izz,1,3,3) + temp
         bj(izz,1,4,1) = bj(izz,1,4,1) + t(izz,6)*t(izz,4)
         bj(izz,1,4,4) = bj(izz,1,4,4) + temp
         bj(izz,1,5,1) = bj(izz,1,5,1) 
     .                 + t(izz,6)*0.5*(t(izz,2)*t(izz,2)+
     .                                 t(izz,3)*t(izz,3)+
     .                                 t(izz,4)*t(izz,4))
         bj(izz,1,5,2) = bj(izz,1,5,2) + temp*t(izz,2)
         bj(izz,1,5,3) = bj(izz,1,5,3) + temp*t(izz,3)
         bj(izz,1,5,4) = bj(izz,1,5,4) + temp*t(izz,4)
         bj(izz,1,5,5) = bj(izz,1,5,5) + t(izz,6)/gm1
 1001    continue
      else
cdir$ ivdep
         do 10011 izz=1,n
         c2 = gamma*t(izz,5)/t(izz,1)
         c = sqrt(c2)
         ekin = 0.5*(t(izz,2)**2 + t(izz,3)**2 + t(izz,4)**2)
         ho = c2/gm1 + ekin
         vmag1 = 2.0*ekin
         vel2 = ccmax(vmag1,avn*uref**2)
         vel = sqrt(ccmin(c2,vel2))
         vel = cprec*vel + (1.-cprec)*c
         thet = (1.0/vel**2 - 1.0/c2)
         temp          = t(izz,6)*t(izz,1)
         bj(izz,1,1,1) = bj(izz,1,1,1) + t(izz,6)
         bj(izz,1,1,5) = bj(izz,1,1,5) + t(izz,6)*thet
         bj(izz,1,2,1) = bj(izz,1,2,1) + t(izz,6)*t(izz,2)
         bj(izz,1,2,2) = bj(izz,1,2,2) + temp
         bj(izz,1,2,5) = bj(izz,1,2,5) + t(izz,6)*thet*t(izz,2)
         bj(izz,1,3,1) = bj(izz,1,3,1) + t(izz,6)*t(izz,3)
         bj(izz,1,3,3) = bj(izz,1,3,3) + temp
         bj(izz,1,3,5) = bj(izz,1,3,5) + t(izz,6)*thet*t(izz,3)
         bj(izz,1,4,1) = bj(izz,1,4,1) + t(izz,6)*t(izz,4)
         bj(izz,1,4,4) = bj(izz,1,4,4) + temp
         bj(izz,1,4,5) = bj(izz,1,4,5) + t(izz,6)*thet*t(izz,4)
         bj(izz,1,5,1) = bj(izz,1,5,1) + t(izz,6)*ekin
         bj(izz,1,5,2) = bj(izz,1,5,2) + temp*t(izz,2)
         bj(izz,1,5,3) = bj(izz,1,5,3) + temp*t(izz,3)
         bj(izz,1,5,4) = bj(izz,1,5,4) + temp*t(izz,4)
         bj(izz,1,5,5) = bj(izz,1,5,5) + t(izz,6)*(1.0/gm1 + thet*ho)
10011   continue
      end if
c
      return
      end
