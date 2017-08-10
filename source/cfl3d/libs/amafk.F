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
      subroutine amafk(i,npl,jdim,kdim,idim,q,ak,bk,ck,dtj,t,nvt,
     .                 dhp,dhm)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Formulate the implicit matrices in the K-direction for
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
      dimension q(jdim,kdim,idim,5),t(nvt,20),dtj(jdim,kdim,idim-1)
      dimension dhp((jdim-1)*npl,kdim,5,5),dhm((jdim-1)*npl,kdim,5,5)
      dimension ak(npl*(jdim-1),kdim-1,5,5),bk(npl*(jdim-1),kdim-1,5,5),
     .          ck(npl*(jdim-1),kdim-1,5,5)
c
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /precond/ cprec,uref,avn
c
c     matrix assembly - interior ponts
c
      kdim1 = kdim-1
      jdim1 = jdim-1
      n     = npl*jdim1*kdim1
      if (abs(ita).eq.1) then
        tfacp1=1.e0
      else
        tfacp1=1.5e0
      end if
      do 2750 k=1,5
      do 2750 l=1,5
cdir$ ivdep
      do 1000 izz=1,n
      bk(izz,1,k,l) = (dhp(izz,2,k,l)-dhm(izz,1,k,l))
      ak(izz,1,k,l) = -dhp(izz,1,k,l)
      ck(izz,1,k,l) =  dhm(izz,2,k,l)
 1000 continue
 2750 continue
c
c      assemble matrix equation - time terms
c
      if (real(cprec) .eq. 0.) then
         do 2730 ipl=1,npl
         ii  = i+ipl-1
         jkv = (ipl-1)*jdim1 + 1
         do 2730 k=1,kdim1
         jk1 = jkv + (k-1)*jdim1*npl
cdir$ ivdep
         do 1001 izz=1,jdim1
         t(izz+jk1-1,1) = q(izz,k,ii,1)
         t(izz+jk1-1,2) = q(izz,k,ii,2)
         t(izz+jk1-1,3) = q(izz,k,ii,3)
         t(izz+jk1-1,4) = q(izz,k,ii,4)
         t(izz+jk1-1,6) = tfacp1*dtj(izz,k,ii)
 1001    continue
 2730    continue
      else
         do 27301 ipl=1,npl
         ii  = i+ipl-1
         jkv = (ipl-1)*jdim1 + 1
         do 27301 k=1,kdim1
         jk1 = jkv + (k-1)*jdim1*npl
cdir$ ivdep
         do 10011 izz=1,jdim1
         t(izz+jk1-1,1) = q(izz,k,ii,1)
         t(izz+jk1-1,2) = q(izz,k,ii,2)
         t(izz+jk1-1,3) = q(izz,k,ii,3)
         t(izz+jk1-1,4) = q(izz,k,ii,4)
         t(izz+jk1-1,5) = q(izz,k,ii,5)
         t(izz+jk1-1,6) = tfacp1*dtj(izz,k,ii)
10011    continue
27301    continue
      end if
c
      if (real(cprec) .eq. 0.) then
cdir$ ivdep
      do 1002 izz=1,n
         temp          = t(izz,6)*t(izz,1)
         bk(izz,1,1,1) = bk(izz,1,1,1) + t(izz,6)
         bk(izz,1,2,1) = bk(izz,1,2,1) + t(izz,6)*t(izz,2)
         bk(izz,1,2,2) = bk(izz,1,2,2) + temp
         bk(izz,1,3,1) = bk(izz,1,3,1) + t(izz,6)*t(izz,3)
         bk(izz,1,3,3) = bk(izz,1,3,3) + temp
         bk(izz,1,4,1) = bk(izz,1,4,1) + t(izz,6)*t(izz,4)
         bk(izz,1,4,4) = bk(izz,1,4,4) + temp
         bk(izz,1,5,1) = bk(izz,1,5,1) 
     .                 + t(izz,6)*.5*(t(izz,2)*t(izz,2)+
     .                                t(izz,3)*t(izz,3)+
     .                                t(izz,4)*t(izz,4))
         bk(izz,1,5,2) = bk(izz,1,5,2) + temp*t(izz,2)
         bk(izz,1,5,3) = bk(izz,1,5,3) + temp*t(izz,3)
         bk(izz,1,5,4) = bk(izz,1,5,4) + temp*t(izz,4)
         bk(izz,1,5,5) = bk(izz,1,5,5) + t(izz,6)/gm1
 1002    continue
      else
cdir$ ivdep
         do 10021 izz=1,n
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
         bk(izz,1,1,1) = bk(izz,1,1,1) + t(izz,6)
         bk(izz,1,1,5) = bk(izz,1,1,5) + t(izz,6)*thet
         bk(izz,1,2,1) = bk(izz,1,2,1) + t(izz,6)*t(izz,2)
         bk(izz,1,2,2) = bk(izz,1,2,2) + temp
         bk(izz,1,2,5) = bk(izz,1,2,5) + t(izz,6)*thet*t(izz,2)
         bk(izz,1,3,1) = bk(izz,1,3,1) + t(izz,6)*t(izz,3)
         bk(izz,1,3,3) = bk(izz,1,3,3) + temp
         bk(izz,1,3,5) = bk(izz,1,3,5) + t(izz,6)*thet*t(izz,3)
         bk(izz,1,4,1) = bk(izz,1,4,1) + t(izz,6)*t(izz,4)
         bk(izz,1,4,4) = bk(izz,1,4,4) + temp
         bk(izz,1,4,5) = bk(izz,1,4,5) + t(izz,6)*thet*t(izz,4)
         bk(izz,1,5,1) = bk(izz,1,5,1) + t(izz,6)*ekin
         bk(izz,1,5,2) = bk(izz,1,5,2) + temp*t(izz,2)
         bk(izz,1,5,3) = bk(izz,1,5,3) + temp*t(izz,3)
         bk(izz,1,5,4) = bk(izz,1,5,4) + temp*t(izz,4)
         bk(izz,1,5,5) = bk(izz,1,5,5) + t(izz,6)*(1./gm1 + thet*ho)
10021    continue
      end if
      return
      end
