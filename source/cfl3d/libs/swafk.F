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
      subroutine swafk(i,npl,jdim,kdim,idim,q,ak,bk,ck,dtj,f,nvt,
     .                 res,imw)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Solve the block 5x5 tridiagonal equations for the 
c     3-factor spatially-split algorithm in the K-direction.
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
      dimension res(jdim,kdim,idim-1,5),dtj(jdim,kdim,idim-1)
      dimension ak(npl*(jdim-1)/(imw+1),(kdim-1)*(imw+1),5,5),
     .          bk(npl*(jdim-1)/(imw+1),(kdim-1)*(imw+1),5,5),
     .          ck(npl*(jdim-1)/(imw+1),(kdim-1)*(imw+1),5,5)
      dimension f(npl*(jdim-1)/(imw+1),(kdim-1)*(imw+1),5)
c
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /precond/ cprec,uref,avn
c
c     k-implicit j-sweep line inversions af
c
c     load rhs (-dq/jdt) into f
c
      kdim1 = kdim-1
      jdim1 = jdim-1
      n     = jdim*kdim1-1
      gm1i  = 1.0/gm1
      if (abs(ita).eq.1) then
        tfacp1=1.e0
      else
        tfacp1=1.5e0
      end if
c
      if (imw.eq.0) then
         if (real(cprec) .eq. 0.) then
            do 2410 ipl=1,npl
            ii = i+ipl-1
cdir$ ivdep
            do 1000 izz=1,n
            temp            = tfacp1*dtj(izz,1,ii)*  q(izz,1,ii,1)
            res(izz,1,ii,1) = tfacp1*dtj(izz,1,ii)*res(izz,1,ii,1)
            res(izz,1,ii,5) = 0.5* res(izz,1,ii,1)
     .                             *(q(izz,1,ii,2)*  q(izz,1,ii,2)
     .                              +q(izz,1,ii,3)*  q(izz,1,ii,3)
     .                              +q(izz,1,ii,4)*  q(izz,1,ii,4))
     .                       +temp*res(izz,1,ii,2)*  q(izz,1,ii,2)
     .                       +temp*res(izz,1,ii,3)*  q(izz,1,ii,3)
     .                       +temp*res(izz,1,ii,4)*  q(izz,1,ii,4)
     .                +tfacp1*gm1i*res(izz,1,ii,5)*dtj(izz,1,ii)
            res(izz,1,ii,2) = temp*res(izz,1,ii,2)+
     .                               q(izz,1,ii,2)*res(izz,1,ii,1)
            res(izz,1,ii,3) = temp*res(izz,1,ii,3)+
     .                               q(izz,1,ii,3)*res(izz,1,ii,1)
            res(izz,1,ii,4) = temp*res(izz,1,ii,4)+
     .                               q(izz,1,ii,4)*res(izz,1,ii,1)
 1000       continue
            do 2410 l=1,5
            js = (ipl-1)*jdim1 + 1
            do 2410 k=1,kdim1
            jk = (k-1)*jdim + 1
cdir$ ivdep
            do 1001 izz=1,jdim1
            f(izz+js-1,k,l) = res(izz+jk-1,1,ii,l)
 1001       continue
 2410       continue
         else
            do 24101 ipl=1,npl
            ii = i+ipl-1
cdir$ ivdep
            do 10001 izz=1,n
c
c           modifications for preconditioning
c 
            c2 = gamma*q(izz,1,ii,5)/q(izz,1,ii,1)
            c = sqrt(c2)
            ekin = 0.5*(q(izz,1,ii,2)**2 + q(izz,1,ii,3)**2 
     .                + q(izz,1,ii,4)**2)
            ho = c2/gm1 + ekin 
            vmag1 = 2.0*ekin
            vel2 = ccmax(vmag1,avn*uref**2)
            vel = sqrt(ccmin(c2,vel2))
            vel = cprec*vel + (1.-cprec)*c
            thet = (1.0/vel**2 - 1.0/c2)
            restmp = tfacp1*dtj(izz,1,ii)*thet*res(izz,1,ii,5)
c
            temp            = tfacp1*dtj(izz,1,ii)*  q(izz,1,ii,1)
            res(izz,1,ii,1) = tfacp1*dtj(izz,1,ii)*res(izz,1,ii,1)
     .                      + restmp
            res(izz,1,ii,5) = res(izz,1,ii,1)*ekin
     .                       +temp*res(izz,1,ii,2)*  q(izz,1,ii,2)
     .                       +temp*res(izz,1,ii,3)*  q(izz,1,ii,3)
     .                       +temp*res(izz,1,ii,4)*  q(izz,1,ii,4)
     .                +tfacp1*gm1i*res(izz,1,ii,5)*dtj(izz,1,ii)
     .                      + restmp*ho
            res(izz,1,ii,2) = temp*res(izz,1,ii,2)+
     .                               q(izz,1,ii,2)*res(izz,1,ii,1)
     .                      + restmp*q(izz,1,ii,2)
            res(izz,1,ii,3) = temp*res(izz,1,ii,3)+
     .                               q(izz,1,ii,3)*res(izz,1,ii,1)
     .                      + restmp*q(izz,1,ii,3)
            res(izz,1,ii,4) = temp*res(izz,1,ii,4)+
     .                               q(izz,1,ii,4)*res(izz,1,ii,1)
     .                      + restmp*q(izz,1,ii,4)
10001       continue
            do 24101 l=1,5
            js = (ipl-1)*jdim1 + 1
            do 24101 k=1,kdim1
            jk = (k-1)*jdim + 1
cdir$ ivdep
            do 10011 izz=1,jdim1
            f(izz+js-1,k,l) = res(izz+jk-1,1,ii,l)
10011       continue
24101       continue
         end if
      else
         if (real(cprec) .eq. 0.) then
            jdh = jdim1/2
            do 2710 ipl=1,npl
            ii  = i+ipl-1
cdir$ ivdep
            do 1002 izz=1,n
            temp            = tfacp1*dtj(izz,1,ii)*  q(izz,1,ii,1)
            res(izz,1,ii,1) = tfacp1*dtj(izz,1,ii)*res(izz,1,ii,1)
            res(izz,1,ii,5) = 0.5 *res(izz,1,ii,1)
     .                            * (q(izz,1,ii,2)*  q(izz,1,ii,2)
     .                              +q(izz,1,ii,3)*  q(izz,1,ii,3)
     .                              +q(izz,1,ii,4)*  q(izz,1,ii,4))
     .                       +temp*res(izz,1,ii,2)*  q(izz,1,ii,2)
     .                       +temp*res(izz,1,ii,3)*  q(izz,1,ii,3)
     .                       +temp*res(izz,1,ii,4)*  q(izz,1,ii,4)
     .                +tfacp1*gm1i*res(izz,1,ii,5)*dtj(izz,1,ii)
            res(izz,1,ii,2) = temp*res(izz,1,ii,2)+
     .                               q(izz,1,ii,2)*res(izz,1,ii,1)
            res(izz,1,ii,3) = temp*res(izz,1,ii,3)+
     .                               q(izz,1,ii,3)*res(izz,1,ii,1)
            res(izz,1,ii,4) = temp*res(izz,1,ii,4)+
     .                               q(izz,1,ii,4)*res(izz,1,ii,1)
 1002       continue
            do 2710 l=1,5
            do 2710 k=1,kdim1
            kq1 = kdim1 + k
            jv  = (ipl-1)*jdh + 1
            kq2 = kdim - k
            jk  = (k-1)*jdim + 1
cdir$ ivdep
            do 1003 izz=1,jdh
            f(izz+jv-1,kq1,l) = res(izz+jk-1,1,ii,l)
 1003       continue
            call q8vrev(jdh,res(jk+jdh,1,ii,l),jdh,f(jv,kq2,l))
 2710       continue
         else
            jdh = jdim1/2
            do 27101 ipl=1,npl
            ii  = i+ipl-1
cdir$ ivdep
            do 10021 izz=1,n
c
c           modifications for preconditioning
c
            c2 = gamma*q(izz,1,ii,5)/q(izz,1,ii,1)
            c = sqrt(c2)
            ekin = 0.5*(q(izz,1,ii,2)**2 + q(izz,1,ii,3)**2 
     .                + q(izz,1,ii,4)**2)
            ho = c2/gm1 + ekin 
            vmag1 = 2.0*ekin
            vel2 = ccmax(vmag1,avn*uref**2)
            vel = sqrt(ccmin(c2,vel2))
            vel = cprec*vel + (1.-cprec)*c
            thet = (1.0/vel**2 - 1.0/c2)
            restmp = tfacp1*dtj(izz,1,ii)*thet*res(izz,1,ii,5)
c
            temp            = tfacp1*dtj(izz,1,ii)*  q(izz,1,ii,1)
            res(izz,1,ii,1) = tfacp1*dtj(izz,1,ii)*res(izz,1,ii,1)
     .                      + restmp
            res(izz,1,ii,5) = 0.5 *res(izz,1,ii,1)
     .                            * (q(izz,1,ii,2)*  q(izz,1,ii,2)
     .                              +q(izz,1,ii,3)*  q(izz,1,ii,3)
     .                              +q(izz,1,ii,4)*  q(izz,1,ii,4))
     .                       +temp*res(izz,1,ii,2)*  q(izz,1,ii,2)
     .                       +temp*res(izz,1,ii,3)*  q(izz,1,ii,3)
     .                       +temp*res(izz,1,ii,4)*  q(izz,1,ii,4)
     .                +tfacp1*gm1i*res(izz,1,ii,5)*dtj(izz,1,ii)
     .                      + restmp*ho
            res(izz,1,ii,2) = temp*res(izz,1,ii,2)+
     .                               q(izz,1,ii,2)*res(izz,1,ii,1)
     .                      + restmp*q(izz,1,ii,2)
            res(izz,1,ii,3) = temp*res(izz,1,ii,3)+
     .                               q(izz,1,ii,3)*res(izz,1,ii,1)
     .                      + restmp*q(izz,1,ii,3)
            res(izz,1,ii,4) = temp*res(izz,1,ii,4)+
     .                               q(izz,1,ii,4)*res(izz,1,ii,1)
     .                      + restmp*q(izz,1,ii,4)
10021       continue
            do 27101 l=1,5
            do 27101 k=1,kdim1
            kq1 = kdim1 + k
            jv  = (ipl-1)*jdh + 1
            kq2 = kdim - k
            jk  = (k-1)*jdim + 1
cdir$ ivdep
            do 10031 izz=1,jdh
            f(izz+jv-1,kq1,l) = res(izz+jk-1,1,ii,l)
10031       continue
            call q8vrev(jdh,res(jk+jdh,1,ii,l),jdh,f(jv,kq2,l))
27101       continue
         end if
      end if
c
c     solve matrix equation
c
      il  = 1
      iu  = kdim1*(imw+1)
      n   = npl*jdim1/(imw+1)
c
      id1 = npl*(jdim-1)/(imw+1)
      id2 = (kdim-1)*(imw+1)
      call bsub(id1,id2,ak,bk,ck,f,1,n,il,iu)
c
c     update delta q
c
      if (imw.eq.0) then
         do 2820 ipl=1,npl
         ii = i+ipl-1
         jv = (ipl-1)*jdim1 + 1
         do 2820 k=1,kdim1
         do 2820 l=1,5
cdir$ ivdep
         do 1008 izz=1,jdim1
         res(izz,k,ii,l) = f(izz+jv-1,k,l)
 1008    continue
 2820    continue
      else
         do 2825 ipl=1,npl
         ii  = i+ipl-1
         jv  = (ipl-1)*jdh + 1
         do 2825 k=1,kdim1
         kq1 = kdim1+k
         kq2 = kdim-k
         do 2825 l=1,5
cdir$ ivdep
         do 1009 izz=1,jdh
         res(izz,k,ii,l) = f(izz+jv-1,kq1,l)
 1009    continue
         call q8vrev(jdh,f(jv,kq2,l),jdh,res(1+jdh,k,ii,l))
 2825    continue
      end if
 3000 continue
      return
      end
