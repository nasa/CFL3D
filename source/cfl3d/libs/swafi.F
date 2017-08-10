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
      subroutine swafi(kp,npl,jdim,kdim,idim,q,ai,bi,ci,dtj,f,nvt,
     .                 res,imw)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Solve the block 5x5 tridiagonal equations for the 
c     3-factor spatially-split algorithm in the I-direction.
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
      dimension ai(npl*(jdim-1)/(imw+1),(idim-1)*(imw+1),5,5),
     .          bi(npl*(jdim-1)/(imw+1),(idim-1)*(imw+1),5,5),
     .          ci(npl*(jdim-1)/(imw+1),(idim-1)*(imw+1),5,5)
      dimension f(npl*(jdim-1)/(imw+1),(idim-1)*(imw+1),10)
c
      common /precond/ cprec,uref,avn
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
c
c     k-implicit j-sweep line inversions af
c
c     load rhs (-dq/jdt) into f
c
      idim1 = idim-1
      jdim1 = jdim-1
      gm1i  = 1.0/gm1
      if (abs(ita).eq.1) then
        tfacp1=1.e0
      else
        tfacp1=1.5e0
      end if
c
      if (imw.eq.0) then
         if (real(cprec) .eq. 0.) then
            do 1000 kpl=1,npl
            kk = kp+kpl-1
            js = (kpl-1)*jdim1 + 1
            do 1000 i=1,idim1
cdir$ ivdep
            do 1000 izz=1,jdim1
            temp            = tfacp1*dtj(izz,kk,i)*  q(izz,kk,i,1)
            f(izz+js-1,i,1) = tfacp1*dtj(izz,kk,i)*res(izz,kk,i,1)
            f(izz+js-1,i,5) = 0.5 *  f(izz+js-1,i,1)
     .                            * (q(izz,kk,i,2)*  q(izz,kk,i,2)
     .                              +q(izz,kk,i,3)*  q(izz,kk,i,3)
     .                              +q(izz,kk,i,4)*  q(izz,kk,i,4))
     .                       +temp*res(izz,kk,i,2)*  q(izz,kk,i,2)
     .                       +temp*res(izz,kk,i,3)*  q(izz,kk,i,3)
     .                       +temp*res(izz,kk,i,4)*  q(izz,kk,i,4)
     .                +tfacp1*gm1i*res(izz,kk,i,5)*dtj(izz,kk,i)
            f(izz+js-1,i,2) = temp*res(izz,kk,i,2)+
     .                               q(izz,kk,i,2)*  f(izz+js-1,i,1)
            f(izz+js-1,i,3) = temp*res(izz,kk,i,3)+
     .                               q(izz,kk,i,3)*  f(izz+js-1,i,1)
            f(izz+js-1,i,4) = temp*res(izz,kk,i,4)+
     .                               q(izz,kk,i,4)*  f(izz+js-1,i,1)
 1000       continue
         else
            do 10001 kpl=1,npl
            kk = kp+kpl-1
            js = (kpl-1)*jdim1 + 1
            do 10001 i=1,idim1
cdir$ ivdep
            do 10001 izz=1,jdim1
c
c           modifications for preconditioning
c
            c2 = gamma*q(izz,kk,i,5)/q(izz,kk,i,1)
            c = sqrt(c2)
            ekin = 0.5*(q(izz,kk,i,2)**2 + q(izz,kk,i,3)**2
     .                + q(izz,kk,i,4)**2)
            ho = c2/gm1 + ekin
            vmag1 = 2.0*ekin
            vel2 = ccmax(vmag1,avn*uref**2)
            vel = sqrt(ccmin(c2,vel2))
            vel = cprec*vel + (1.-cprec)*c
            thet = (1.0/vel**2 - 1.0/c2)
            restmp = tfacp1*dtj(izz,kk,i)*thet*res(izz,kk,i,5)
c
            temp            = tfacp1*dtj(izz,kk,i)*  q(izz,kk,i,1)
            f(izz+js-1,i,1) = tfacp1*dtj(izz,kk,i)*res(izz,kk,i,1)
     .                      + restmp
            f(izz+js-1,i,5) = f(izz+js-1,i,1)*ekin
     .                       +temp*res(izz,kk,i,2)*  q(izz,kk,i,2)
     .                       +temp*res(izz,kk,i,3)*  q(izz,kk,i,3)
     .                       +temp*res(izz,kk,i,4)*  q(izz,kk,i,4)
     .                +tfacp1*gm1i*res(izz,kk,i,5)*dtj(izz,kk,i)
     .                      + restmp*ho
            f(izz+js-1,i,2) = temp*res(izz,kk,i,2)+
     .                               q(izz,kk,i,2)*  f(izz+js-1,i,1)
     .                      + restmp*q(izz,kk,i,2)
            f(izz+js-1,i,3) = temp*res(izz,kk,i,3)+
     .                               q(izz,kk,i,3)*  f(izz+js-1,i,1)
     .                      + restmp*q(izz,kk,i,3)
            f(izz+js-1,i,4) = temp*res(izz,kk,i,4)+
     .                               q(izz,kk,i,4)*  f(izz+js-1,i,1)
     .                      + restmp*q(izz,kk,i,4)
10001       continue
         end if
      else
         if (real(cprec) .eq. 0.) then
            jdh = jdim1/2
            idh = idim1*2
            do 2411 kpl=1,npl
            kk  = kp+kpl-1
            js  = (kpl-1)*jdh + 1
            do 2411 i=1,idim1
            iq1 = i
            iq2 = idh+1-i
cdir$ ivdep
            do 1001 izz=1,jdh
            temp             =tfacp1*
     .                         dtj(izz+jdh,kk,i)*  q(izz+jdh,kk,i,1)
            f(izz+js-1,iq1,1)=tfacp1*
     .                         dtj(izz+jdh,kk,i)*res(izz+jdh,kk,i,1)
            f(izz+js-1,iq1,5)=0.5*f(izz+js-1,iq1,1)
     .                           *(q(izz+jdh,kk,i,2)*q(izz+jdh,kk,i,2)
     .                            +q(izz+jdh,kk,i,3)*q(izz+jdh,kk,i,3)
     .                            +q(izz+jdh,kk,i,4)*q(izz+jdh,kk,i,4))
     .                       +temp*res(izz+jdh,kk,i,2)*q(izz+jdh,kk,i,2)
     .                       +temp*res(izz+jdh,kk,i,3)*q(izz+jdh,kk,i,3)
     .                       +temp*res(izz+jdh,kk,i,4)*q(izz+jdh,kk,i,4)
     .                +tfacp1*gm1i*res(izz+jdh,kk,i,5)*dtj(izz+jdh,kk,i)
            f(izz+js-1,iq1,2)=temp*res(izz+jdh,kk,i,2)
     .                              +q(izz+jdh,kk,i,2)*f(izz+js-1,iq1,1)
            f(izz+js-1,iq1,3)=temp*res(izz+jdh,kk,i,3)
     .                              +q(izz+jdh,kk,i,3)*f(izz+js-1,iq1,1)
            f(izz+js-1,iq1,4)=temp*res(izz+jdh,kk,i,4)
     .                              +q(izz+jdh,kk,i,4)*f(izz+js-1,iq1,1)
 1001       continue
            call q8vrev(jdh,f(js,iq1,1),jdh,f(js,iq2,1))
            call q8vrev(jdh,f(js,iq1,2),jdh,f(js,iq2,2))
            call q8vrev(jdh,f(js,iq1,3),jdh,f(js,iq2,3))
            call q8vrev(jdh,f(js,iq1,4),jdh,f(js,iq2,4))
            call q8vrev(jdh,f(js,iq1,5),jdh,f(js,iq2,5))
cdir$ ivdep
            do 1002 izz=1,jdh
            temp              = tfacp1*dtj(izz,kk,i)*  q(izz,kk,i,1)
            f(izz+js-1,iq1,1) = tfacp1*dtj(izz,kk,i)*res(izz,kk,i,1)
            f(izz+js-1,iq1,5) = 0.5 *  f(izz+js-1,iq1,1)
     .                              * (q(izz,kk,i,2)*q(izz,kk,i,2)
     .                                +q(izz,kk,i,3)*q(izz,kk,i,3)
     .                                +q(izz,kk,i,4)*q(izz,kk,i,4))
     .                         +temp*res(izz,kk,i,2)*q(izz,kk,i,2)
     .                         +temp*res(izz,kk,i,3)*q(izz,kk,i,3)
     .                         +temp*res(izz,kk,i,4)*q(izz,kk,i,4)
     .                  +tfacp1*gm1i*res(izz,kk,i,5)*dtj(izz,kk,i)
            f(izz+js-1,iq1,2) = temp*res(izz,kk,i,2)
     .                                +q(izz,kk,i,2)*f(izz+js-1,iq1,1)
            f(izz+js-1,iq1,3) = temp*res(izz,kk,i,3)
     .                                +q(izz,kk,i,3)*f(izz+js-1,iq1,1)
            f(izz+js-1,iq1,4) = temp*res(izz,kk,i,4)
     .                                +q(izz,kk,i,4)*f(izz+js-1,iq1,1)
 1002       continue
 2411       continue
         else
            jdh = jdim1/2
            idh = idim1*2
            do 24111 kpl=1,npl
            kk  = kp+kpl-1
            js  = (kpl-1)*jdh + 1
            do 24111 i=1,idim1
            iq1 = i
            iq2 = idh+1-i
cdir$ ivdep
            do 10011 izz=1,jdh
c
c           modifications for preconditioning
c 
            c2 = gamma*q(izz+jdh,kk,i,5)/q(izz+jdh,kk,i,1)
            c = sqrt(c2)
            ekin = 0.5*(q(izz+jdh,kk,i,2)**2 + q(izz+jdh,kk,i,3)**2 
     .                + q(izz+jdh,kk,i,4)**2)
            ho = c2/gm1 + ekin 
            vmag1 = 2.0*ekin
            vel2 = ccmax(vmag1,avn*uref**2)
            vel = sqrt(ccmin(c2,vel2))
            vel = cprec*vel + (1.-cprec)*c
            thet = (1.0/vel**2 - 1.0/c2)
            restmp = tfacp1*dtj(izz+jdh,kk,i)*thet*res(izz+jdh,kk,i,5)
c
            temp             =tfacp1*
     .                          dtj(izz+jdh,kk,i)*  q(izz+jdh,kk,i,1)
            f(izz+js-1,iq1,1)=tfacp1*
     .                        dtj(izz+jdh,kk,i)*res(izz+jdh,kk,i,1)
     .                       +restmp
            f(izz+js-1,iq1,5)=f(izz+js-1,iq1,1)*ekin
     .                       +temp*res(izz+jdh,kk,i,2)*q(izz+jdh,kk,i,2)
     .                       +temp*res(izz+jdh,kk,i,3)*q(izz+jdh,kk,i,3)
     .                       +temp*res(izz+jdh,kk,i,4)*q(izz+jdh,kk,i,4)
     .                +tfacp1*gm1i*res(izz+jdh,kk,i,5)*dtj(izz+jdh,kk,i)
     .                       + restmp*ho
            f(izz+js-1,iq1,2)=temp*res(izz+jdh,kk,i,2)
     .                              +q(izz+jdh,kk,i,2)*f(izz+js-1,iq1,1)
     .                       +restmp*q(izz+jdh,kk,i,2)
            f(izz+js-1,iq1,3)=temp*res(izz+jdh,kk,i,3)
     .                              +q(izz+jdh,kk,i,3)*f(izz+js-1,iq1,1)
     .                       +restmp*q(izz+jdh,kk,i,3)
            f(izz+js-1,iq1,4)=temp*res(izz+jdh,kk,i,4)
     .                              +q(izz+jdh,kk,i,4)*f(izz+js-1,iq1,1)
     .                       +restmp*q(izz+jdh,kk,i,4)
10011       continue
            call q8vrev(jdh,f(js,iq1,1),jdh,f(js,iq2,1))
            call q8vrev(jdh,f(js,iq1,2),jdh,f(js,iq2,2))
            call q8vrev(jdh,f(js,iq1,3),jdh,f(js,iq2,3))
            call q8vrev(jdh,f(js,iq1,4),jdh,f(js,iq2,4))
            call q8vrev(jdh,f(js,iq1,5),jdh,f(js,iq2,5))
cdir$ ivdep
            do 10021 izz=1,jdh
c
c           modifications for preconditioning
c
            c2 = gamma*q(izz,kk,i,5)/q(izz,kk,i,1)
            c = sqrt(c2)
            ekin = 0.5*(q(izz,kk,i,2)**2 + q(izz,kk,i,3)**2 
     .                + q(izz,kk,i,4)**2)
            ho = c2/gm1 + ekin 
            vmag1 = 2.0*ekin
            vel2 = ccmax(vmag1,avn*uref**2)
            vel = sqrt(ccmin(c2,vel2))
            vel = cprec*vel + (1.-cprec)*c
            thet = (1.0/vel**2 - 1.0/c2)
            restmp = tfacp1*dtj(izz,kk,i)*thet*res(izz,kk,i,5)
c
            temp              = tfacp1*dtj(izz,kk,i)*  q(izz,kk,i,1)
            f(izz+js-1,iq1,1) = tfacp1*dtj(izz,kk,i)*res(izz,kk,i,1)
     .                        + restmp  
            f(izz+js-1,iq1,5) = f(izz+js-1,iq1,1)*ekin
     .                         +temp*res(izz,kk,i,2)*q(izz,kk,i,2)
     .                         +temp*res(izz,kk,i,3)*q(izz,kk,i,3)
     .                         +temp*res(izz,kk,i,4)*q(izz,kk,i,4)
     .                  +tfacp1*gm1i*res(izz,kk,i,5)*dtj(izz,kk,i)
     .                        + restmp*ho  
            f(izz+js-1,iq1,2) = temp*res(izz,kk,i,2)
     .                                +q(izz,kk,i,2)*f(izz+js-1,iq1,1)
     .                        + restmp*q(izz,kk,i,2)
            f(izz+js-1,iq1,3) = temp*res(izz,kk,i,3)
     .                                +q(izz,kk,i,3)*f(izz+js-1,iq1,1)
     .                        + restmp*q(izz,kk,i,3)
            f(izz+js-1,iq1,4) = temp*res(izz,kk,i,4)
     .                                +q(izz,kk,i,4)*f(izz+js-1,iq1,1)
     .                        + restmp*q(izz,kk,i,4)
10021       continue
24111       continue
         end if
      end if
c
c     solve matrix equation
c
      il  = 1
      iu  = idim1*(imw+1)
      n   = npl*jdim1/(imw+1)
c
      id1 = npl*(jdim-1)/(imw+1)
      id2 = (idim-1)*(imw+1)
      call bsub(id1,id2,ai,bi,ci,f,1,n,il,iu)
c
c     update delta q
c
      if (imw.eq.0) then
         do 2820 kpl=1,npl
         kk = kp+kpl-1
         jv = (kpl-1)*jdim1 + 1
         do 2820 i=1,idim1
         do 2820 l=1,5
cdir$ ivdep
         do 1007 izz=1,jdim1
         res(izz,kk,i,l) = f(izz+jv-1,i,l)
 1007    continue
 2820    continue
      else
         do 2821 kpl=1,npl
         kk  = kp+kpl-1
         jv  = (kpl-1)*jdh + 1
         do 2821 i=1,idim1
         iq1 = i
         iq2 = idh+1-i
         do 2821 l=1,5
cdir$ ivdep
         do 1008 izz=1,jdh
         res(izz,kk,i,l) = f(izz+jv-1,iq1,l)
 1008    continue
         call q8vrev(jdh,f(jv,iq2,l),jdh,res(jdh+1,kk,i,l))
 2821    continue
      end if
 3000 continue
      return
      end
