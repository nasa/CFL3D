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
      subroutine amafi(kp,npl,jdim,kdim,idim,q,ai,bi,ci,dtj,t,nvt,
     .                 dfp,dfm,imw)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Formulate the implicit matrices in the I-direction for
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
      dimension dfp(npl*(jdim-1),idim,5,5),dfm(npl*(jdim-1),idim,5,5)
      dimension t(nvt,20),dtj(jdim,kdim,idim-1)
      dimension ai(npl*(jdim-1)/(imw+1),(idim-1)*(imw+1),5,5),
     .          bi(npl*(jdim-1)/(imw+1),(idim-1)*(imw+1),5,5),
     .          ci(npl*(jdim-1)/(imw+1),(idim-1)*(imw+1),5,5)
c
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /precond/ cprec,uref,avn
c
c     matrix assembly - interior ponts
c
      idim1 = idim-1
      jdim1 = jdim-1
      n     = npl*jdim1*idim1
      if (abs(ita).eq.1) then
        tfacp1=1.e0
      else
        tfacp1=1.5e0
      end if
c
      if (imw.eq.0) then
c
         do 2750 m=1,5
         do 2750 l=1,5
cdir$ ivdep
         do 1000 izz=1,n
         bi(izz,1,m,l) = (dfp(izz,2,m,l)-dfm(izz,1,m,l))
         ai(izz,1,m,l) = -dfp(izz,1,m,l)
         ci(izz,1,m,l) =  dfm(izz,2,m,l)
 1000    continue
 2750    continue
c
c      assemble matrix equation - time terms
c
         if (real(cprec) .eq. 0.) then
            do 2730 kpl=1,npl
            kk  = kp+kpl-1
            jiv = (kpl-1)*jdim1 + 1
            do 2730 i=1,idim1
            ji1 = jiv + (i-1)*jdim1*npl
cdir$ ivdep
            do 1001 izz=1,jdim1
            t(izz+ji1-1,1) = q(izz,kk,i,1)
            t(izz+ji1-1,2) = q(izz,kk,i,2)
            t(izz+ji1-1,3) = q(izz,kk,i,3)
            t(izz+ji1-1,4) = q(izz,kk,i,4)
            t(izz+ji1-1,6) = tfacp1*dtj(izz,kk,i)
 1001       continue
 2730       continue
         else
            do 27301 kpl=1,npl
            kk  = kp+kpl-1
            jiv = (kpl-1)*jdim1 + 1
            do 27301 i=1,idim1
            ji1 = jiv + (i-1)*jdim1*npl
cdir$ ivdep
            do 10011 izz=1,jdim1
            t(izz+ji1-1,1) = q(izz,kk,i,1)
            t(izz+ji1-1,2) = q(izz,kk,i,2)
            t(izz+ji1-1,3) = q(izz,kk,i,3)
            t(izz+ji1-1,4) = q(izz,kk,i,4)
            t(izz+ji1-1,5) = q(izz,kk,i,5)
            t(izz+ji1-1,6) = tfacp1*dtj(izz,kk,i)
10011       continue
27301       continue
         end if
c
         if (real(cprec) .eq. 0.) then
cdir$ ivdep
            do 1002 izz=1,n
            temp          = t(izz,6)*t(izz,1)
            bi(izz,1,1,1) = bi(izz,1,1,1) + t(izz,6)
            bi(izz,1,2,1) = bi(izz,1,2,1) + t(izz,6)*t(izz,2)
            bi(izz,1,2,2) = bi(izz,1,2,2) + temp
            bi(izz,1,3,1) = bi(izz,1,3,1) + t(izz,6)*t(izz,3)
            bi(izz,1,3,3) = bi(izz,1,3,3) + temp
            bi(izz,1,4,1) = bi(izz,1,4,1) + t(izz,6)*t(izz,4)
            bi(izz,1,4,4) = bi(izz,1,4,4) + temp
            bi(izz,1,5,1) = bi(izz,1,5,1) + 
     .                      t(izz,6)*.5*(t(izz,2)*t(izz,2)+
     .                                   t(izz,3)*t(izz,3)+
     .                                   t(izz,4)*t(izz,4))
            bi(izz,1,5,2) = bi(izz,1,5,2) + temp*t(izz,2)
            bi(izz,1,5,3) = bi(izz,1,5,3) + temp*t(izz,3)
            bi(izz,1,5,4) = bi(izz,1,5,4) + temp*t(izz,4)
            bi(izz,1,5,5) = bi(izz,1,5,5) + t(izz,6)/gm1
 1002       continue
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
            bi(izz,1,1,1) = bi(izz,1,1,1) + t(izz,6)
            bi(izz,1,1,5) = bi(izz,1,1,5) + t(izz,6)*thet
            bi(izz,1,2,1) = bi(izz,1,2,1) + t(izz,6)*t(izz,2)
            bi(izz,1,2,2) = bi(izz,1,2,2) + temp
            bi(izz,1,2,5) = bi(izz,1,2,5) + t(izz,6)*thet*t(izz,2)
            bi(izz,1,3,1) = bi(izz,1,3,1) + t(izz,6)*t(izz,3)
            bi(izz,1,3,3) = bi(izz,1,3,3) + temp
            bi(izz,1,3,5) = bi(izz,1,3,5) + t(izz,6)*thet*t(izz,3)
            bi(izz,1,4,1) = bi(izz,1,4,1) + t(izz,6)*t(izz,4)
            bi(izz,1,4,4) = bi(izz,1,4,4) + temp
            bi(izz,1,4,5) = bi(izz,1,4,5) + t(izz,6)*thet*t(izz,4)
            bi(izz,1,5,1) = bi(izz,1,5,1) + t(izz,6)*ekin
            bi(izz,1,5,2) = bi(izz,1,5,2) + temp*t(izz,2)
            bi(izz,1,5,3) = bi(izz,1,5,3) + temp*t(izz,3)
            bi(izz,1,5,4) = bi(izz,1,5,4) + temp*t(izz,4)
            bi(izz,1,5,5) = bi(izz,1,5,5) + t(izz,6)*(1./gm1 + thet*ho)
10021       continue
         end if
c
      else
c
         jdh = jdim1/2
         mm  = jdh
         idh = idim1*2
c
         if (real(cprec) .eq. 0.) then
            do 4000 kpl=1,npl
            kk  = kp+kpl-1
            jv1 = (kpl-1)*jdh + 1
            jv2 = (kpl-1)*jdim1 + 1
            do 4000 i=1,idim1
            iq  = i
            do 3500 k=1,5
            do 3500 l=1,5
cdir$ ivdep
            do 1003 izz=1,mm
            bi(izz+jv1-1,iq,k,l) = (dfp(izz+jv2-1,i+1,k,l)
     .                             -dfm(izz+jv2-1,i,k,l))
            ai(izz+jv1-1,iq,k,l) = -dfp(izz+jv2-1,i,k,l)
            ci(izz+jv1-1,iq,k,l) =  dfm(izz+jv2-1,i+1,k,l)
 1003       continue
 3500       continue
            jv3 = jv1+(iq-1)*jdh*npl
cdir$ ivdep
            do 1004 izz=1,mm
            t(izz+jv3-1,1) = q(izz,kk,i,1)
            t(izz+jv3-1,2) = q(izz,kk,i,2)
            t(izz+jv3-1,3) = q(izz,kk,i,3)
            t(izz+jv3-1,4) = q(izz,kk,i,4)
            t(izz+jv3-1,6) = dtj(izz,kk,i)
 1004       continue
 4000       continue
         else
            do 40001 kpl=1,npl
            kk  = kp+kpl-1
            jv1 = (kpl-1)*jdh + 1
            jv2 = (kpl-1)*jdim1 + 1
            do 40001 i=1,idim1
            iq  = i
            do 35001 k=1,5
            do 35001 l=1,5
cdir$ ivdep
            do 10031 izz=1,mm
            bi(izz+jv1-1,iq,k,l) = (dfp(izz+jv2-1,i+1,k,l)
     .                             -dfm(izz+jv2-1,i,k,l))
            ai(izz+jv1-1,iq,k,l) = -dfp(izz+jv2-1,i,k,l)
            ci(izz+jv1-1,iq,k,l) =  dfm(izz+jv2-1,i+1,k,l)
10031       continue
35001       continue
            jv3 = jv1+(iq-1)*jdh*npl
cdir$ ivdep
            do 10041 izz=1,mm
            t(izz+jv3-1,1) = q(izz,kk,i,1)
            t(izz+jv3-1,2) = q(izz,kk,i,2)
            t(izz+jv3-1,3) = q(izz,kk,i,3)
            t(izz+jv3-1,4) = q(izz,kk,i,4)
            t(izz+jv3-1,5) = q(izz,kk,i,5)
            t(izz+jv3-1,6) = dtj(izz,kk,i)
10041       continue
40001       continue
         end if
c
         do 6000 kpl=1,npl
         kk  = kp+kpl-1
         jv1 = (kpl-1)*jdh + 1
         jv2 = (kpl-1)*jdim1 + 1 + jdh
         do 6000 i=1,idim1
         iq  = idh+1-i
         do 5500 k=1,5
         do 5500 l=1,5
cdir$ ivdep
        do 1005 izz=1,mm
         t(izz,7) = (dfp(izz+jv2-1,i+1,k,l)-dfm(izz+jv2-1,i,k,l))
 1005    continue
         call q8vrev(mm,t(1,7),mm,bi(jv1,iq,k,l))
         call q8vrev(mm,dfp(jv2,i,k,l),mm,ci(jv1,iq,k,l))
cdir$ ivdep
         do 1006 izz=1,mm
         ci(izz+jv1-1,iq,k,l) = -ci(izz+jv1-1,iq,k,l)
 1006    continue
         call q8vrev(mm,dfm(jv2,i+1,k,l),mm,ai(jv1,iq,k,l))
 5500    continue
         jv3 = jv1+(iq-1)*jdh*npl
         call q8vrev(mm,q(jdh+1,kk,i,1),mm,t(jv3,1))
         call q8vrev(mm,q(jdh+1,kk,i,2),mm,t(jv3,2))
         call q8vrev(mm,q(jdh+1,kk,i,3),mm,t(jv3,3))
         call q8vrev(mm,q(jdh+1,kk,i,4),mm,t(jv3,4))
         call q8vrev(mm,dtj(jdh+1,kk,i),mm,t(jv3,6))
 6000    continue
c
c      assemble matrix equation - time terms
c
         if (real(cprec) .eq. 0.) then
cdir$ ivdep
            do 1007 izz=1,n
            t(izz,6)      = tfacp1*t(izz,6)
            temp          = t(izz,6)*t(izz,1)
            bi(izz,1,1,1) = bi(izz,1,1,1) + t(izz,6)
            bi(izz,1,2,1) = bi(izz,1,2,1) + t(izz,6)*t(izz,2)
            bi(izz,1,2,2) = bi(izz,1,2,2) + temp
            bi(izz,1,3,1) = bi(izz,1,3,1) + t(izz,6)*t(izz,3)
            bi(izz,1,3,3) = bi(izz,1,3,3) + temp
            bi(izz,1,4,1) = bi(izz,1,4,1) + t(izz,6)*t(izz,4)
            bi(izz,1,4,4) = bi(izz,1,4,4) + temp
            bi(izz,1,5,1) = bi(izz,1,5,1)
     .                    + t(izz,6)*.5*(t(izz,2)*t(izz,2)+
     .                                   t(izz,3)*t(izz,3)+
     .                                   t(izz,4)*t(izz,4))
            bi(izz,1,5,2) = bi(izz,1,5,2) + temp*t(izz,2)
            bi(izz,1,5,3) = bi(izz,1,5,3) + temp*t(izz,3)
            bi(izz,1,5,4) = bi(izz,1,5,4) + temp*t(izz,4)
            bi(izz,1,5,5) = bi(izz,1,5,5) + t(izz,6)/gm1
 1007       continue
         else
cdir$ ivdep
            do 10071 izz=1,n
            c2 = gamma*t(izz,5)/t(izz,1)
            c = sqrt(c2)
            ekin = 0.5*(t(izz,2)**2 + t(izz,3)**2 + t(izz,4)**2)
            ho = c2/gm1 + ekin
            vmag1 = 2.0*ekin
            vel2 = ccmax(vmag1,avn*uref**2)
            vel = sqrt(ccmin(c2,vel2))
            vel = cprec*vel + (1.-cprec)*c
            thet = (1.0/vel**2 - 1.0/c2)
            t(izz,6)      = tfacp1*t(izz,6)
            temp          = t(izz,6)*t(izz,1)
            bi(izz,1,1,1) = bi(izz,1,1,1) + t(izz,6)
            bi(izz,1,1,5) = bi(izz,1,1,5) + t(izz,6)*thet
            bi(izz,1,2,1) = bi(izz,1,2,1) + t(izz,6)*t(izz,2)
            bi(izz,1,2,2) = bi(izz,1,2,2) + temp
            bi(izz,1,2,5) = bi(izz,1,2,5) + t(izz,6)*thet*t(izz,2)
            bi(izz,1,3,1) = bi(izz,1,3,1) + t(izz,6)*t(izz,3)
            bi(izz,1,3,3) = bi(izz,1,3,3) + temp
            bi(izz,1,3,5) = bi(izz,1,3,5) + t(izz,6)*thet*t(izz,3)
            bi(izz,1,4,1) = bi(izz,1,4,1) + t(izz,6)*t(izz,4)
            bi(izz,1,4,4) = bi(izz,1,4,4) + temp
            bi(izz,1,4,5) = bi(izz,1,4,5) + t(izz,6)*thet*t(izz,4)
            bi(izz,1,5,1) = bi(izz,1,5,1) + t(izz,6)*ekin
            bi(izz,1,5,2) = bi(izz,1,5,2) + temp*t(izz,2)
            bi(izz,1,5,3) = bi(izz,1,5,3) + temp*t(izz,3)
            bi(izz,1,5,4) = bi(izz,1,5,4) + temp*t(izz,4)
            bi(izz,1,5,5) = bi(izz,1,5,5) + t(izz,6)*(1./gm1 + thet*ho)
10071       continue
         end if
c
      end if
c
      return
      end
