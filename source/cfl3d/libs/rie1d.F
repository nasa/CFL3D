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
      subroutine rie1d(jvdim,t,jv,cl)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Determine far-field boundary data using quasi 1-d 
c     characteristic relations.
c     Modified for Weiss-Smith preconditioning by J.R. Edwards, NCSU
c       cprec = 0 ---> original code used
c             > 0 ---> modified code used
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension t(jvdim,23)
c
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .        tur10(7)
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /fsum/ sref,cref,bref,xmc,ymc,zmc
      common /precond/ cprec,uref,avn
c
c     unsteady quasi-1-d characteristic boundary conditions
c     for inflow/outflow, including a point vortex (point vortex in
c     2D only)
c
c     storage locations 21,22,23 in the t array contain the far-field
c     u, w and speed of sound values
c
      if (real(cprec).eq.0.) then
c
c ---- use Riemann invariants to determine conditions
c
      uf = u0
      vf = v0
      wf = w0
      cf = c0
      do izz=1,jv
        t(izz,21)=u0
        t(izz,22)=w0
        t(izz,23)=c0
      enddo
c   Point vortex hardwire addition:
c   Must be 2-D!  Must be x-z plane!  Uses xmc, zmc for location for vortex
      if (iipv .eq. 1) then
        xcenter=xmc
        zcenter=zmc
        cosa=cos(alpha)
        sina=sin(alpha)
        pi=acos(-1.)
        do izz=1,jv
          xa=t(izz,18)-xcenter
          za=t(izz,19)-zcenter
          call ccf(xa,za,cosa,sina,cl,xmach,t(izz,21),
     +             t(izz,22),t(izz,23),pi)
        enddo
      end if
c
cdir$ ivdep
      do 3117 izz=1,jv
      t(izz,10) = gamma*t(izz,5)/t(izz,1)
 3117 continue
c
      x2gm1=2.e0/gm1
cdir$ ivdep
      do 3118 izz=1,jv
      t(izz,11) = t(izz,6)*t(izz,2)+t(izz,7)*t(izz,3)+t(izz,8)*t(izz,4)
      t(izz,12) = sqrt(t(izz,10))
      t(izz,12) = t(izz,11)+x2gm1*t(izz,12)
      t(izz,13) = t(izz,6)*t(izz,21)+t(izz,7)*vf+t(izz,8)*t(izz,22)
      t(izz,14) = t(izz,13)-x2gm1*t(izz,23)
 3118 continue
c
      x4gm1 = 0.25e0*gm1
cdir$ ivdep
      do 3119 izz=1,jv
      t(izz,15) = .5e0*(t(izz,12)+t(izz,14))
      t(izz,16) = x4gm1*(t(izz,12)-t(izz,14))
      t(izz,17) = t(izz,15)-t(izz,13)
c  put unsteady metrics into t(5)
      t(izz,5)  = t(izz,20)
      t(izz,18) = t(izz,21)+t(izz,6)*t(izz,17)
      t(izz,19) = vf+t(izz,7)*t(izz,17)
      t(izz,20) = t(izz,22)+t(izz,8)*t(izz,17)
 3119 continue
c
      ent0 = gamma*p0/rho0**gamma
cdir$ ivdep
      do 3120 izz=1,jv
      t(izz,9)  = ent0
      t(izz,12) = 1.e0/(t(izz,1)**gm1)
      t(izz,13) = t(izz,15)+t(izz,5)
      t(izz,17) = ccvmgt(t(izz,15)-t(izz,11),t(izz,17),
     .                  (real(t(izz,13)).ge.0.e0))
      t(izz,18) = ccvmgt(t(izz,2)+t(izz,6)*t(izz,17),t(izz,18),
     .                  (real(t(izz,13)).ge.0.e0))
      t(izz,19) = ccvmgt(t(izz,3)+t(izz,7)*t(izz,17),t(izz,19),
     .                  (real(t(izz,13)).ge.0.e0))
      t(izz,20) = ccvmgt(t(izz,4)+t(izz,8)*t(izz,17),t(izz,20),
     .                  (real(t(izz,13)).ge.0.e0))
      t(izz,9)  = ccvmgt(t(izz,10)*t(izz,12),t(izz,9),
     .                  (real(t(izz,13)).ge.0.e0))
 3120 continue
c
      xgm1 = 1.e0/gm1
cdir$ ivdep
      do 3121 izz=1,jv
      t(izz,16) = t(izz,16)*t(izz,16)
      t(izz,1)  = (t(izz,16)/t(izz,9))**xgm1
      t(izz,2)  = t(izz,18)
      t(izz,3)  = t(izz,19)
      t(izz,4)  = t(izz,20)
      t(izz,5)  = t(izz,1)*t(izz,16)/gamma
 3121 continue
c
      else
c
c     use simpler approach (DOES NOT INCLUDE MOVING BOUNDARY
c                           PROVISION AS YET
      do izz=1,jv
         ubar = t(izz,6)*t(izz,2)+t(izz,7)*t(izz,3)+t(izz,8)*t(izz,4)
c
c        outflow boundary
c
         if (real(ubar) .gt. 0.0) then
            t(izz,5) = p0 
         else
            t(izz,1) = rho0*t(izz,5)/p0
            t(izz,2) = u0
            t(izz,3) = v0
            t(izz,4) = w0
         endif
      enddo
c
      endif
c
      return
      end
