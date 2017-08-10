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
      subroutine fluxm(ax,ay,az,are,at,qe,f,n,t,jkpro,nvtq,
     .                 nou,bou,nbuf,ibufdim)
c
c     $Id$
c
c***********************************************************************
c     Computes "negative" parts of the fluxes using the flux-vector-
c     splitting method of van Leer.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension ax(jkpro),ay(jkpro),az(jkpro),are(jkpro),at(jkpro),
     .          t(nvtq,19),qe(nvtq,5),f(nvtq,5)
c
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
c
      g1 = 1.e0/(gamma*gamma-1.e0)
      g2 =  ggm1*0.5e0
      c1 = 1.e0/gamma
      c2 =  ggm1*g1
      c3 =  2.e0*g1
      c4 = 2.e0*c1
      c5 = 4.e0*g1
      c6 = 2.e0*g1*gm1
      c7 = 2.e0-gamma
      c8 = .5e0*gm1
      zeroc = 0.
      onec  = 1.
c
cdir$ ivdep
      do 1000 izz=1,n
      t(izz,2) = ax(izz)*qe(izz,2)+ay(izz)*qe(izz,3)+az(izz)*qe(izz,4)
     .          +at(izz)
      t(izz,3) = 1.e0/qe(izz,1)
      t(izz,1) = t(izz,2)*qe(izz,1)
c
      t(izz,4) = qe(izz,2)*qe(izz,2)+qe(izz,3)*qe(izz,3)+
     .           qe(izz,4)*qe(izz,4)
      t(izz,5) = gamma*qe(izz,5)*t(izz,3)
 1000 continue
c
      dprint   = 0.e0
cdir$ ivdep
      do 1001 izz=1,n
      t(izz,6) = 0.e0
      t(izz,6) = ccvmgt(onec,t(izz,6),(real(t(izz,5)).lt.0.e0))
 1001 continue
c
      dprint   = q8sdot(n,t(1,6),n,t(1,6))
      if (real(dprint).gt.0.e0) then
      nou(1) = min(nou(1)+1,ibufdim)
      write(bou(nou(1),1),1441) real(dprint),n,jkpro
 1441 format(1x,e20.10,2i10)
c     do 7812 in=1,n
c     if (real(t(in,5)).gt.0.e0) go to 7812
c     nou(1) = min(nou(1)+1,ibufdim)
c     write(bou(nou(1),1),7811) in,n,(real(qe(in,l)),l=1,5),
c    .                          real(t(in,5))
c7811 format(1x,3hvfp,2i6,6e11.4)
c7812 continue
      end if
c
      if (real(dprint).gt.0.e0) then
cdir$ ivdep
      do 1002 izz=1,n
      t(izz,6)  = t(izz,5)
      qe(izz,1) = ccvmgt(f(izz,1),qe(izz,1),(real(t(izz,6)).lt.0.e0))
      qe(izz,2) = ccvmgt(f(izz,2),qe(izz,2),(real(t(izz,6)).lt.0.e0))
      qe(izz,3) = ccvmgt(f(izz,3),qe(izz,3),(real(t(izz,6)).lt.0.e0))
      qe(izz,4) = ccvmgt(f(izz,4),qe(izz,4),(real(t(izz,6)).lt.0.e0))
      qe(izz,5) = ccvmgt(f(izz,5),qe(izz,5),(real(t(izz,6)).lt.0.e0))
      t(izz,2)  = ccvmgt(ax(izz)*f(izz,2)+ay(izz)*f(izz,3)+
     .            az(izz)*f(izz,4),t(izz,2),(real(t(izz,6)).lt.0.e0))
      t(izz,3)  = ccvmgt(1.e0/qe(izz,1),t(izz,3),
     .            (real(t(izz,6)).lt.0.e0))
      t(izz,1)  = ccvmgt(t(izz,2)*qe(izz,1),t(izz,1),
     .            (real(t(izz,6)).lt.0.e0))
c
      t(izz,4)  = ccvmgt(qe(izz,2)*qe(izz,2)+qe(izz,3)*qe(izz,3)+
     .            qe(izz,4)*qe(izz,4),t(izz,4),(real(t(izz,6)).lt.0.e0))
      t(izz,5)  = ccvmgt(gamma*qe(izz,5)*t(izz,3),t(izz,5),
     .                 (real(t(izz,6)).lt.0.e0))
 1002 continue
      end if
c
cdir$ ivdep
      do 1003 izz=1,n
      t(izz,6)  = sqrt(t(izz,5))
      qe(izz,5) = qe(izz,5)/gm1+0.5e0*t(izz,4)*qe(izz,1)
      t(izz,7)  = t(izz,2)/t(izz,6)
c
      t(izz,8)  = .5e0*t(izz,7)-0.5e0
      f(izz,1)  = -are(izz)*qe(izz,1)*t(izz,6)*t(izz,8)*t(izz,8)
c
      t(izz,8)  = -c1*(t(izz,2)+2.e0*t(izz,6))
      f(izz,2)  = ax(izz)*t(izz,8)+qe(izz,2)
      f(izz,3)  = ay(izz)*t(izz,8)+qe(izz,3)
      f(izz,4)  = az(izz)*t(izz,8)+qe(izz,4)
c
      f(izz,5)  = c2*t(izz,2)*t(izz,8)+c3*t(izz,5)+.5e0*t(izz,4)
     .               -at(izz)*t(izz,8)
c
      f(izz,2)  = f(izz,1)*f(izz,2)
      f(izz,3)  = f(izz,1)*f(izz,3)
      f(izz,4)  = f(izz,1)*f(izz,4)
      f(izz,5)  = f(izz,1)*f(izz,5)
c
      t(izz,12) = ccabs(t(izz,7))
      t(izz,11) = 0.e0
      t(izz,11) = ccvmgt(onec,t(izz,11),(real(t(izz,12)).ge.+1.e0))
 1003 continue
      skip      = q8sdot(n,t(1,11),n,t(1,11))
      if (real(skip).lt.0.5e0) return
c
cdir$ ivdep
      do 1004 izz=1,n
      t(izz,11) = ccvmgt(ax(izz)*are(izz),t(izz,11),
     .            (real(t(izz,7)).le.-1.e0))
      t(izz,12) = ccvmgt(ay(izz)*are(izz),t(izz,12),
     .            (real(t(izz,7)).le.-1.e0))
      t(izz,13) = ccvmgt(az(izz)*are(izz),t(izz,13),
     .            (real(t(izz,7)).le.-1.e0))
      t(izz,14) = ccvmgt(at(izz)*are(izz),t(izz,14),
     .            (real(t(izz,7)).le.-1.e0))
      t(izz,2)  = ccvmgt(t(izz,2)*are(izz),t(izz,2),
     .            (real(t(izz,7)).le.-1.e0))
c
      f(izz,1) = ccvmgt(t(izz,2)*qe(izz,1),f(izz,1),
     .           (real(t(izz,7)).le.-1.e0))
      t(izz,6) = ccvmgt(qe(izz,1)*t(izz,5)*c1,t(izz,6),
     .           (real(t(izz,7)).le.-1.e0))
      f(izz,2) = ccvmgt(f(izz,1)*qe(izz,2)+t(izz,6)*t(izz,11),f(izz,2),
     .           (real(t(izz,7)).le.-1.e0))
      f(izz,3) = ccvmgt(f(izz,1)*qe(izz,3)+t(izz,6)*t(izz,12),f(izz,3),
     .           (real(t(izz,7)).le.-1.e0))
      f(izz,4) = ccvmgt(f(izz,1)*qe(izz,4)+t(izz,6)*t(izz,13),f(izz,4),
     .           (real(t(izz,7)).le.-1.e0))
      f(izz,5) = ccvmgt(t(izz,2)*(qe(izz,5)+
     .           t(izz,6))-t(izz,14)*t(izz,6),f(izz,5),
     .           (real(t(izz,7)).le.-1.e0))
c
      f(izz,1) = ccvmgt(zeroc,f(izz,1),(real(t(izz,7)).ge.+1.e0))
      f(izz,2) = ccvmgt(zeroc,f(izz,2),(real(t(izz,7)).ge.+1.e0))
      f(izz,3) = ccvmgt(zeroc,f(izz,3),(real(t(izz,7)).ge.+1.e0))
      f(izz,4) = ccvmgt(zeroc,f(izz,4),(real(t(izz,7)).ge.+1.e0))
      f(izz,5) = ccvmgt(zeroc,f(izz,5),(real(t(izz,7)).ge.+1.e0))
 1004 continue
      return
      end
