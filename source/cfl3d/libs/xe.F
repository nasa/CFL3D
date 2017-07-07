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
      subroutine xe(jdim,kdim,nsub,l,x,y,z,xmid,ymid,zmid,xmide,ymide,
     .              zmide,jcell,kcell,xc,yc,zc,xie,eta,imiss,ifit,ic0,
     .              nou,bou,nbuf,ibufdim,myid)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Select proper coordinates to use for inversion (can only
c     use 2 of the 3 available (x,y,z)).
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension x(jdim,kdim,nsub),y(jdim,kdim,nsub),z(jdim,kdim,nsub)
      dimension xmid(jdim,kdim,nsub),ymid(jdim,kdim,nsub),
     .          zmid(jdim,kdim,nsub),a(3)
      dimension xmide(jdim,kdim,nsub),ymide(jdim,kdim,nsub),
     .          zmide(jdim,kdim,nsub)
c
      common/areas/ap(3),imaxa
c
      idum1 = 0
      idum2 = 0
      idum3 = 0
      idum4 = 0
      dum1  = 0.
      dum2  = 0.
      dum3  = 0.
c
      x1 = x(jcell,kcell,l)
      y1 = y(jcell,kcell,l)
      z1 = z(jcell,kcell,l)
      x2 = x(jcell+1,kcell,l)
      y2 = y(jcell+1,kcell,l)
      z2 = z(jcell+1,kcell,l)
      x4 = x(jcell,kcell+1,l)
      y4 = y(jcell,kcell+1,l)
      z4 = z(jcell,kcell+1,l)
      x3 = x(jcell+1,kcell+1,l)
      y3 = y(jcell+1,kcell+1,l)
      z3 = z(jcell+1,kcell+1,l)
      x5 = xmid(jcell,kcell,l)
      y5 = ymid(jcell,kcell,l)
      z5 = zmid(jcell,kcell,l)
      x6 = xmid(jcell,kcell+1,l)
      y6 = ymid(jcell,kcell+1,l)
      z6 = zmid(jcell,kcell+1,l)
      x7 = xmide(jcell,kcell,l)
      y7 = ymide(jcell,kcell,l)
      z7 = zmide(jcell,kcell,l)
      x8 = xmide(jcell+1,kcell,l)
      y8 = ymide(jcell+1,kcell,l)
      z8 = zmide(jcell+1,kcell,l)
c
c     compute normalized directed areas/unit normals of "from" cell 
c     (corresponding values for "to" cell are stored in common/areas/)
c
      call direct(x5,x6,x7,x8,y5,y6,y7,y8,z5,z6,z7,z8,
     .                   a1,a2,a3,itoss,nou,bou,nbuf,ibufdim)
      a(1)=a1
      a(2)=a2
      a(3)=a3
c
c     project current "from" cell node points onto plane defined by "to" cell
c
c     check inner product of unit normals for compatable orientation of "to" 
c     and "from" cells before projection.  If not compatable, exit to try 
c     a nearby "from" cell 
c
      tol = .1
      prod = ap(1)*a(1) + ap(2)*a(2) + ap(3)*a(3)
c
      if(abs(real(prod)) .gt. real(tol))then
        if(ic0.eq.0)then
          call project(xc,yc,zc,x1,y1,z1,ap(1),ap(2),ap(3),x1,y1,z1)      
          call project(xc,yc,zc,x2,y2,z2,ap(1),ap(2),ap(3),x2,y2,z2)      
          call project(xc,yc,zc,x3,y3,z3,ap(1),ap(2),ap(3),x3,y3,z3)      
          call project(xc,yc,zc,x4,y4,z4,ap(1),ap(2),ap(3),x4,y4,z4)      
          call project(xc,yc,zc,x5,y5,z5,ap(1),ap(2),ap(3),x5,y5,z5)      
          call project(xc,yc,zc,x6,y6,z6,ap(1),ap(2),ap(3),x6,y6,z6)      
          call project(xc,yc,zc,x7,y7,z7,ap(1),ap(2),ap(3),x7,y7,z7)      
          call project(xc,yc,zc,x8,y8,z8,ap(1),ap(2),ap(3),x8,y8,z8)      
        end if
      else
c
c     call trace(8,idum1,idum2,idum3,idum4,dum1,dum2,dum3)
c
        imiss = 1
        xie = xie + 1.
        eta = eta + 1.
        return
      end if
c  
c     use best coodinates to for inversion in "to" cell... toss out
c     equation for which directed area ap(i) is maximum.  
c     special cases:
c      a) two values of ap(i) are equal, and greater than the third.  In this
c         case, the equation corresponding to the first max is thrown out.
c      b) all three values of ap(i) are equal.  In this case, the x-equation
c         is tossed out, and inversion performed using y-z coordinates.   
c
      itoss=imaxa
c
      if (itoss.eq.1) then
c
c     use only y and z equations
c
c     call trace(9,l,jcell,kcell,idum4,dum1,dum2,dum3)
c
      call xe2(y1,y2,y3,y4,y5,y6,y7,y8,yc,z1,z2,z3,z4,z5,
     .         z6,z7,z8,zc,xie,eta,imiss,ifit,nou,bou,nbuf,
     .         ibufdim,myid)
      end if
c
      if (itoss.eq.2) then
c
c     use only x and z equations
c
c     call trace(10,l,jcell,kcell,idum4,dum1,dum2,dum3)
c
      call xe2(x1,x2,x3,x4,x5,x6,x7,x8,xc,z1,z2,z3,z4,z5,
     .         z6,z7,z8,zc,xie,eta,imiss,ifit,nou,bou,nbuf,
     .         ibufdim,myid)
      end if
c
      if (itoss.eq.3) then
c
c     use only x and y equations
c
c     call trace(11,l,jcell,kcell,idum4,dum1,dum2,dum3)
c
      call xe2(x1,x2,x3,x4,x5,x6,x7,x8,xc,y1,y2,y3,y4,y5,
     .         y6,y7,y8,yc,xie,eta,imiss,ifit,nou,bou,nbuf,
     .         ibufdim,myid)
      end if
      return
      end
