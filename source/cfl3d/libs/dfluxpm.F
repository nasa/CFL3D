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
      subroutine dfluxpm(s1,s2,ax,ay,az,area,at,q,df,n,jkpro,nvtq,ipm)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Compute "positive" or "negative" parts of the flux
c     Jacobians using the flux-vector-splitting method of van Leer.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension ax(jkpro),ay(jkpro),az(jkpro),at(jkpro),area(jkpro)
      dimension q(nvtq,5),df(n,5,5),s1(jkpro),s2(jkpro)
c
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
c
      gam1   = gamma/gm1
      gam21i = 1.0/(gamma*gamma-1.0)
      zeroc  = 0.
c
      if (ipm.gt.0) sign = +1.0
      if (ipm.lt.0) sign = -1.0
c
cdir$ ivdep
      do 1000 l=1,n
      q2    = q(l,2)*q(l,2)+q(l,3)*q(l,3)+q(l,4)*q(l,4)
      q22   = 0.5*q2
      epp   = gam1*q(l,5)+q(l,1)*q22
      a     = sqrt(gamma*q(l,5)/q(l,1))
      a2    = a*a
      am2   = 2.0*a
      am2i  = 1.0/am2
      ubar  = ax(l)*q(l,2)+ay(l)*q(l,3)+az(l)*q(l,4)+at(l)
      t1    = at(l)/gamma
c
      xm    = sign*ubar/a
      s1(l) = ccvmgt(area(l),zeroc,(abs(real(xm)).lt.1.0))
      s2(l) = ccvmgt(area(l),zeroc,(real(xm).ge.1.0))
c
      ua    =  ubar+sign*a
      ub    = -ubar+sign*a
c
      u1    = ubar*am2i
      u2    = am2-u1*ub
      u3    = ubar/gp1
c
      ua1   = s1(l)*0.25*ua
      ua2   = sign*ua*am2i
      ua3   = 3.0*ua2-1.0
      ua4   = 1.0-ua2
      ua5   = s1(l)*q(l,1)*ua2
      ua6   = ua1/a2
      ua7   = gamma*ua6
      ua8   = 0.5*ua5
      ua8x  = ua8*ax(l)
      ua8y  = ua8*ay(l)
      ua8z  = ua8*az(l)
      ua8a  = ua8*ua
c
      ub1   = 3.0*ub/gamma
      ub2   = ub1*ax(l)
      ub3   = ub2+2.0*q(l,2)
      ub4   = ub1*ay(l)
      ub5   = ub4+2.0*q(l,3)
      ub6   = ub1*az(l)
      ub7   = ub6+2.0*q(l,4)
      ub8   = 2.0*a2/gm1+ub*(4.*u3-3.*t1)+q2
c
      f1    = s2(l)*ubar
      f2    = s2(l)*ax(l)
      f3    = s2(l)*ay(l)
      f4    = s2(l)*az(l)
      f5    = f1*q(l,1)
      f6    = f2*q(l,1)
      f7    = f3*q(l,1)
      f8    = f4*q(l,1)
c
      df(l,1,1) = ua1*ua3                             + f1
      df(l,1,2) = ua5*ax(l)                           + f6
      df(l,1,3) = ua5*ay(l)                           + f7
      df(l,1,4) = ua5*az(l)                           + f8
      df(l,1,5) = ua7*ua4
c
      df(l,2,1) = ua1*(sign*ub2*u1+q(l,2)*ua3)        + f1*q(l,2)
      df(l,2,2) = ua8x*ub3+ua8a                       + f6*q(l,2)+f5
      df(l,2,3) = ua8y*ub3                            + f7*q(l,2)
      df(l,2,4) = ua8z*ub3                            + f8*q(l,2)
      df(l,2,5) = ua6*(sign*ax(l)*u2+gamma*q(l,2)*ua4)+ f2
c
      df(l,3,1) = ua1*(sign*ub4*u1+q(l,3)*ua3)        + f1*q(l,3)
      df(l,3,2) = ua8x*ub5                            + f6*q(l,3)
      df(l,3,3) = ua8y*ub5+ua8a                       + f7*q(l,3)+f5
      df(l,3,4) = ua8z*ub5                            + f8*q(l,3)
      df(l,3,5) = ua6*(sign*ay(l)*u2+gamma*q(l,3)*ua4)+ f3
c
      df(l,4,1) = ua1*(sign*ub6*u1+q(l,4)*ua3)        + f1*q(l,4)
      df(l,4,2) = ua8x*ub7                            + f6*q(l,4)
      df(l,4,3) = ua8y*ub7                            + f7*q(l,4)
      df(l,4,4) = ua8z*ub7+ua8a                       + f8*q(l,4)+f5
      df(l,4,5) = ua6*(sign*az(l)*u2+gamma*q(l,4)*ua4)+ f4
c
      df(l,5,1) = ua1*(sign*ub*(3.0*u1*(u3-t1)-
     .            a*gam21i)+q22*ua3)                  + f1*q22
      df(l,5,2) = ua8x*ub8+q(l,2)*ua8a                + f2*epp+f5*q(l,2)
      df(l,5,3) = ua8y*ub8+q(l,3)*ua8a                + f3*epp+f5*q(l,3)
      df(l,5,4) = ua8z*ub8+q(l,4)*ua8a                + f4*epp+f5*q(l,4)
      df(l,5,5) = ua7*(sign*(u3-t1)*u2+
     .            a*gam21i*(am2+sign*ua)+q22*ua4)  + f1*gam1-s2(l)*at(l)
 1000 continue
      return
      end
