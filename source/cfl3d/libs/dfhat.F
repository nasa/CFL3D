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
      subroutine dfhat(ax,ay,az,area,at,q,df,nn,nvtq,ipm)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Compute a Jacobian matrix with respect to primitive
c     variables at the cell interface.  The Jacobian evaluation is
c     approximate, being taken as either A+ or A- (T lambda Tinv), and
c     is computed with metric terms from the interface and dependent
c     variables from the cell centers.
c     Modified for Weiss-Smith preconditioning by J.R. Edwards, NCSU
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension ax(nn),ay(nn),az(nn),area(nn),at(nn)
      dimension q(nvtq,5),df(nn,5,5)
c
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /precond/ cprec,uref,avn
      common /entfix/ epsa_l,epsa_r
      common /zero/iexp
c
c     10.**(-iexp) is machine zero
      zero    = 10.**(-iexp)
      epsa_l  = 2.*epsa_r
c
      gm1i  = 1.0/gm1
c
      sign = 1.e0
      if (ipm.lt.0) sign = -1.e0
c
cdir$ ivdep
      do 1000 n=1,nn
      rhoi  = 1.0/q(n,1)
      c2    = gamma*q(n,5)*rhoi
      c     = sqrt(c2) 
      c2i   = 1.0/c2
      q22   = 0.5*(q(n,2)*q(n,2)+q(n,3)*q(n,3)+q(n,4)*q(n,4))
      h     = c2*gm1i+q22
      ubar  = q(n,2)*ax(n)+q(n,3)*ay(n)+q(n,4)*az(n)+at(n)
      ar2   = 0.5*area(n)
#if defined V5
c
c
c     limit eigenvalues a la Harten and Gnoffo (NASA TP-2953)
c
      cc    = c
      uu    = ccabs(q(n,2))
      vv    = ccabs(q(n,3))
      ww    = ccabs(q(n,4))
      epsaa = epsa_l*(cc + uu + vv + ww)
      epsbb = 0.25/ccmax(epsaa,zero)
      epscc = 2.00*epsaa
      t17   = ccabs(ubar-c)
      t18   = ccabs(ubar)
      t19   = ccabs(ubar+c)
      if (real(t18).lt.real(epscc)) t18 = t18*t18*epsbb + epsaa
      if (real(t17).lt.real(epscc)) t17 = t17*t17*epsbb + epsaa
      if (real(t19).lt.real(epscc)) t19 = t19*t19*epsbb + epsaa
c
      e1    = ar2*( ubar   +sign*t18)
      e2    = ar2*((ubar+c)+sign*t19)
      e3    = ar2*((ubar-c)+sign*t17)
c
      ruc1  = q(n,1)*e1
      ruc2  = 0.5*q(n,1)*(e2-e3)/c
      ruc3  = 0.5*q(n,1)*(e2+e3-2.0*e1)
c
      rhor2 = rhoi *ruc2
      axr3  = ax(n)*ruc3
      ayr3  = ay(n)*ruc3
      azr3  = az(n)*ruc3
#else
c 
c     preconditioning additions
c
      vmag1 =  2.0*q22
      vel2 = ccmax(vmag1,avn*uref**2)
      vel = sqrt(ccmin(c2,vel2))
      vel = cprec*vel + (1.-cprec)*c
      xm2 = (vel/c)**2   
      xmave = ubar/c
      tt1 = 0.5*(1.+xm2)   
      tt2 = 0.5*sqrt(xmave**2*(1.-xm2)**2 + 4.0*xm2)
c
      e1u    = ubar
      e2u    = tt1*ubar+tt2*c
      e3u    = tt1*ubar-tt2*c
c
c      limit eigenvalues a la Harten and Gnoffo (NASA TP-2953)
c
      cc    = c
      uu    = ccabs(q(n,2))
      vv    = ccabs(q(n,3))
      ww    = ccabs(q(n,4))
      epsaa = epsa_l*(cc + uu + vv + ww)
      epsbb = 0.25/ccmax(epsaa,zero)
      epscc = 2.00*epsaa
      t17   = ccabs(ubar-c)
      t18   = ccabs(ubar)
      t19   = ccabs(ubar+c)
      if (real(t18).lt.real(epscc)) t18 = t18*t18*epsbb + epsaa
      if (real(t17).lt.real(epscc)) t17 = t17*t17*epsbb + epsaa
      if (real(t19).lt.real(epscc)) t19 = t19*t19*epsbb + epsaa
c
      e1    = ar2*(e1u + sign*t18)
      e2    = ar2*(e2u + sign*t19)
      e3    = ar2*(e3u + sign*t17)
      fplus = (e2u-e1u)/(xm2*c)
      fmins = -(e3u-e1u)/(xm2*c)
      fsum = 2.0/(fplus + fmins)/xm2
c
      ruc1  = q(n,1)*e1
      ruc2  = 0.5*q(n,1)*fsum*(e2-e3)/c
      ruc3  = 0.5*q(n,1)*(fsum*(fplus*e2+fmins*e3)-2.0*e1)
      ruc4  = 0.5*q(n,1)*(fsum*xm2*(fmins*e2+fplus*e3)-2.0*e1)
c
      rhor2 = rhoi*ruc2*xm2*fplus*fmins
      axr3  = ax(n)*ruc4
      ayr3  = ay(n)*ruc4
      azr3  = az(n)*ruc4
#endif
c 
      df(n,1,1) = e1
      df(n,1,2) = ax(n)*ruc2
      df(n,1,3) = ay(n)*ruc2
      df(n,1,4) = az(n)*ruc2
      df(n,1,5) = c2i*rhoi*ruc3
c
      df(n,2,1) =      q(n,2)*df(n,1,1)
      df(n,2,2) = ruc1+q(n,2)*df(n,1,2)+ax(n)*axr3
      df(n,2,3) =      q(n,2)*df(n,1,3)+ax(n)*ayr3
      df(n,2,4) =      q(n,2)*df(n,1,4)+ax(n)*azr3
      df(n,2,5) =      q(n,2)*df(n,1,5)+ax(n)*rhor2 
c
      df(n,3,1) =      q(n,3)*df(n,1,1)
      df(n,3,2) =      q(n,3)*df(n,1,2)+ay(n)*axr3
      df(n,3,3) = ruc1+q(n,3)*df(n,1,3)+ay(n)*ayr3
      df(n,3,4) =      q(n,3)*df(n,1,4)+ay(n)*azr3
      df(n,3,5) =      q(n,3)*df(n,1,5)+ay(n)*rhor2
c
      df(n,4,1) =      q(n,4)*df(n,1,1)
      df(n,4,2) =      q(n,4)*df(n,1,2)+az(n)*axr3
      df(n,4,3) =      q(n,4)*df(n,1,3)+az(n)*ayr3
      df(n,4,4) = ruc1+q(n,4)*df(n,1,4)+az(n)*azr3
      df(n,4,5) =      q(n,4)*df(n,1,5)+az(n)*rhor2
c
      df(n,5,1) =   df(n,1,1)*q22
      df(n,5,2) = h*df(n,1,2)+(ubar-at(n))*axr3+q(n,2)*ruc1
      df(n,5,3) = h*df(n,1,3)+(ubar-at(n))*ayr3+q(n,3)*ruc1
      df(n,5,4) = h*df(n,1,4)+(ubar-at(n))*azr3+q(n,4)*ruc1
      df(n,5,5) = h*df(n,1,5)+(ubar-at(n))*rhor2+df(n,1,1)*gm1i
 1000 continue
      return
      end
