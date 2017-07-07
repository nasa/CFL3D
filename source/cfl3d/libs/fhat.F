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
      subroutine fhat(ax,ay,az,area,at,f,qr,ql,n,nvtq)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Compute Roe's generalized flux at the interface
c     given the left and right states at the interface.
c     Extension to Weiss-Smith's preconditioning: J.R. Edwards, July,1998
c       cprec = 0 ---> original code used
c             > 0 ---> modified code used
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension ax(n),ay(n),az(n),area(n),at(n),
     .          qr(nvtq,5),ql(nvtq,5),f(nvtq,5)
c
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /precond/ cprec,uref,avn
      common /entfix/ epsa_l,epsa_r
      common /zero/iexp
c
c     Notes on the entropy fix implemented here and in dfhat, diagj,
c     diagk, and diagi:
c
c     for flow conditions in which the minimum eigenvalue would be
c     less than 2*epsa, epsa nominally forces the minimum eignvalue
c     to be 2*maximum eigenvalue = 2*(|ubar|+c); _l and _r simply
c     denote the limiting values used for the LHS terms and the RHS
c     terms
c
c     below (and on LHS), rather than using |ubar|+c, we use
c     (|u| + |v| + |w| + c), which should add a bit more dissipation.
c     this is slighly different from the implementation in LAURA, 
c     where (|ubar| + |vbar| + |wbar| + c) is used.
c     
c     take LHS limiter = 2 x RHS limiter
c
c     10.**(-iexp) is machine zero
c
      zero    = 10.**(-iexp)
      epsa_l  = 2.*epsa_r
c
c      delta q across faces t(1-5)
c
      x1   = gamma/gm1
      c1   = 1.e0/gm1
c
      if (real(cprec) .eq. 0.) then
cdir$ ivdep
         do 1000 i=1,n 
         t1 = qr(i,1)-ql(i,1)
         t2 = qr(i,2)-ql(i,2)
         t3 = qr(i,3)-ql(i,3)
         t4 = qr(i,4)-ql(i,4)
c
c         pressure and enthalpy
c
         t16     = 1.e0/qr(i,1)
         t5      = qr(i,5)
         qr(i,5) = x1*qr(i,5)*t16+.5e0*(qr(i,2)*qr(i,2)+qr(i,3)*qr(i,3)
     .                                 +qr(i,4)*qr(i,4)) 
c
         t15     = 1.e0/ql(i,1)
         t19     = ql(i,5) 
         ql(i,5) = x1*ql(i,5)*t15+.5e0*(ql(i,2)*ql(i,2)+ql(i,3)*ql(i,3)
     .                                 +ql(i,4)*ql(i,4)) 
c
c         unsplit contributions  f(r)+f(l)
c
         t18 = ax(i)*qr(i,2)+ay(i)*qr(i,3)+az(i)*qr(i,4)+at(i) 
         t17 = ax(i)*ql(i,2)+ay(i)*ql(i,3)+az(i)*ql(i,4)+at(i) 
         t6  = t18*qr(i,1)
         t7  = t17*ql(i,1)
         f1  = t6+t7
         f2  = t6*qr(i,2)+t7*ql(i,2)
         f3  = t6*qr(i,3)+t7*ql(i,3)
         f4  = t6*qr(i,4)+t7*ql(i,4)
         f5  = t6*qr(i,5)+t7*ql(i,5)
         t8  = t5+t19 
         f2  = f2+ax(i)*t8
         f3  = f3+ay(i)*t8
         f4  = f4+az(i)*t8
         f5  = f5-at(i)*t8
c
c         roe averaged variables
c
         t6 = qr(i,1)*t15
         t7 = sqrt(t6) 
         t6 = 1.e0/(1.e0+t7) 
         t8 = t7*t6
c
c        average density
c
         qr(i,1) = ql(i,1)*t7
c
c         u,v,w,h average
c
         t9  = ql(i,2)*t6+qr(i,2)*t8
         t10 = ql(i,3)*t6+qr(i,3)*t8 
         t11 = ql(i,4)*t6+qr(i,4)*t8 
         t12 = ql(i,5)*t6+qr(i,5)*t8 
c
c         extract sound speed
c
         t6  = (t9*t9+t10*t10+t11*t11)*0.5e0
         t7  = gm1*(t12-t6) 
         t8  = sqrt(t7) 
c
         t13 = t9*ax(i)+t10*ay(i)+t11*az(i)
c
c         the variables in t are as follows
c         1-4        delta q1-q4
c         6          q2a
c         7          c2a
c         8          ca
c         9,10,11,12 ua,va,wa,ha
c         13         ubara
c         14,15,16   alpha(1),alpha(4),alpha(5)
c         18,19,17   lambda(1),lambda(4),lambda(5)
c
c         rhoa*delta(ubar) , delta(p)/c2a
c
         ql(i,1) = qr(i,1)*(t18-t17) 
         ql(i,2) = (t5-t19)/t7 
c
         t18 = t13+at(i)
         t18 = ccabs(t18)
         t19 = t13+at(i)+t8
         t19 = ccabs(t19)
         t17 = t13+at(i)-t8
         t17 = ccabs(t17)
c
c        limit eigenvalues a la Harten and Gnoffo (NASA TP-2953)
c
         if (real(epsa_r) .gt.0.) then
            cc    = t8
            uu    = ccabs(t9)
            vv    = ccabs(t10)
            ww    = ccabs(t11)
            epsaa = epsa_r*(cc + uu + vv + ww)
            epsbb = 0.25/ccmax(epsaa,zero)
            epscc = 2.00*epsaa
            if (real(t18).lt.real(epscc)) t18 = t18*t18*epsbb + epsaa
            if (real(t17).lt.real(epscc)) t17 = t17*t17*epsbb + epsaa
            if (real(t19).lt.real(epscc)) t19 = t19*t19*epsbb + epsaa
         end if
c
         t14 = t18*(t1-ql(i,2))
         t15 = .5e0*(ql(i,2)+ql(i,1)/t8)
         t16 = (ql(i,2)-t15)*t17
         t15 = t15*t19
c
         qr(i,2) = t18*(t2*qr(i,1)-ax(i)*ql(i,1))
         qr(i,3) = t18*(t3*qr(i,1)-ay(i)*ql(i,1))
         qr(i,4) = t18*(t4*qr(i,1)-az(i)*ql(i,1))
         qr(i,5) = t9*qr(i,2)+t10*qr(i,3)+t11*qr(i,4)
c
         ql(i,1) = t14+t15+t16
         ql(i,2) = t8*(t15-t16)
c
         f1 = f1-ql(i,1) 
         f2 = f2-ql(i,1)*t9-ax(i)*ql(i,2)-qr(i,2)
         f3 = f3-ql(i,1)*t10-ay(i)*ql(i,2)-qr(i,3) 
         f4 = f4-ql(i,1)*t11-az(i)*ql(i,2)-qr(i,4) 
         f5 = f5-ql(i,1)*t12-t13*ql(i,2)-qr(i,5)+t7*c1*t14
c
c         include factor one-half area
c
         t7     = .5e0*area(i) 
         f(i,1) = t7*f1
         f(i,2) = t7*f2
         f(i,3) = t7*f3
         f(i,4) = t7*f4
         f(i,5) = t7*f5
 1000    continue
      else
cdir$ ivdep
         do 10001 i=1,n 
         t1 = qr(i,1)-ql(i,1)
         t2 = qr(i,2)-ql(i,2)
         t3 = qr(i,3)-ql(i,3)
         t4 = qr(i,4)-ql(i,4)
c
c         pressure and enthalpy
c
         t16     = 1.e0/qr(i,1)
         t5      = qr(i,5)
         qr(i,5) = x1*qr(i,5)*t16+.5e0*(qr(i,2)*qr(i,2)+qr(i,3)*qr(i,3)
     .                                 +qr(i,4)*qr(i,4)) 
c
         t15     = 1.e0/ql(i,1)
         t19     = ql(i,5) 
         ql(i,5) = x1*ql(i,5)*t15+.5e0*(ql(i,2)*ql(i,2)+ql(i,3)*ql(i,3)
     .                                 +ql(i,4)*ql(i,4)) 
c
c         unsplit contributions  f(r)+f(l)
c
         t18 = ax(i)*qr(i,2)+ay(i)*qr(i,3)+az(i)*qr(i,4)+at(i) 
         t17 = ax(i)*ql(i,2)+ay(i)*ql(i,3)+az(i)*ql(i,4)+at(i) 
         t6  = t18*qr(i,1)
         t7  = t17*ql(i,1)
         f1  = t6+t7
         f2  = t6*qr(i,2)+t7*ql(i,2)
         f3  = t6*qr(i,3)+t7*ql(i,3)
         f4  = t6*qr(i,4)+t7*ql(i,4)
         f5  = t6*qr(i,5)+t7*ql(i,5)
         t8  = t5+t19 
         f2  = f2+ax(i)*t8
         f3  = f3+ay(i)*t8
         f4  = f4+az(i)*t8
         f5  = f5-at(i)*t8
         t8t = t8
         delp = t5-t19
c
c         roe averaged variables
c
         t6 = qr(i,1)*t15
         t7 = sqrt(t6) 
         t6 = 1.e0/(1.e0+t7) 
         t8 = t7*t6
c
c        average density
c
         qr(i,1) = ql(i,1)*t7
c
c         u,v,w,h average
c
         t9  = ql(i,2)*t6+qr(i,2)*t8
         t10 = ql(i,3)*t6+qr(i,3)*t8 
         t11 = ql(i,4)*t6+qr(i,4)*t8 
         t12 = ql(i,5)*t6+qr(i,5)*t8 
c
c         extract sound speed
c
         t6  = (t9*t9+t10*t10+t11*t11)*0.5e0
         t7  = gm1*(t12-t6) 
         t8  = sqrt(t7) 
c
         t13 = t9*ax(i)+t10*ay(i)+t11*az(i)
c    
c         compute preconditioning parameters
c
         vmag1 =  2.0*t6 + 2.0*t7*ccabs(delp)/t8t
         vel2 = ccmax(vmag1,avn*uref**2)
         vel = sqrt(ccmin(t7,vel2))
         vel = cprec*vel + (1.-cprec)*t8
         xm2 = (vel/t8)**2   
         xmave = t13/t8
         tt1 = 0.5*(1.+xm2)   
         tt2 = 0.5*sqrt(xmave**2*(1.-xm2)**2 + 4.0*xm2)
c
c         the variables in t are as follows
c         1-4        delta q1-q4
c         6          q2a
c         7          c2a
c         8          ca
c         9,10,11,12 ua,va,wa,ha
c         13         ubara
c         14,15,16   alpha(1),alpha(4),alpha(5)
c         18,19,17   lambda(1),lambda(4),lambda(5)
c
c         rhoa*delta(ubar) , delta(p)/c2a
c
         ql(i,1) = qr(i,1)*(t18-t17) 
         ql(i,2) = (t5-t19)/t7 
c
         t18 = t13+at(i)
         t19 = tt1*t13+tt2*t8+at(i) 
         t17 = tt1*t13-tt2*t8+at(i)
         fplus = (t19-t18)/(xm2*t8)
         fmins = -(t17-t18)/(xm2*t8)
         fsum = 2.0/(fplus + fmins)/xm2
         t18 = ccabs(t18)
         t19 = ccabs(t19)
         t17 = ccabs(t17)
c
c        limit eigenvalues a la Harten and Gnoffo (NASA TP-2953)
c
         if (real(epsa_r) .gt.0.) then
            cc    = t8
            uu    = ccabs(t9)
            vv    = ccabs(t10)
            ww    = ccabs(t11)
            epsaa = epsa_r*(cc + uu + vv + ww)
            epsbb = 0.25/ccmax(epsaa,zero)
            epscc = 2.00*epsaa
            if (real(t18).lt.real(epscc)) t18 = t18*t18*epsbb + epsaa
            if (real(t17).lt.real(epscc)) t17 = t17*t17*epsbb + epsaa
            if (real(t19).lt.real(epscc)) t19 = t19*t19*epsbb + epsaa
         end if
c
         t14 = t18*(t1-ql(i,2))
         t15 = .5e0*t19*(fplus*ql(i,2)+ql(i,1)/t8) 
         t16 = .5e0*t17*(fmins*ql(i,2)-ql(i,1)/t8) 
c
         qr(i,2) = t18*(t2*qr(i,1)-ax(i)*ql(i,1))
         qr(i,3) = t18*(t3*qr(i,1)-ay(i)*ql(i,1))
         qr(i,4) = t18*(t4*qr(i,1)-az(i)*ql(i,1))
         qr(i,5) = t9*qr(i,2)+t10*qr(i,3)+t11*qr(i,4)
c
         ql(i,1) = t14+fsum*(t15+t16)
         ql(i,2) = fsum*t8*xm2*(fmins*t15-fplus*t16)
c
         f1 = f1-ql(i,1) 
         f2 = f2-ql(i,1)*t9-ax(i)*ql(i,2)-qr(i,2)
         f3 = f3-ql(i,1)*t10-ay(i)*ql(i,2)-qr(i,3) 
         f4 = f4-ql(i,1)*t11-az(i)*ql(i,2)-qr(i,4) 
         f5 = f5-ql(i,1)*t12-t13*ql(i,2)-qr(i,5)+t7*c1*t14
c
c         include factor one-half area
c
         t7     = .5e0*area(i) 
         f(i,1) = t7*f1
         f(i,2) = t7*f2
         f(i,3) = t7*f3
         f(i,4) = t7*f4
         f(i,5) = t7*f5
10001    continue
      end if
c
      return
      end
