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
      subroutine ffluxv(k,npl,jdim,kdim,idim,res,q,qi0,si,vol,t,nvtq,
     .                  wi0,vist3d,vmui,vi0,bci,zksav,ti0,cmuv,voli0,
     .                  nou,bou,nbuf,ibufdim,iadv,nummem,ux)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Calculate right-hand residual contributions in the 
c     I-direction due to the viscous terms.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension si(jdim,kdim,idim,5),vol(jdim,kdim,idim-1)
      dimension q(jdim,kdim,idim,5),qi0(jdim,kdim,5,4)
      dimension res(jdim,kdim,idim-1,5),vist3d(jdim,kdim,idim),
     .          vi0(jdim,kdim,1,4),voli0(jdim,kdim,4)
      dimension t(nvtq,32),vmui(jdim-1,kdim-1,2),wi0(npl*(jdim-1),22)
      dimension bci(jdim,kdim,2)
      dimension zksav(jdim,kdim,idim,nummem),ti0(jdim,kdim,nummem,4)
      dimension cmuv(jdim-1,kdim-1,idim-1)
      dimension ux(jdim-1,kdim-1,idim-1,9)
c
      common /easmv/ c10,c11,c2,c3,c4,c5,sigk1,cmuc1,ieasm_type
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /fluid2/ pr,prt,cbar
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .        tur10(7)
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /reyue/ reue,tinf,ivisc(3)
      common /sklton/ isklton
      common /constit/ i_nonlin,c_nonlin,snonlin_lim,i_tauijs,i_qcr2000,
     .        i_qcr2013,i_qcr2013v
      real :: coef_eddy
c
c     coefficient used to turn off eddy-viscosity based
c     turbulent stress calculation for ivmx==70 (full
c     Reynolds stress model)  - xxiao
c
      coef_eddy = 1.0
      if(ivisc(1)>=70) coef_eddy=0.0
      if(i_tauijs .eq. 1) coef_eddy=0.
      ccr2=2.5
c
      idim1  =  idim-1
      jdim1  =  jdim-1
c
c n  : number of cell centers for npl planes
c l0 : number of cell interfaces for npl planes
c jv : number of cell centers (and interfaces) on a i=constant plane
c
      n      =  npl*jdim1*idim1
      l0     =  npl*jdim1*idim
      jv     =  npl*jdim1
      js     =  jv+1
      nn     =  n-jv
c
      xmre   =  xmach/reue
      gpr    =  gamma/pr
      gm1pr  =  gm1*pr
      prtr   =  pr/prt
c     Durbin TCFD 1991 near-wall limiter (0=off)
c     (using it tends to delay separation - generally not desired!)
      idurbinlim=0
c
      if (isklton.gt.0 .and. k.eq.1 .and. iadv.ge.0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1444)
         if (i_tauijs .eq. 1) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''      using tau_ij method'')')
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''(note: not done accurately at '',
     +       ''boundaries for multiblock cases...'')')
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''currently one sided; need ux '',
     +       ''derivs from across boundaries)'')')
         end if
         if (i_qcr2000 .eq. 1) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''      using QCR2000'')')
         end if
         if (i_qcr2013 .eq. 1) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''      using QCR2013 (limited)'')')
         end if
         if (i_qcr2013v .eq. 1) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''      using QCR2013-V'')')
         end if
         if (ivisc(1).eq.11 .or. ivisc(1).eq.12 .or.
     +       ivisc(1).eq.13 .or. ivisc(1).eq.14 .or.
     +       i_nonlin.ne.0) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),1445)
         end if
         if (i_nonlin.ne.0) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''      snonlin_lim='',e19.8)')
     +       snonlin_lim
         end if
      end if
 1444 format(52h   computing viscous  fluxes, I-direction - central ,
     .12hdifferencing)
 1445 format(37h      nonlinear constitutive relation)
c
c******************************************************************************
c
c  store selected cell centered data in t array
c  t(1)    : 1/density
c  t(16)   : c*c/(gm1*pr)
c  t(18)   : u
c  t(19)   : v
c  t(20)   : w
c  t(22)   : density
c
      l1            =  jdim*idim-1
      do 10 kpl=1,npl
      kk            =  k+kpl-1
      do 50 i=1,idim1
      jiv           =  (i-1)*npl*jdim1 + (kpl-1)*jdim1 + 1
      do 1002 j=1,jdim1
      t(j+jiv-1,22) = q(j,kk,i,1)
      t(j+jiv-1,18) = q(j,kk,i,2)
      t(j+jiv-1,19) = q(j,kk,i,3)
      t(j+jiv-1,20) = q(j,kk,i,4)
      t(j+jiv-1,16) = q(j,kk,i,5) 
 1002 continue
   50 continue
   10 continue
c
c******************************************************************************
c
c  store dependent variables on i=0 and i=idim boundaries in wi0
c
      do 1127 m=0,10,10
      if (m.eq.0) then
        mm=1
        mb=1
      else
        mm=3
        mb=2
      end if
c
      jmin = 1
      jmax = jdim-1
      kmin = k
      kmax = k+npl-1
c
      do 6120 kk=kmin,kmax
      jc         = (kk-k)*jdim1 + jmin-1
      do 6120 jj=jmin,jmax
      jc         = jc + 1
      wi0(jc,m+12) = qi0(jj,kk,1,mm)
      wi0(jc,m+8)  = qi0(jj,kk,2,mm)
      wi0(jc,m+9)  = qi0(jj,kk,3,mm)
      wi0(jc,m+10) = qi0(jj,kk,4,mm) 
      wi0(jc,m+6)  = qi0(jj,kk,5,mm) 
c     m+5 flag indicates whether ghost cell or interface values stored in wi0
      wi0(jc,m+5)  = bci(jj,kk,mb)
 6120 continue
 6400 continue
 1127 continue
c
c******************************************************************************
c
c  t(17)/wi0(m+7)   : (u*u+v*v+w*w)/2
c  t(1)/wi0(m+1)    : 1/density
c  t(16)/wi0(m+6)   : c*c/(gm1*pr)
c  note: m=0 gives values at i=0, m=10 gives values at i=idim
c
      do 7005 j=1,n
      t(j,17) = 0.5e0*( t(j,18)*t(j,18)
     .                + t(j,19)*t(j,19)
     .                + t(j,20)*t(j,20) ) 
      t(j,21) = t(j,17)*t(j,22)+t(j,16)/gm1 
      t(j,1)  = 1.e0/t(j,22) 
      t(j,16) = gpr*t(j,16)*t(j,1)/gm1
 7005 continue
      do 2121 m=0,10,10
      do 1005 j=1,jv
      wi0(j,m+7) = 0.5e0*( wi0(j,m+8)*wi0(j,m+8)
     .                  + wi0(j,m+9)*wi0(j,m+9)
     .                  + wi0(j,m+10)*wi0(j,m+10) )
      wi0(j,m+11) = wi0(j,m+7)*wi0(j,m+12)+wi0(j,m+6)/gm1 
      wi0(j,m+1)  = 1.e0/wi0(j,m+12) 
      wi0(j,m+6) = gpr*wi0(j,m+6)*wi0(j,m+1)/gm1
 1005 continue
 2121 continue
c
c******************************************************************************
c
c t(7) : laminar viscosity at cell centers (via sutherland relation)
c
      c2b    =  cbar/tinf
      c2bp   =  c2b+1.e0
c
      do 1006 j=1,n
      t5     =  gm1pr*t(j,16)
      t6     =  sqrt(t5) 
      t(j,7) =  c2bp*t5*t6/(c2b+t5)
 1006 continue
c
c******************************************************************************
c
c t(14) : laminar viscosity values at cell interfaces
c
c     interior interfaces
      do 1007 j=1,nn
      t(j+js-1,14) = ( t(j,7)+t(j+js-1,7) )*.5e0
 1007 continue
c
c     i=0 and i=idim interfaces
cdir$ ivdep
      do 1008 izz=1,jv
      ab=1.+wi0(izz,5)
      bb=1.-wi0(izz,5)
      wi05        = gm1pr*0.5*(ab*wi0(izz,6)+bb*t(izz,16))
      wi06        = sqrt(wi05) 
      t(izz,14)   = c2bp*wi05*wi06/(c2b+wi05) 
      ab=1.+wi0(izz,15)
      bb=1.-wi0(izz,15)
      wi05        = gm1pr*0.5*(ab*wi0(izz,16)+bb*t(izz+n-jv,16))
      wi06        = sqrt(wi05) 
      t(izz+n,14) = c2bp*wi05*wi06/(c2b+wi05) 
 1008 continue
c
c******************************************************************************
c
c t(15) : average jacobian (inverse volume) at cell interface
c
      n1  = jdim*idim1 + 1
      l1  = jdim*idim-1
      do 200 kpl=1,npl
      kk  = k+kpl-1
      do 200 i=1,idim
      jiv = (i-1)*npl*jdim1 + (kpl-1)*jdim1 
      ji  = (i-1)*jdim 
      if (i.eq.1) then 
c     inverse volume at i=0 interface
        do 7013 j=1,jdim1 
        t(j+jiv,15) = 2./(voli0(j,kk,1)+vol(j,kk,1))
 7013   continue
      else if (i.eq.idim) then 
c     inverse volume at i=idim interface
        do 7113 j=1,jdim1 
        t(j+jiv,15) = 2./(voli0(j,kk,3)+vol(j,kk,idim1))
 7113   continue
      else
c     inverse volume at interior interfaces
        do 1013 j=1,jdim1 
        t(j+jiv,15) = 2./(vol(j,kk,i)+vol(j,kk,i-1))
 1013   continue
      end if
c
c******************************************************************************
c
c  t(25) - t(27) : components of grad(xie)
c  t1    : grad(xie)/J
c  t(25) : d(xie)/dx
c  t(26) : d(xie)/dy
c  t(27) : d(xie)/dz
c
      do 1014 j=1,jdim1 
      t1         =  si(j,kk,i,4)*t(j+jiv,15) 
      t(j+jiv,25) =  si(j,kk,i,1)*t1
      t(j+jiv,26) =  si(j,kk,i,2)*t1
      t(j+jiv,27) =  si(j,kk,i,3)*t1
 1014 continue
  200 continue
c
c******************************************************************************
c
c  t(6)-t(10) : gradients at cell interfaces
c  t(6)       : d( c*c/(gm1*pr) )/d(xie)
c  t(7)       : d( (u*u+v*v+w*w)/2 )/d(xie)
c  t(8)       : d(u)/d(xie)
c  t(9)       : d(v)/d(xie)
c  t(10)      : d(w)/d(xie)
c
      do 1000 kq=1,5
      k1           = kq+5
      k2           = kq+15
cdir$ ivdep
c     interior interfaces
      do 1015 j=1,nn
      t(j+js-1,k1) = t(j+js-1,k2)-t(j,k2)
 1015 continue
c     interfaces at i=0/idim
      do 1016 j=1,jv
      z1 = 1.+wi0(j,5)
      z2 = 1.+wi0(j,15)
      t(j,k1) = z1*(t(j,k2) - wi0(j,k1))
      t(j+n,k1) = z2*(wi0(j,k2) - t(j+n-jv,k2))
 1016 continue
 1000 continue
c
c******************************************************************************
c
c  t(24) : turbulent viscosity at interfaces (=0 for laminar flow)
c
      iviscc = ivisc(1)
      if (iviscc.gt.1) then
      do 1401 kpl=1,npl
      kk = k+kpl-1
      do 1401 i=1,idim
      jiv=(i-1)*npl*jdim1 + (kpl-1)*jdim1
c     i=0 interface
      if (i.eq.1) then
        do 1017 j=1,jdim1
          t(j+jiv,24)=bci(j,kk,1)*vi0(j,kk,1,1)+
     +     (1.-bci(j,kk,1))*0.5*(vi0(j,kk,1,1)+vist3d(j,kk,1))
 1017   continue
c     i=idim interface
      else if (i.eq.idim) then
        do 1117 j=1,jdim1
          t(j+jiv,24)=bci(j,kk,2)*vi0(j,kk,1,3)+
     +     (1.-bci(j,kk,2))*0.5*(vi0(j,kk,1,3)+vist3d(j,kk,idim1))
 1117   continue
      else
c     interior interfaces
        do 1217 j=1,jdim1
          t(j+jiv,24)=0.5*(vist3d(j,kk,i)+vist3d(j,kk,i-1))
 1217   continue
      end if
 1401 continue
      else
c     laminar
      do 1018 izz=1,l0
        t(izz,24)=0.
 1018 continue
      end if
c
      do 1022 j=1,l0
c  ratio of turbulent to laminar viscosity
      t25     = t(j,24)/t(j,14) 
      t24     = 1.e0+t25*coef_eddy
c
c******************************************************************************
c
c  t(23) : ratio of laminar Prandtl number to total Prandtl number
c
      t(j,23) = (1.e0+prtr*t25)/t24
c
c******************************************************************************
c
c  t(14): (mach/re)*viscosity/J
c
      t(j,14) = xmre*t24*t(j,14)/t(j,15)
c
c******************************************************************************
c  t(5) : t(14)*(grad(zeta))**2
c
      t(j,5)  = t(j,14)*(t(j,25)*t(j,25)
     .                  +t(j,26)*t(j,26)
     .                  +t(j,27)*t(j,27))
 1022 continue
c
c******************************************************************************
c
c  t(11) : t(14)*((d(xie)/dx * d(u)/d(xie))
c               + (d(xie)/dy * d(v)/d(xie))
c               + (d(xie)/dz * d(w)/d(xie)))/3
c
      tr1 = 1.e0/3.e0
      do 1023 j=1,l0
      t(j,11) = tr1*t(j,14)*( t(j,25)*t(j,8)
     .                      + t(j,26)*t(j,9)
     .                      + t(j,27)*t(j,10) )
 1023 continue
c
c******************************************************************************
c
c  store off vmui....viscosity coefficient (laminar + turbulent) on i=0 and
c  i=idim boundaries for later use in force integration
      do 1923 kpl=1,npl
      kk=k+kpl-1
cdir$ ivdep
      do 1922 j=1,jdim1
      izz = (kpl-1)*jdim1+j
      vmui(j,kk,1) = t(izz,14)*t(izz,15)/xmre
      vmui(j,kk,2) = t(izz+n,14)*t(izz+n,15)/xmre
 1922 continue
 1923 continue
c
      if (iadv .lt. 0) return
c
c******************************************************************************
c
c  t(12) : contravariant velocity U at cell interfaces
c
c     interior interfaces
      do 1024 j=1,nn
      t(j+js-1,12) = 0.5e0*( t(j+js-1,25)*(t(j+js-1,18)+t(j,18))
     .                     + t(j+js-1,26)*(t(j+js-1,19)+t(j,19))
     .                     + t(j+js-1,27)*(t(j+js-1,20)+t(j,20)) )
 1024 continue
c     interfaces at i=0 and i=idim
      do 1025 izz=1,jv
      ab=1.+wi0(izz,15)
      bb=1.-wi0(izz,15)
      t(izz+n,12)=0.5*(t(izz+n,25)*(ab*wi0(izz,18)+bb*t(izz+n-jv,18))
     .                +t(izz+n,26)*(ab*wi0(izz,19)+bb*t(izz+n-jv,19))
     .                +t(izz+n,27)*(ab*wi0(izz,20)+bb*t(izz+n-jv,20)))
      ab=1.+wi0(izz,5)
      bb=1.-wi0(izz,5)
      t(izz,12) = 0.5*(t(izz,25)*(ab*wi0(izz,8) +bb*t(izz,18))
     .                +t(izz,26)*(ab*wi0(izz,9) +bb*t(izz,19))
     .                +t(izz,27)*(ab*wi0(izz,10)+bb*t(izz,20)))
 1025 continue
c
c******************************************************************************
c
c  calculate fluxes
c
      do 3000 l=2,4
c
c  viscous terms at interfaces (momentum eqs.)
c
      do 1026 j=1,l0
      t(j,16) = t(j,5)*t(j,l+6) + t(j,l+23)*t(j,11)
 1026 continue
c
c  (-)viscous flux =  Fv(k-1/2)-Fv(k+1/2) (momentum eqs.)
c
      do 1027 j=1,n
      t(j,l)  = -t(j+js-1,16) + t(j,16)
 1027 continue
 3000 continue
c
c  viscous terms at interfaces (energy eq.)
c
      do 1028 j=1,l0
      t(j,16) = t(j,5)*(t(j,23)*t(j,6)+t(j,7))
     .        + t(j,12)*t(j,11)
 1028 continue
c
c  (-)viscous flux = Fv(k-1/2)-Fv(k+1/2) (energy eq.)
c
      do 1029 j=1,n
      t(j,5) = -t(j+js-1,16) + t(j,16)
 1029 continue
c
c******************************************************************************
c
c  calculate residuals
c
      do 400 kpl=1,npl
      kk  = k+kpl-1
      do 400 i=1,idim1
      jiv = (i-1)*jdim1*npl + (kpl-1)*jdim1 + 1
      do 400 l=2,5
      do 1030 j=1,jdim1
      res(j,kk,i,l) = res(j,kk,i,l) + t(j+jiv-1,l)
 1030 continue
  400 continue
c  Add tau_ij's for general method:
      if (i_tauijs .eq. 1 .and. ivisc(1) .lt. 70) then
        re=reue/xmach
        do 1196 kpl=1,npl
          kk=k+kpl-1
          jiv2=(kpl-1)*jdim1
          do 1195 i=1,idim
            jiv=(i-1)*npl*jdim1 + (kpl-1)*jdim1
            ip1=min(i+1,idim-1)
            im2=max(i-2,1)
            do 1194 j=1,jdim1
              if (i.eq.1) then
                rho=0.5*(q(j,kk,i,1)+qi0(j,kk,1,1))*(1.-wi0(j+jiv2,5))+
     +              qi0(j,kk,1,1)*wi0(j+jiv2,5)
                u=0.5*(q(j,kk,i,2)+qi0(j,kk,2,1))*(1.-wi0(j+jiv2,5))+
     +              qi0(j,kk,2,1)*wi0(j+jiv2,5)
                v=0.5*(q(j,kk,i,3)+qi0(j,kk,3,1))*(1.-wi0(j+jiv2,5))+
     +              qi0(j,kk,3,1)*wi0(j+jiv2,5)
                w=0.5*(q(j,kk,i,4)+qi0(j,kk,4,1))*(1.-wi0(j+jiv2,5))+
     +              qi0(j,kk,4,1)*wi0(j+jiv2,5)
             v3d=0.5*(vist3d(j,kk,i)+vi0(j,kk,1,1))*(1.-wi0(j+jiv2,5))+
     +              vi0(j,kk,1,1)*wi0(j+jiv2,5)
c               use 1st order interpolation at end
                ux1=1.5*ux(j,kk,i,1)-0.5*ux(j,kk,ip1,1)
                ux2=1.5*ux(j,kk,i,2)-0.5*ux(j,kk,ip1,2)
                ux3=1.5*ux(j,kk,i,3)-0.5*ux(j,kk,ip1,3)
                ux4=1.5*ux(j,kk,i,4)-0.5*ux(j,kk,ip1,4)
                ux5=1.5*ux(j,kk,i,5)-0.5*ux(j,kk,ip1,5)
                ux6=1.5*ux(j,kk,i,6)-0.5*ux(j,kk,ip1,6)
                ux7=1.5*ux(j,kk,i,7)-0.5*ux(j,kk,ip1,7)
                ux8=1.5*ux(j,kk,i,8)-0.5*ux(j,kk,ip1,8)
                ux9=1.5*ux(j,kk,i,9)-0.5*ux(j,kk,ip1,9)
              else if (i .eq. idim) then
                rho=0.5*(q(j,kk,i-1,1)+qi0(j,kk,1,3))*
     +              (1.-wi0(j+jiv2,15))+
     +              qi0(j,kk,1,3)*wi0(j+jiv2,15)
                u=0.5*(q(j,kk,i-1,2)+qi0(j,kk,2,3))*(1.-wi0(j+jiv2,15))+
     +              qi0(j,kk,2,3)*wi0(j+jiv2,15)
                v=0.5*(q(j,kk,i-1,3)+qi0(j,kk,3,3))*(1.-wi0(j+jiv2,15))+
     +              qi0(j,kk,3,3)*wi0(j+jiv2,15)
                w=0.5*(q(j,kk,i-1,4)+qi0(j,kk,4,3))*(1.-wi0(j+jiv2,15))+
     +              qi0(j,kk,4,3)*wi0(j+jiv2,15)
           v3d=0.5*(vist3d(j,kk,i-1)+vi0(j,kk,1,3))*(1.-wi0(j+jiv2,15))+
     +              vi0(j,kk,1,3)*wi0(j+jiv2,15)
c               use 1st order interpolation at end
                ux1=1.5*ux(j,kk,i-1,1)-0.5*ux(j,kk,im2,1)
                ux2=1.5*ux(j,kk,i-1,2)-0.5*ux(j,kk,im2,2)
                ux3=1.5*ux(j,kk,i-1,3)-0.5*ux(j,kk,im2,3)
                ux4=1.5*ux(j,kk,i-1,4)-0.5*ux(j,kk,im2,4)
                ux5=1.5*ux(j,kk,i-1,5)-0.5*ux(j,kk,im2,5)
                ux6=1.5*ux(j,kk,i-1,6)-0.5*ux(j,kk,im2,6)
                ux7=1.5*ux(j,kk,i-1,7)-0.5*ux(j,kk,im2,7)
                ux8=1.5*ux(j,kk,i-1,8)-0.5*ux(j,kk,im2,8)
                ux9=1.5*ux(j,kk,i-1,9)-0.5*ux(j,kk,im2,9)
              else
                rho=0.5*(q(j,kk,i,1)+q(j,kk,i-1,1))
                u=0.5*(q(j,kk,i,2)+q(j,kk,i-1,2))
                v=0.5*(q(j,kk,i,3)+q(j,kk,i-1,3))
                w=0.5*(q(j,kk,i,4)+q(j,kk,i-1,4))
                v3d=0.5*(vist3d(j,kk,i)+vist3d(j,kk,i-1))
                ux1=0.5*(ux(j,kk,i,1)+ux(j,kk,i-1,1))
                ux2=0.5*(ux(j,kk,i,2)+ux(j,kk,i-1,2))
                ux3=0.5*(ux(j,kk,i,3)+ux(j,kk,i-1,3))
                ux4=0.5*(ux(j,kk,i,4)+ux(j,kk,i-1,4))
                ux5=0.5*(ux(j,kk,i,5)+ux(j,kk,i-1,5))
                ux6=0.5*(ux(j,kk,i,6)+ux(j,kk,i-1,6))
                ux7=0.5*(ux(j,kk,i,7)+ux(j,kk,i-1,7))
                ux8=0.5*(ux(j,kk,i,8)+ux(j,kk,i-1,8))
                ux9=0.5*(ux(j,kk,i,9)+ux(j,kk,i-1,9))
              end if
              s11=ux1
              s22=ux5
              s33=ux9
              s11star=s11-((s11+s22+s33)/3.)
              s22star=s22-((s11+s22+s33)/3.)
              s33star=s33-((s11+s22+s33)/3.)
              s12=0.5*(ux2+ux4)
              s13=0.5*(ux3+ux7)
              s23=0.5*(ux6+ux8)
              w12=0.5*(ux2-ux4)
              w13=0.5*(ux3-ux7)
              w23=0.5*(ux6-ux8)
c             ignoring -2/3*rho*k*delij term
              t11=-2.0*v3d*s11star
              t22=-2.0*v3d*s22star
              t33=-2.0*v3d*s33star
              t12=-2.0*v3d*s12
              t13=-2.0*v3d*s13
              t23=-2.0*v3d*s23
              if (i_qcr2000 .eq. 1 .or. i_qcr2013 .eq. 1 .or.
     +            i_qcr2013v .eq. 1) then
                denom=sqrt(ux1*ux1 + ux2*ux2 + ux3*ux3 +
     +                     ux4*ux4 + ux5*ux5 + ux6*ux6 +
     +                     ux7*ux7 + ux8*ux8 + ux9*ux9)+1.e-20
                o11=0.
                o22=0.
                o33=0.
                o12=2.*w12/denom
                o13=2.*w13/denom
                o23=2.*w23/denom
                t11p=t11-0.3*( o11*t11+o12*t12+o13*t13+
     +                         o11*t11+o12*t12+o13*t13)
                t22p=t22-0.3*(-o12*t12+o22*t22+o23*t23+
     +                        -o12*t12+o22*t22+o23*t23)
                t33p=t33-0.3*(-o13*t13-o23*t23+o33*t33+
     +                        -o13*t13-o23*t23+o33*t33)
                t12p=t12-0.3*( o11*t12+o12*t22+o13*t23+
     +                        -o12*t11+o22*t12+o23*t13)
                t13p=t13-0.3*( o11*t13+o12*t23+o13*t33+
     +                        -o13*t11-o23*t12+o33*t13)
                t23p=t23-0.3*(-o12*t13+o22*t23+o23*t33+
     +                        -o13*t12-o23*t22+o33*t23)
                t11=t11p
                t22=t22p
                t33=t33p
                t12=t12p
                t13=t13p
                t23=t23p
              end if
              if (i_qcr2013 .eq. 1) then
                xis2013 = s11star*s11star + s22star*s22star + 
     +                    s33star*s33star +
     +                    2.*s12*s12 + 2.*s13*s13 + 2.*s23*s23
                s_mod = sqrt(2.*xis2013)
                wis = 2.*w12*w12 + 2.*w13*w13 + 2.*w23*w23
c               limiting is typically necessary
                w_mod = sqrt(2.*wis)
                s_mod = ccmin(s_mod , 1.2*w_mod)
                t11=t11+ccr2*v3d*s_mod
                t22=t22+ccr2*v3d*s_mod
                t33=t33+ccr2*v3d*s_mod
              end if
              if (i_qcr2013v .eq. 1) then
                wis = 2.*w12*w12 + 2.*w13*w13 + 2.*w23*w23
                w_mod = sqrt(2.*wis)
                t11=t11+ccr2*v3d*w_mod
                t22=t22+ccr2*v3d*w_mod
                t33=t33+ccr2*v3d*w_mod
              end if
              t(j+jiv,28)=xmre*(t(j+jiv,25)*t11+t(j+jiv,26)*t12+
     +                          t(j+jiv,27)*t13)/t(j+jiv,15)
              t(j+jiv,29)=xmre*(t(j+jiv,25)*t12+t(j+jiv,26)*t22+
     +                          t(j+jiv,27)*t23)/t(j+jiv,15)
              t(j+jiv,30)=xmre*(t(j+jiv,25)*t13+t(j+jiv,26)*t23+
     +                          t(j+jiv,27)*t33)/t(j+jiv,15)
              t(j+jiv,31)=xmre/t(j+jiv,15)*(
     +           t(j+jiv,25)*(u*t11 + v*t12 + w*t13)
     +          +t(j+jiv,26)*(u*t12 + v*t22 + w*t23)
     +          +t(j+jiv,27)*(u*t13 + v*t23 + w*t33))
 1194       continue
 1195     continue
 1196   continue
        do 1171 l=2,5
        do 1170 izz=1,n
c   Take positive derivative, because want negative of all above for correct
c   sign (this term is going into -d/dxj(Hv) term)
          t(izz,l)=t(izz+js-1,l+26)-t(izz,l+26)
 1170   continue
 1171   continue
        do 1168 kpl=1,npl
          kk=k+kpl-1
          do 1168 i=1,idim1
            jiv=(i-1)*jdim1*npl + (kpl-1)*jdim1+1
            do 1168 l=2,5
              do 1167 izz=1,jdim1
                res(izz,kk,i,l)=res(izz,kk,i,l)+t(izz+jiv-1,l)
 1167         continue
 1168   continue
      end if
c  Add tau_ij's for RSMs:
      if (ivisc(1) .ge. 70) then
        re=reue/xmach
        do 8296 kpl=1,npl
          kk=k+kpl-1
          jiv2=(kpl-1)*jdim1
          do 8295 i=1,idim
            jiv=(i-1)*npl*jdim1 + (kpl-1)*jdim1
            if (i.eq.1) then
              do 8294 j=1,jdim1
              rho=0.5*(q(j,kk,i,1)+qi0(j,kk,1,1))*(1.-wi0(j+jiv2,5))+
     +            qi0(j,kk,1,1)*wi0(j+jiv2,5)
              u=0.5*(q(j,kk,i,2)+qi0(j,kk,2,1))*(1.-wi0(j+jiv2,5))+
     +            qi0(j,kk,2,1)*wi0(j+jiv2,5)
              v=0.5*(q(j,kk,i,3)+qi0(j,kk,3,1))*(1.-wi0(j+jiv2,5))+
     +            qi0(j,kk,3,1)*wi0(j+jiv2,5)
              w=0.5*(q(j,kk,i,4)+qi0(j,kk,4,1))*(1.-wi0(j+jiv2,5))+
     +            qi0(j,kk,4,1)*wi0(j+jiv2,5)
              t11 =-rho*0.5*(zksav(j,kk,i,1)+ti0(j,kk,1,1))*re
              t22 =-rho*0.5*(zksav(j,kk,i,2)+ti0(j,kk,2,1))*re
              t33 =-rho*0.5*(zksav(j,kk,i,3)+ti0(j,kk,3,1))*re
              t12 =-rho*0.5*(zksav(j,kk,i,4)+ti0(j,kk,4,1))*re
              t23 =-rho*0.5*(zksav(j,kk,i,5)+ti0(j,kk,5,1))*re
              t13 =-rho*0.5*(zksav(j,kk,i,6)+ti0(j,kk,6,1))*re
              t(j+jiv,28)=xmre*(t(j+jiv,25)*t11+t(j+jiv,26)*t12+
     +                          t(j+jiv,27)*t13)/t(j+jiv,15)
              t(j+jiv,29)=xmre*(t(j+jiv,25)*t12+t(j+jiv,26)*t22+
     +                          t(j+jiv,27)*t23)/t(j+jiv,15)
              t(j+jiv,30)=xmre*(t(j+jiv,25)*t13+t(j+jiv,26)*t23+
     +                          t(j+jiv,27)*t33)/t(j+jiv,15)
              t(j+jiv,31)=xmre/t(j+jiv,15)*(
     +           t(j+jiv,25)*(u*t11 + v*t12 + w*t13)
     +          +t(j+jiv,26)*(u*t12 + v*t22 + w*t23)
     +          +t(j+jiv,27)*(u*t13 + v*t23 + w*t33))
 8294         continue
            else if(i.eq.idim) then
              do 8293 j=1,jdim1
              rho=0.5*(q(j,kk,i-1,1)+qi0(j,kk,1,3))*(1.-wi0(j+jiv2,15))+
     +            qi0(j,kk,1,3)*wi0(j+jiv2,15)
              u=0.5*(q(j,kk,i-1,2)+qi0(j,kk,2,3))*(1.-wi0(j+jiv2,15))+
     +            qi0(j,kk,2,3)*wi0(j+jiv2,15)
              v=0.5*(q(j,kk,i-1,3)+qi0(j,kk,3,3))*(1.-wi0(j+jiv2,15))+
     +            qi0(j,kk,3,3)*wi0(j+jiv2,15)
              w=0.5*(q(j,kk,i-1,4)+qi0(j,kk,4,3))*(1.-wi0(j+jiv2,15))+
     +            qi0(j,kk,4,3)*wi0(j+jiv2,15)
              t11 =-rho*0.5*(zksav(j,kk,i-1,1)+ti0(j,kk,1,3))*re
              t22 =-rho*0.5*(zksav(j,kk,i-1,2)+ti0(j,kk,2,3))*re
              t33 =-rho*0.5*(zksav(j,kk,i-1,3)+ti0(j,kk,3,3))*re
              t12 =-rho*0.5*(zksav(j,kk,i-1,4)+ti0(j,kk,4,3))*re
              t23 =-rho*0.5*(zksav(j,kk,i-1,5)+ti0(j,kk,5,3))*re
              t13 =-rho*0.5*(zksav(j,kk,i-1,6)+ti0(j,kk,6,3))*re
              t(j+jiv,28)=xmre*(t(j+jiv,25)*t11+t(j+jiv,26)*t12+
     +                          t(j+jiv,27)*t13)/t(j+jiv,15)
              t(j+jiv,29)=xmre*(t(j+jiv,25)*t12+t(j+jiv,26)*t22+
     +                          t(j+jiv,27)*t23)/t(j+jiv,15)
              t(j+jiv,30)=xmre*(t(j+jiv,25)*t13+t(j+jiv,26)*t23+
     +                          t(j+jiv,27)*t33)/t(j+jiv,15)
              t(j+jiv,31)=xmre/t(j+jiv,15)*(
     +           t(j+jiv,25)*(u*t11 + v*t12 + w*t13)
     +          +t(j+jiv,26)*(u*t12 + v*t22 + w*t23)
     +          +t(j+jiv,27)*(u*t13 + v*t23 + w*t33))
 8293         continue
            else
              do 8292 j=1,jdim1
              rho=0.5*(q(j,kk,i,1)+q(j,kk,i-1,1))
              t11=-rho*0.5*(zksav(j,kk,i,1)+zksav(j,kk,i-1,1))*re
              t22=-rho*0.5*(zksav(j,kk,i,2)+zksav(j,kk,i-1,2))*re
              t33=-rho*0.5*(zksav(j,kk,i,3)+zksav(j,kk,i-1,3))*re
              t12=-rho*0.5*(zksav(j,kk,i,4)+zksav(j,kk,i-1,4))*re
              t23=-rho*0.5*(zksav(j,kk,i,5)+zksav(j,kk,i-1,5))*re
              t13=-rho*0.5*(zksav(j,kk,i,6)+zksav(j,kk,i-1,6))*re
              t(j+jiv,28)=xmre*(t(j+jiv,25)*t11+t(j+jiv,26)*t12+
     +                          t(j+jiv,27)*t13)/t(j+jiv,15)
              t(j+jiv,29)=xmre*(t(j+jiv,25)*t12+t(j+jiv,26)*t22+
     +                          t(j+jiv,27)*t23)/t(j+jiv,15)
              t(j+jiv,30)=xmre*(t(j+jiv,25)*t13+t(j+jiv,26)*t23+
     +                          t(j+jiv,27)*t33)/t(j+jiv,15)
              t(j+jiv,31)=xmre/t(j+jiv,15)*(
     +           t(j+jiv,25)*(0.5*(q(j,kk,i,2)+q(j,kk,i-1,2))*t11
     +                       +0.5*(q(j,kk,i,3)+q(j,kk,i-1,3))*t12
     +                       +0.5*(q(j,kk,i,4)+q(j,kk,i-1,4))*t13)
     +          +t(j+jiv,26)*(0.5*(q(j,kk,i,2)+q(j,kk,i-1,2))*t12
     +                       +0.5*(q(j,kk,i,3)+q(j,kk,i-1,3))*t22
     +                       +0.5*(q(j,kk,i,4)+q(j,kk,i-1,4))*t23)
     +          +t(j+jiv,27)*(0.5*(q(j,kk,i,2)+q(j,kk,i-1,2))*t13
     +                       +0.5*(q(j,kk,i,3)+q(j,kk,i-1,3))*t23
     +                       +0.5*(q(j,kk,i,4)+q(j,kk,i-1,4))*t33))
 8292         continue
            end if
 8295     continue
 8296   continue
        do 8271 l=2,5
        do 8270 izz=1,n
c   Take positive derivative, because want negative of all above for correct
c   sign (this term is going into -d/dxj(Hv) term)
          t(izz,l)=t(izz+js-1,l+26)-t(izz,l+26)
 8270   continue
 8271   continue
        do 8268 kpl=1,npl
          kk=k+kpl-1
          do 8268 i=1,idim1
            jiv=(i-1)*jdim1*npl + (kpl-1)*jdim1+1
            do 8268 l=2,5
              do 8267 izz=1,jdim1
                res(izz,kk,i,l)=res(izz,kk,i,l)+t(izz+jiv-1,l)
 8267         continue
 8268   continue
      end if
c
      if(i_tauijs .ne. 1 .and. (ivisc(1) .eq. 11 .or. ivisc(1) .eq. 12)
     .    .and. level .ge. lglobal) then
c f7=factor for determining whether it's a k-omega or k-epsilon formulation
c   f7=0 for k-epsilon, 1 for k-omega
        f7=0.
        if (ivisc(1) .eq. 12) f7=1.
        re=reue/xmach
        gg=1./(c10+c5-1.)
c   Add nonlinear terms to RHS
        do 1276 kpl=1,npl
          kk=k+kpl-1
          jiv2=(kpl-1)*jdim1
          do 1275 i=1,idim
            jiv=(i-1)*npl*jdim1 + (kpl-1)*jdim1
            if (i.eq.1) then
              do 1274 j=1,jdim1
c            Determine Sij and Wij values (at cell interface)
              s11 = t(j+jiv,8) *t(j+jiv,25)
              s22 = t(j+jiv,9) *t(j+jiv,26)
              s33 = t(j+jiv,10)*t(j+jiv,27)
              s12 = 0.5*(t(j+jiv,8) *t(j+jiv,26) +
     +                   t(j+jiv,9) *t(j+jiv,25))
              s13 = 0.5*(t(j+jiv,8) *t(j+jiv,27) +
     +                   t(j+jiv,10)*t(j+jiv,25))
              s23 = 0.5*(t(j+jiv,9) *t(j+jiv,27) +
     +                   t(j+jiv,10)*t(j+jiv,26))
              w12 = 0.5*(t(j+jiv,8) *t(j+jiv,26) -
     +                   t(j+jiv,9) *t(j+jiv,25))
              w13 = 0.5*(t(j+jiv,8) *t(j+jiv,27) -
     +                   t(j+jiv,10)*t(j+jiv,25))
              w23 = 0.5*(t(j+jiv,9) *t(j+jiv,27) -
     +                   t(j+jiv,10)*t(j+jiv,26))
              xis = s11*s11 + s22*s22 + s33*s33 +
     +              2.*s12*s12 + 2.*s13*s13 + 2.*s23*s23
              wis = 2.*w12*w12 + 2.*w13*w13 + 2.*w23*w23
c            "Constants" are function of c3, c4, and gg:
              alpa1=(2.-c4)/2.*gg
              alpa2=(2.-c3)*gg
c            Find tauij values:
              rho=0.5*(q(j,kk,i,1)+qi0(j,kk,1,1))*(1.-wi0(j+jiv2,5))+
     +            qi0(j,kk,1,1)*wi0(j+jiv2,5)
              u=0.5*(q(j,kk,i,2)+qi0(j,kk,2,1))*(1.-wi0(j+jiv2,5))+
     +            qi0(j,kk,2,1)*wi0(j+jiv2,5)
              v=0.5*(q(j,kk,i,3)+qi0(j,kk,3,1))*(1.-wi0(j+jiv2,5))+
     +            qi0(j,kk,3,1)*wi0(j+jiv2,5)
              w=0.5*(q(j,kk,i,4)+qi0(j,kk,4,1))*(1.-wi0(j+jiv2,5))+
     +            qi0(j,kk,4,1)*wi0(j+jiv2,5)
              zksav1=0.5*(zksav(j,kk,i,1)+ti0(j,kk,1,1))
              zksav1=ccmax(zksav1,tur10(1))
              zksav2=0.5*(zksav(j,kk,i,2)+ti0(j,kk,2,1))
              bnum=zksav2*(1.-f7)+f7
              eta=(2.-c3)**2*(gg*gg/4.)*xis*(bnum/(zksav1*re))**2
              squig=(2.-c4)**2*(gg*gg/4.)*wis*(bnum/(zksav1*re))**2
              eta=ccmincr(eta,10.)
              squig=ccmincr(squig,10.)
              factre=(3.*(1.+eta)+.2e-8*(eta*eta*eta +
     +          squig*squig*squig))/
     +               (3.*(1.+eta)+   .2*(eta*eta*eta +
     +          squig*squig*squig))
              t11 = 2.*rho*zksav2*re/3.
     +             -4.*alpa1*t(j+jiv,24)*factre*bnum/
     +                 zksav1*(-s12*w12 - s13*w13)/re
     +           +2.*alpa2*t(j+jiv,24)*factre*bnum/zksav1*
     +                 (s11*s11 + s12*s12 + s13*s13 - 0.33333*xis)/re
              xlnpt=2.*t(j+jiv,24)*(s11-(s11+s22+s33)/3.)
              t11=ccmax(t11,xlnpt)
              t11 = t11*(1.-bci(j,kk,1))
              t22 = 2.*rho*zksav2*re/3.
     +             -4.*alpa1*t(j+jiv,24)*factre*bnum/
     +                 zksav1*( s12*w12 - s23*w23)/re
     +           +2.*alpa2*t(j+jiv,24)*factre*bnum/zksav1*
     +                 (s22*s22 + s12*s12 + s23*s23 - 0.33333*xis)/re
              xlnpt=2.*t(j+jiv,24)*(s22-(s11+s22+s33)/3.)
              t22=ccmax(t22,xlnpt)
              t22 = t22*(1.-bci(j,kk,1))
              t33 = 2.*rho*zksav2*re/3.
     +             -4.*alpa1*t(j+jiv,24)*factre*bnum/
     +                 zksav1*( s23*w23 + s13*w13)/re
     +           +2.*alpa2*t(j+jiv,24)*factre*bnum/zksav1*
     +                 (s33*s33 + s23*s23 + s13*s13 - 0.33333*xis)/re
              xlnpt=2.*t(j+jiv,24)*(s33-(s11+s22+s33)/3.)
              t33=ccmax(t33,xlnpt)
              t33 = t33*(1.-bci(j,kk,1))
              t12 =-2.*alpa1*t(j+jiv,24)*factre*bnum/zksav1*
     +                 ( s11*w12 - s22*w12 - s13*w23 - s23*w13)/re
     +           +2.*alpa2*t(j+jiv,24)*factre*bnum/zksav1*
     +                 (s11*s12 + s12*s22 + s13*s23)/re
              t12 = t12*(1.-bci(j,kk,1))
              t13 =-2.*alpa1*t(j+jiv,24)*factre*bnum/zksav1*
     +                 ( s11*w13 + s12*w23 - s23*w12 - s33*w13)/re
     +           +2.*alpa2*t(j+jiv,24)*factre*bnum/zksav1*
     +                 (s11*s13 + s12*s23 + s13*s33)/re
              t13 = t13*(1.-bci(j,kk,1))
              t23 =-2.*alpa1*t(j+jiv,24)*factre*bnum/zksav1*
     +                 ( s12*w13 + s22*w23 + s13*w12 - s33*w23)/re
     +           +2.*alpa2*t(j+jiv,24)*factre*bnum/zksav1*
     +                 (s12*s13 + s22*s23 + s23*s33)/re
              t23 = t23*(1.-bci(j,kk,1))
c
              t(j+jiv,28)=xmre*(t(j+jiv,25)*t11+t(j+jiv,26)*t12+
     +                          t(j+jiv,27)*t13)/t(j+jiv,15)
              t(j+jiv,29)=xmre*(t(j+jiv,25)*t12+t(j+jiv,26)*t22+
     +                          t(j+jiv,27)*t23)/t(j+jiv,15)
              t(j+jiv,30)=xmre*(t(j+jiv,25)*t13+t(j+jiv,26)*t23+
     +                          t(j+jiv,27)*t33)/t(j+jiv,15)
              t(j+jiv,31)=xmre/t(j+jiv,15)*(
     +           t(j+jiv,25)*(u*t11 + v*t12 + w*t13)
     +          +t(j+jiv,26)*(u*t12 + v*t22 + w*t23)
     +          +t(j+jiv,27)*(u*t13 + v*t23 + w*t33))
c   to get mulam, need mulam=t14*t15/xmre-t24
c   (currently uses eddy viscosity mu_t (t(24)), even though diffusion term 
c   in k-eqn uses cmuc1*rho*k/omega or cmuc1*rho*k^2/epsilon)
     +         -xmre/t(j+jiv,15)*(t(j+jiv,14)*
     +          t(j+jiv,15)/xmre-t(j+jiv,24)+sigk1*t(j+jiv,24))*
     +          (t(j+jiv,25)**2 + t(j+jiv,26)**2 + t(j+jiv,27)**2)*
     +          (zksav(j,kk,i,2)-ti0(j,kk,2,1))
 1274         continue
            else if(i.eq.idim) then
              do 1273 j=1,jdim1
c            Determine Sij and Wij values (at cell interface)
              s11 = t(j+jiv,8) *t(j+jiv,25)
              s22 = t(j+jiv,9) *t(j+jiv,26)
              s33 = t(j+jiv,10)*t(j+jiv,27)
              s12 = 0.5*(t(j+jiv,8) *t(j+jiv,26) +
     +                   t(j+jiv,9) *t(j+jiv,25))
              s13 = 0.5*(t(j+jiv,8) *t(j+jiv,27) +
     +                   t(j+jiv,10)*t(j+jiv,25))
              s23 = 0.5*(t(j+jiv,9) *t(j+jiv,27) +
     +                   t(j+jiv,10)*t(j+jiv,26))
              w12 = 0.5*(t(j+jiv,8) *t(j+jiv,26) -
     +                   t(j+jiv,9) *t(j+jiv,25))
              w13 = 0.5*(t(j+jiv,8) *t(j+jiv,27) -
     +                   t(j+jiv,10)*t(j+jiv,25))
              w23 = 0.5*(t(j+jiv,9) *t(j+jiv,27) -
     +                   t(j+jiv,10)*t(j+jiv,26))
              xis = s11*s11 + s22*s22 + s33*s33 +
     +              2.*s12*s12 + 2.*s13*s13 + 2.*s23*s23
              wis = 2.*w12*w12 + 2.*w13*w13 + 2.*w23*w23
c            "Constants" are function of c3, c4, and gg:
              alpa1=(2.-c4)/2.*gg
              alpa2=(2.-c3)*gg
c            Find tauij values:
              rho=0.5*(q(j,kk,i-1,1)+qi0(j,kk,1,3))*(1.-wi0(j+jiv2,15))+
     +            qi0(j,kk,1,3)*wi0(j+jiv2,15)
              u=0.5*(q(j,kk,i-1,2)+qi0(j,kk,2,3))*(1.-wi0(j+jiv2,15))+
     +            qi0(j,kk,2,3)*wi0(j+jiv2,15)
              v=0.5*(q(j,kk,i-1,3)+qi0(j,kk,3,3))*(1.-wi0(j+jiv2,15))+
     +            qi0(j,kk,3,3)*wi0(j+jiv2,15)
              w=0.5*(q(j,kk,i-1,4)+qi0(j,kk,4,3))*(1.-wi0(j+jiv2,15))+
     +            qi0(j,kk,4,3)*wi0(j+jiv2,15)
              zksav1=0.5*(zksav(j,kk,i-1,1)+ti0(j,kk,1,3))
              zksav1=ccmax(zksav1,tur10(1))
              zksav2=0.5*(zksav(j,kk,i-1,2)+ti0(j,kk,2,3))
              bnum=zksav2*(1.-f7)+f7
              eta=(2.-c3)**2*(gg*gg/4.)*xis*(bnum/(zksav1*re))**2
              squig=(2.-c4)**2*(gg*gg/4.)*wis*(bnum/(zksav1*re))**2
              eta=ccmincr(eta,10.)
              squig=ccmincr(squig,10.)
              factre=(3.*(1.+eta)+.2e-8*(eta*eta*eta +
     +          squig*squig*squig))/
     +               (3.*(1.+eta)+   .2*(eta*eta*eta +
     +          squig*squig*squig))
              t11 = 2.*rho*zksav2*re/3.
     +             -4.*alpa1*t(j+jiv,24)*factre*bnum/
     +                 zksav1*(-s12*w12 - s13*w13)/re
     +           +2.*alpa2*t(j+jiv,24)*factre*bnum/zksav1*
     +                 (s11*s11 + s12*s12 + s13*s13 - 0.33333*xis)/re
              xlnpt=2.*t(j+jiv,24)*(s11-(s11+s22+s33)/3.)
              t11=ccmax(t11,xlnpt)
              t11 = t11*(1.-bci(j,kk,2))
              t22 = 2.*rho*zksav2*re/3.
     +             -4.*alpa1*t(j+jiv,24)*factre*bnum/
     +                 zksav1*( s12*w12 - s23*w23)/re
     +           +2.*alpa2*t(j+jiv,24)*factre*bnum/zksav1*
     +                 (s22*s22 + s12*s12 + s23*s23 - 0.33333*xis)/re
              xlnpt=2.*t(j+jiv,24)*(s22-(s11+s22+s33)/3.)
              t22=ccmax(t22,xlnpt)
              t22 = t22*(1.-bci(j,kk,2))
              t33 = 2.*rho*zksav2*re/3.
     +             -4.*alpa1*t(j+jiv,24)*factre*bnum/
     +                 zksav1*( s23*w23 + s13*w13)/re
     +           +2.*alpa2*t(j+jiv,24)*factre*bnum/zksav1*
     +                 (s33*s33 + s23*s23 + s13*s13 - 0.33333*xis)/re
              xlnpt=2.*t(j+jiv,24)*(s33-(s11+s22+s33)/3.)
              t33=ccmax(t33,xlnpt)
              t33 = t33*(1.-bci(j,kk,2))
              t12 =-2.*alpa1*t(j+jiv,24)*factre*bnum/zksav1*
     +                 ( s11*w12 - s22*w12 - s13*w23 - s23*w13)/re
     +           +2.*alpa2*t(j+jiv,24)*factre*bnum/zksav1*
     +                 (s11*s12 + s12*s22 + s13*s23)/re
              t12 = t12*(1.-bci(j,kk,2))
              t13 =-2.*alpa1*t(j+jiv,24)*factre*bnum/zksav1*
     +                 ( s11*w13 + s12*w23 - s23*w12 - s33*w13)/re
     +           +2.*alpa2*t(j+jiv,24)*factre*bnum/zksav1*
     +                 (s11*s13 + s12*s23 + s13*s33)/re
              t13 = t13*(1.-bci(j,kk,2))
              t23 =-2.*alpa1*t(j+jiv,24)*factre*bnum/zksav1*
     +                 ( s12*w13 + s22*w23 + s13*w12 - s33*w23)/re
     +           +2.*alpa2*t(j+jiv,24)*factre*bnum/zksav1*
     +                 (s12*s13 + s22*s23 + s23*s33)/re
              t23 = t23*(1.-bci(j,kk,2))
c
              t(j+jiv,28)=xmre*(t(j+jiv,25)*t11+t(j+jiv,26)*t12+
     +                          t(j+jiv,27)*t13)/t(j+jiv,15)
              t(j+jiv,29)=xmre*(t(j+jiv,25)*t12+t(j+jiv,26)*t22+
     +                          t(j+jiv,27)*t23)/t(j+jiv,15)
              t(j+jiv,30)=xmre*(t(j+jiv,25)*t13+t(j+jiv,26)*t23+
     +                          t(j+jiv,27)*t33)/t(j+jiv,15)
              t(j+jiv,31)=xmre/t(j+jiv,15)*(
     +           t(j+jiv,25)*(u*t11 + v*t12 + w*t13)
     +          +t(j+jiv,26)*(u*t12 + v*t22 + w*t23)
     +          +t(j+jiv,27)*(u*t13 + v*t23 + w*t33))
c   to get mulam, need mulam=t14*t15/xmre-t24
c   (currently uses eddy viscosity mu_t (t(24)), even though diffusion term 
c   in k-eqn uses cmuc1*rho*k/omega or cmuc1*rho*k^2/epsilon)
     +         -xmre/t(j+jiv,15)*(t(j+jiv,14)*
     +          t(j+jiv,15)/xmre-t(j+jiv,24)+sigk1*t(j+jiv,24))*
     +          (t(j+jiv,25)**2 + t(j+jiv,26)**2 + t(j+jiv,27)**2)*
     +          (ti0(j,kk,2,3)-zksav(j,kk,i-1,2))
 1273         continue
            else
              do 1272 j=1,jdim1
c            Determine Sij and Wij values (at cell interfaces):
              s11 = t(j+jiv,8) *t(j+jiv,25)
              s22 = t(j+jiv,9) *t(j+jiv,26)
              s33 = t(j+jiv,10)*t(j+jiv,27)
              s12 = 0.5*(t(j+jiv,8) *t(j+jiv,26) +
     +                   t(j+jiv,9) *t(j+jiv,25))
              s13 = 0.5*(t(j+jiv,8) *t(j+jiv,27) +
     +                   t(j+jiv,10)*t(j+jiv,25))
              s23 = 0.5*(t(j+jiv,9) *t(j+jiv,27) +
     +                   t(j+jiv,10)*t(j+jiv,26))
              w12 = 0.5*(t(j+jiv,8) *t(j+jiv,26) -
     +                   t(j+jiv,9) *t(j+jiv,25))
              w13 = 0.5*(t(j+jiv,8) *t(j+jiv,27) -
     +                   t(j+jiv,10)*t(j+jiv,25))
              w23 = 0.5*(t(j+jiv,9) *t(j+jiv,27) -
     +                   t(j+jiv,10)*t(j+jiv,26))
              xis = s11*s11 + s22*s22 + s33*s33 +
     +              2.*s12*s12 + 2.*s13*s13 + 2.*s23*s23
              wis = 2.*w12*w12 + 2.*w13*w13 + 2.*w23*w23
c            "Constants" are function of c3, c4, and gg:
              alpa1=(2.-c4)/2.*gg
              alpa2=(2.-c3)*gg
c            Find tauij values:
              rho=0.5*(q(j,kk,i,1)+q(j,kk,i-1,1))
              zksav1=0.5*(zksav(j,kk,i,1)+zksav(j,kk,i-1,1))
              zksav1=ccmax(zksav1,tur10(1))
              zksav2=0.5*(zksav(j,kk,i,2)+zksav(j,kk,i-1,2))
              bnum=zksav2*(1.-f7)+f7
              eta=(2.-c3)**2*(gg*gg/4.)*xis*(bnum/(zksav1*re))**2
              squig=(2.-c4)**2*(gg*gg/4.)*wis*(bnum/(zksav1*re))**2
              eta=ccmincr(eta,10.)
              squig=ccmincr(squig,10.)
              factre=(3.*(1.+eta)+.2e-8*(eta*eta*eta +
     +          squig*squig*squig))/
     +               (3.*(1.+eta)+   .2*(eta*eta*eta +
     +          squig*squig*squig))
              t11 = 2.*rho*zksav2*re/3.
     +             -4.*alpa1*t(j+jiv,24)*factre*bnum/
     +                 zksav1*(-s12*w12 - s13*w13)/re
     +           +2.*alpa2*t(j+jiv,24)*factre*bnum/zksav1*
     +                 (s11*s11 + s12*s12 + s13*s13 - 0.33333*xis)/re
              xlnpt=2.*t(j+jiv,24)*(s11-(s11+s22+s33)/3.)
              t11=ccmax(t11,xlnpt)
              t22 = 2.*rho*zksav2*re/3.
     +             -4.*alpa1*t(j+jiv,24)*factre*bnum/
     +                 zksav1*( s12*w12 - s23*w23)/re
     +           +2.*alpa2*t(j+jiv,24)*factre*bnum/zksav1*
     +                 (s22*s22 + s12*s12 + s23*s23 - 0.33333*xis)/re
              xlnpt=2.*t(j+jiv,24)*(s22-(s11+s22+s33)/3.)
              t22=ccmax(t22,xlnpt)
              t33 = 2.*rho*zksav2*re/3.
     +             -4.*alpa1*t(j+jiv,24)*factre*bnum/
     +                 zksav1*( s23*w23 + s13*w13)/re
     +           +2.*alpa2*t(j+jiv,24)*factre*bnum/zksav1*
     +                 (s33*s33 + s23*s23 + s13*s13 - 0.33333*xis)/re
              xlnpt=2.*t(j+jiv,24)*(s33-(s11+s22+s33)/3.)
              t33=ccmax(t33,xlnpt)
              t12 =-2.*alpa1*t(j+jiv,24)*factre*bnum/zksav1*
     +                 ( s11*w12 - s22*w12 - s13*w23 - s23*w13)/re
     +           +2.*alpa2*t(j+jiv,24)*factre*bnum/zksav1*
     +                 (s11*s12 + s12*s22 + s13*s23)/re
              t13 =-2.*alpa1*t(j+jiv,24)*factre*bnum/zksav1*
     +                 ( s11*w13 + s12*w23 - s23*w12 - s33*w13)/re
     +           +2.*alpa2*t(j+jiv,24)*factre*bnum/zksav1*
     +                 (s11*s13 + s12*s23 + s13*s33)/re
              t23 =-2.*alpa1*t(j+jiv,24)*factre*bnum/zksav1*
     +                 ( s12*w13 + s22*w23 + s13*w12 - s33*w23)/re
     +           +2.*alpa2*t(j+jiv,24)*factre*bnum/zksav1*
     +                 (s12*s13 + s22*s23 + s23*s33)/re
c
              t(j+jiv,28)=xmre*(t(j+jiv,25)*t11+t(j+jiv,26)*t12+
     +                          t(j+jiv,27)*t13)/t(j+jiv,15)
              t(j+jiv,29)=xmre*(t(j+jiv,25)*t12+t(j+jiv,26)*t22+
     +                          t(j+jiv,27)*t23)/t(j+jiv,15)
              t(j+jiv,30)=xmre*(t(j+jiv,25)*t13+t(j+jiv,26)*t23+
     +                          t(j+jiv,27)*t33)/t(j+jiv,15)
              t(j+jiv,31)=xmre/t(j+jiv,15)*(
     +           t(j+jiv,25)*(0.5*(q(j,kk,i,2)+q(j,kk,i-1,2))*t11
     +                       +0.5*(q(j,kk,i,3)+q(j,kk,i-1,3))*t12
     +                       +0.5*(q(j,kk,i,4)+q(j,kk,i-1,4))*t13)
     +          +t(j+jiv,26)*(0.5*(q(j,kk,i,2)+q(j,kk,i-1,2))*t12
     +                       +0.5*(q(j,kk,i,3)+q(j,kk,i-1,3))*t22
     +                       +0.5*(q(j,kk,i,4)+q(j,kk,i-1,4))*t23)
     +          +t(j+jiv,27)*(0.5*(q(j,kk,i,2)+q(j,kk,i-1,2))*t13
     +                       +0.5*(q(j,kk,i,3)+q(j,kk,i-1,3))*t23
     +                       +0.5*(q(j,kk,i,4)+q(j,kk,i-1,4))*t33))
c   to get mulam, need mulam=t14*t15/xmre-t24
c   (currently uses eddy viscosity mu_t (t(24)), even though diffusion term 
c   in k-eqn uses cmuc1*rho*k/omega or cmuc1*rho*k^2/epsilon)
     +         -xmre/t(j+jiv,15)*(t(j+jiv,14)*
     +          t(j+jiv,15)/xmre-t(j+jiv,24)+sigk1*t(j+jiv,24))*
     +          (t(j+jiv,25)**2 + t(j+jiv,26)**2 + t(j+jiv,27)**2)*
     +          (zksav(j,kk,i,2)-zksav(j,kk,i-1,2))
 1272         continue
            end if
 1275     continue
 1276   continue
      else if(i_tauijs .ne. 1 .and. (ivisc(1) .eq. 13 .or.
     +  ivisc(1) .eq. 14) .and. level .ge. lglobal) then
c f7=factor for determining whether it's a k-omega or k-epsilon formulation
c   f7=0 for k-epsilon, 1 for k-omega
        f7=0.
        if (ivisc(1) .eq. 14) f7=1.
        re=reue/xmach
              al10 = 0.5*c10-1.
              al1 = 2.*(0.5*c11+1.)
c          2-line improvement for wake
              if (ieasm_type .eq. 0 .or. ieasm_type .eq. 3 .or.
     +            ieasm_type .eq. 4) then
                al10=al10+1.8864
                al1=al1-2.
              end if
              al2 = 0.5*c2 -2./3.
              al3 = 0.5*c3 -1.
              al4 = 0.5*c4 -1.
c   Add nonlinear terms to RHS
        do 1286 kpl=1,npl
          kk=k+kpl-1
          jiv2=(kpl-1)*jdim1
          do 1285 i=1,idim
            jiv=(i-1)*npl*jdim1 + (kpl-1)*jdim1
            if (i.eq.1) then
              do 1284 j=1,jdim1
c            Determine Sij and Wij values (at cell interface)
              s11 = t(j+jiv,8) *t(j+jiv,25)
              s22 = t(j+jiv,9) *t(j+jiv,26)
              s33 = t(j+jiv,10)*t(j+jiv,27)
              s12 = 0.5*(t(j+jiv,8) *t(j+jiv,26) +
     +                   t(j+jiv,9) *t(j+jiv,25))
              s13 = 0.5*(t(j+jiv,8) *t(j+jiv,27) +
     +                   t(j+jiv,10)*t(j+jiv,25))
              s23 = 0.5*(t(j+jiv,9) *t(j+jiv,27) +
     +                   t(j+jiv,10)*t(j+jiv,26))
              w12 = 0.5*(t(j+jiv,8) *t(j+jiv,26) -
     +                   t(j+jiv,9) *t(j+jiv,25))
              w13 = 0.5*(t(j+jiv,8) *t(j+jiv,27) -
     +                   t(j+jiv,10)*t(j+jiv,25))
              w23 = 0.5*(t(j+jiv,9) *t(j+jiv,27) -
     +                   t(j+jiv,10)*t(j+jiv,26))
              xis = s11*s11 + s22*s22 + s33*s33 +
     +              2.*s12*s12 + 2.*s13*s13 + 2.*s23*s23
c             Girimaji JFM 2000 fix to c4
              if (ieasm_type .eq. 4) then
                wis = 2.*w12*w12 + 2.*w13*w13 + 2.*w23*w23
                eta1_girimaji=xis/(xis+wis)
                if (real(eta1_girimaji) .lt. 0.5) then
                  c4new=2.0-((2.0-c4)*
     +                  (eta1_girimaji/(1.-eta1_girimaji))**0.75)
                else
                  c4new=c4
                end if
                al4 = 0.5*c4new -1.
              end if
c            "Constants" are function of c3, c4, and gg:
c            Find tauij values:
              rho=0.5*(q(j,kk,i,1)+qi0(j,kk,1,1))*(1.-wi0(j+jiv2,5))+
     +            qi0(j,kk,1,1)*wi0(j+jiv2,5)
              u=0.5*(q(j,kk,i,2)+qi0(j,kk,2,1))*(1.-wi0(j+jiv2,5))+
     +            qi0(j,kk,2,1)*wi0(j+jiv2,5)
              v=0.5*(q(j,kk,i,3)+qi0(j,kk,3,1))*(1.-wi0(j+jiv2,5))+
     +            qi0(j,kk,3,1)*wi0(j+jiv2,5)
              w=0.5*(q(j,kk,i,4)+qi0(j,kk,4,1))*(1.-wi0(j+jiv2,5))+
     +            qi0(j,kk,4,1)*wi0(j+jiv2,5)
              zksav1=0.5*(zksav(j,kk,i,1)+ti0(j,kk,1,1))
              zksav2=0.5*(zksav(j,kk,i,2)+ti0(j,kk,2,1))
              cmuu = cmuv(j,kk,i)
              bnum=zksav2*(1.-f7)+f7
c             Durbin TCFD 1991 near-wall limiter
              if (idurbinlim .ne. 0 .and. (ieasm_type .eq. 3 .or.
     +            ieasm_type .eq. 4)) then
                tau=bnum/zksav1
                fmu=t(j+jiv,14)*t(j+jiv,15)/xmre-t(j+jiv,24)
                zkolmog=6.*sqrt(fmu*bnum/(rho*zksav2*zksav1))
                tau=ccmax(tau,zkolmog)
                zksav1=bnum/tau
              end if
              zksav1=ccmax(zksav1,tur10(1))
              eta1=xis*(bnum/zksav1/re)**2
              alpa1 = -al4/(al10-eta1*al1*cmuu)
              alpa2 = -2.*al3/(al10-eta1*cmuu*al1)
              t11 = 2.*rho*zksav2*re/3.
     +             -4.*alpa1*t(j+jiv,24)*bnum/
     +                 zksav1*(-s12*w12 - s13*w13)/re
     +           +2.*alpa2*t(j+jiv,24)*bnum/zksav1*
     +                 (s11*s11 + s12*s12 + s13*s13 - 0.33333*xis)/re
              xlnpt=2.*t(j+jiv,24)*(s11-(s11+s22+s33)/3.)
              t11=ccmax(t11,xlnpt)
              t11 = t11*(1.-bci(j,kk,1))
              t22 = 2.*rho*zksav2*re/3.
     +             -4.*alpa1*t(j+jiv,24)*bnum/
     +                 zksav1*( s12*w12 - s23*w23)/re
     +           +2.*alpa2*t(j+jiv,24)*bnum/zksav1*
     +                 (s22*s22 + s12*s12 + s23*s23 - 0.33333*xis)/re
              xlnpt=2.*t(j+jiv,24)*(s22-(s11+s22+s33)/3.)
              t22=ccmax(t22,xlnpt)
              t22 = t22*(1.-bci(j,kk,1))
              t33 = 2.*rho*zksav2*re/3.
     +             -4.*alpa1*t(j+jiv,24)*bnum/
     +                 zksav1*( s23*w23 + s13*w13)/re
     +           +2.*alpa2*t(j+jiv,24)*bnum/zksav1*
     +                 (s33*s33 + s23*s23 + s13*s13 - 0.33333*xis)/re
              xlnpt=2.*t(j+jiv,24)*(s33-(s11+s22+s33)/3.)
              t33=ccmax(t33,xlnpt)
              t33 = t33*(1.-bci(j,kk,1))
              t12 =-2.*alpa1*t(j+jiv,24)*bnum/zksav1*
     +                 ( s11*w12 - s22*w12 - s13*w23 - s23*w13)/re
     +           +2.*alpa2*t(j+jiv,24)*bnum/zksav1*
     +                 (s11*s12 + s12*s22 + s13*s23)/re
              t12 = t12*(1.-bci(j,kk,1))
              t13 =-2.*alpa1*t(j+jiv,24)*bnum/zksav1*
     +                 ( s11*w13 + s12*w23 - s23*w12 - s33*w13)/re
     +           +2.*alpa2*t(j+jiv,24)*bnum/zksav1*
     +                 (s11*s13 + s12*s23 + s13*s33)/re
              t13 = t13*(1.-bci(j,kk,1))
              t23 =-2.*alpa1*t(j+jiv,24)*bnum/zksav1*
     +                 ( s12*w13 + s22*w23 + s13*w12 - s33*w23)/re
     +           +2.*alpa2*t(j+jiv,24)*bnum/zksav1*
     +                 (s12*s13 + s22*s23 + s23*s33)/re
              t23 = t23*(1.-bci(j,kk,1))
c
              t(j+jiv,28)=xmre*(t(j+jiv,25)*t11+t(j+jiv,26)*t12+
     +                          t(j+jiv,27)*t13)/t(j+jiv,15)
              t(j+jiv,29)=xmre*(t(j+jiv,25)*t12+t(j+jiv,26)*t22+
     +                          t(j+jiv,27)*t23)/t(j+jiv,15)
              t(j+jiv,30)=xmre*(t(j+jiv,25)*t13+t(j+jiv,26)*t23+
     +                          t(j+jiv,27)*t33)/t(j+jiv,15)
              xmut=f7*t(j+jiv,24)+(1.-f7)*cmuc1*rho*zksav2**2/zksav1
              t(j+jiv,31)=xmre/t(j+jiv,15)*(
     +           t(j+jiv,25)*(u*t11 + v*t12 + w*t13)
     +          +t(j+jiv,26)*(u*t12 + v*t22 + w*t23)
     +          +t(j+jiv,27)*(u*t13 + v*t23 + w*t33))
c   to get mulam, need mulam=t14*t15/xmre-t24
     +         -xmre/t(j+jiv,15)*(t(j+jiv,14)*
     +          t(j+jiv,15)/xmre-t(j+jiv,24)+sigk1*xmut)*
     +          (t(j+jiv,25)**2 + t(j+jiv,26)**2 + t(j+jiv,27)**2)*
     +          (zksav(j,kk,i,2)-ti0(j,kk,2,1))
 1284         continue
            else if(i.eq.idim) then
              do 1283 j=1,jdim1
c            Determine Sij and Wij values (at cell interface)
              s11 = t(j+jiv,8) *t(j+jiv,25)
              s22 = t(j+jiv,9) *t(j+jiv,26)
              s33 = t(j+jiv,10)*t(j+jiv,27)
              s12 = 0.5*(t(j+jiv,8) *t(j+jiv,26) +
     +                   t(j+jiv,9) *t(j+jiv,25))
              s13 = 0.5*(t(j+jiv,8) *t(j+jiv,27) +
     +                   t(j+jiv,10)*t(j+jiv,25))
              s23 = 0.5*(t(j+jiv,9) *t(j+jiv,27) +
     +                   t(j+jiv,10)*t(j+jiv,26))
              w12 = 0.5*(t(j+jiv,8) *t(j+jiv,26) -
     +                   t(j+jiv,9) *t(j+jiv,25))
              w13 = 0.5*(t(j+jiv,8) *t(j+jiv,27) -
     +                   t(j+jiv,10)*t(j+jiv,25))
              w23 = 0.5*(t(j+jiv,9) *t(j+jiv,27) -
     +                   t(j+jiv,10)*t(j+jiv,26))
              xis = s11*s11 + s22*s22 + s33*s33 +
     +              2.*s12*s12 + 2.*s13*s13 + 2.*s23*s23
c             Girimaji JFM 2000 fix to c4
              if (ieasm_type .eq. 4) then
                wis = 2.*w12*w12 + 2.*w13*w13 + 2.*w23*w23
                eta1_girimaji=xis/(xis+wis)
                if (real(eta1_girimaji) .lt. 0.5) then
                  c4new=2.0-((2.0-c4)*
     +                  (eta1_girimaji/(1.-eta1_girimaji))**0.75)
                else
                  c4new=c4
                end if
                al4 = 0.5*c4new -1.
              end if
c            "Constants" are function of c3, c4, and gg:
c            Find tauij values:
              rho=0.5*(q(j,kk,i-1,1)+qi0(j,kk,1,3))*(1.-wi0(j+jiv2,15))+
     +            qi0(j,kk,1,3)*wi0(j+jiv2,15)
              u=0.5*(q(j,kk,i-1,2)+qi0(j,kk,2,3))*(1.-wi0(j+jiv2,15))+
     +            qi0(j,kk,2,3)*wi0(j+jiv2,15)
              v=0.5*(q(j,kk,i-1,3)+qi0(j,kk,3,3))*(1.-wi0(j+jiv2,15))+
     +            qi0(j,kk,3,3)*wi0(j+jiv2,15)
              w=0.5*(q(j,kk,i-1,4)+qi0(j,kk,4,3))*(1.-wi0(j+jiv2,15))+
     +            qi0(j,kk,4,3)*wi0(j+jiv2,15)
              zksav1=0.5*(zksav(j,kk,i-1,1)+ti0(j,kk,1,3))
              zksav2=0.5*(zksav(j,kk,i-1,2)+ti0(j,kk,2,3))
              cmuu=cmuv(j,kk,i-1)
              bnum=zksav2*(1.-f7)+f7
c             Durbin TCFD 1991 near-wall limiter
              if (idurbinlim .ne. 0 .and. (ieasm_type .eq. 3 .or.
     +            ieasm_type .eq. 4)) then
                tau=bnum/zksav1
                fmu=t(j+jiv,14)*t(j+jiv,15)/xmre-t(j+jiv,24)
                zkolmog=6.*sqrt(fmu*bnum/(rho*zksav2*zksav1))
                tau=ccmax(tau,zkolmog)
                zksav1=bnum/tau
              end if
              zksav1=ccmax(zksav1,tur10(1))
              eta1=xis*(bnum/zksav1/re)**2
              alpa1 = -al4/(al10-eta1*al1*cmuu)
              alpa2 = -2.*al3/(al10-eta1*cmuu*al1)
              t11 = 2.*rho*zksav2*re/3.
     +             -4.*alpa1*t(j+jiv,24)*bnum/
     +                 zksav1*(-s12*w12 - s13*w13)/re
     +           +2.*alpa2*t(j+jiv,24)*bnum/zksav1*
     +                 (s11*s11 + s12*s12 + s13*s13 - 0.33333*xis)/re
              xlnpt=2.*t(j+jiv,24)*(s11-(s11+s22+s33)/3.)
              t11=ccmax(t11,xlnpt)
              t11 = t11*(1.-bci(j,kk,2))
              t22 = 2.*rho*zksav2*re/3.
     +             -4.*alpa1*t(j+jiv,24)*bnum/
     +                 zksav1*( s12*w12 - s23*w23)/re
     +           +2.*alpa2*t(j+jiv,24)*bnum/zksav1*
     +                 (s22*s22 + s12*s12 + s23*s23 - 0.33333*xis)/re
              xlnpt=2.*t(j+jiv,24)*(s22-(s11+s22+s33)/3.)
              t22=ccmax(t22,xlnpt)
              t22 = t22*(1.-bci(j,kk,2))
              t33 = 2.*rho*zksav2*re/3.
     +             -4.*alpa1*t(j+jiv,24)*bnum/
     +                 zksav1*( s23*w23 + s13*w13)/re
     +           +2.*alpa2*t(j+jiv,24)*bnum/zksav1*
     +                 (s33*s33 + s23*s23 + s13*s13 - 0.33333*xis)/re
              xlnpt=2.*t(j+jiv,24)*(s33-(s11+s22+s33)/3.)
              t33=ccmax(t33,xlnpt)
              t33 = t33*(1.-bci(j,kk,2))
              t12 =-2.*alpa1*t(j+jiv,24)*bnum/zksav1*
     +                 ( s11*w12 - s22*w12 - s13*w23 - s23*w13)/re
     +           +2.*alpa2*t(j+jiv,24)*bnum/zksav1*
     +                 (s11*s12 + s12*s22 + s13*s23)/re
              t12 = t12*(1.-bci(j,kk,2))
              t13 =-2.*alpa1*t(j+jiv,24)*bnum/zksav1*
     +                 ( s11*w13 + s12*w23 - s23*w12 - s33*w13)/re
     +           +2.*alpa2*t(j+jiv,24)*bnum/zksav1*
     +                 (s11*s13 + s12*s23 + s13*s33)/re
              t13 = t13*(1.-bci(j,kk,2))
              t23 =-2.*alpa1*t(j+jiv,24)*bnum/zksav1*
     +                 ( s12*w13 + s22*w23 + s13*w12 - s33*w23)/re
     +           +2.*alpa2*t(j+jiv,24)*bnum/zksav1*
     +                 (s12*s13 + s22*s23 + s23*s33)/re
              t23 = t23*(1.-bci(j,kk,2))
c
              t(j+jiv,28)=xmre*(t(j+jiv,25)*t11+t(j+jiv,26)*t12+
     +                          t(j+jiv,27)*t13)/t(j+jiv,15)
              t(j+jiv,29)=xmre*(t(j+jiv,25)*t12+t(j+jiv,26)*t22+
     +                          t(j+jiv,27)*t23)/t(j+jiv,15)
              t(j+jiv,30)=xmre*(t(j+jiv,25)*t13+t(j+jiv,26)*t23+
     +                          t(j+jiv,27)*t33)/t(j+jiv,15)
              xmut=f7*t(j+jiv,24)+(1.-f7)*cmuc1*rho*zksav2**2/zksav1
              t(j+jiv,31)=xmre/t(j+jiv,15)*(
     +           t(j+jiv,25)*(u*t11 + v*t12 + w*t13)
     +          +t(j+jiv,26)*(u*t12 + v*t22 + w*t23)
     +          +t(j+jiv,27)*(u*t13 + v*t23 + w*t33))
c   to get mulam, need mulam=t14*t15/xmre-t24
     +         -xmre/t(j+jiv,15)*(t(j+jiv,14)*
     +          t(j+jiv,15)/xmre-t(j+jiv,24)+sigk1*xmut)*
     +          (t(j+jiv,25)**2 + t(j+jiv,26)**2 + t(j+jiv,27)**2)*
     +          (ti0(j,kk,2,3)-zksav(j,kk,i-1,2))
 1283         continue
            else
              do 1282 j=1,jdim1
c            Determine Sij and Wij values (at cell interfaces):
              s11 = t(j+jiv,8) *t(j+jiv,25)
              s22 = t(j+jiv,9) *t(j+jiv,26)
              s33 = t(j+jiv,10)*t(j+jiv,27)
              s12 = 0.5*(t(j+jiv,8) *t(j+jiv,26) +
     +                   t(j+jiv,9) *t(j+jiv,25))
              s13 = 0.5*(t(j+jiv,8) *t(j+jiv,27) +
     +                   t(j+jiv,10)*t(j+jiv,25))
              s23 = 0.5*(t(j+jiv,9) *t(j+jiv,27) +
     +                   t(j+jiv,10)*t(j+jiv,26))
              w12 = 0.5*(t(j+jiv,8) *t(j+jiv,26) -
     +                   t(j+jiv,9) *t(j+jiv,25))
              w13 = 0.5*(t(j+jiv,8) *t(j+jiv,27) -
     +                   t(j+jiv,10)*t(j+jiv,25))
              w23 = 0.5*(t(j+jiv,9) *t(j+jiv,27) -
     +                   t(j+jiv,10)*t(j+jiv,26))
              xis = s11*s11 + s22*s22 + s33*s33 +
     +              2.*s12*s12 + 2.*s13*s13 + 2.*s23*s23
c             Girimaji JFM 2000 fix to c4
              if (ieasm_type .eq. 4) then
                wis = 2.*w12*w12 + 2.*w13*w13 + 2.*w23*w23
                eta1_girimaji=xis/(xis+wis)
                if (real(eta1_girimaji) .lt. 0.5) then
                  c4new=2.0-((2.0-c4)*
     +                  (eta1_girimaji/(1.-eta1_girimaji))**0.75)
                else
                  c4new=c4
                end if
                al4 = 0.5*c4new -1.
              end if
c            "Constants" are function of c3, c4, and gg:
c            Find tauij values:
              rho=0.5*(q(j,kk,i,1)+q(j,kk,i-1,1))
              zksav1=0.5*(zksav(j,kk,i,1)+zksav(j,kk,i-1,1))
              zksav2=0.5*(zksav(j,kk,i,2)+zksav(j,kk,i-1,2))
              cmuu=0.5*(cmuv(j,kk,i)+cmuv(j,kk,i-1))
              bnum=zksav2*(1.-f7)+f7
c             Durbin TCFD 1991 near-wall limiter
              if (idurbinlim .ne. 0 .and. (ieasm_type .eq. 3 .or.
     +            ieasm_type .eq. 4)) then
                tau=bnum/zksav1
                fmu=t(j+jiv,14)*t(j+jiv,15)/xmre-t(j+jiv,24)
                zkolmog=6.*sqrt(fmu*bnum/(rho*zksav2*zksav1))
                tau=ccmax(tau,zkolmog)
                zksav1=bnum/tau
              end if
              zksav1=ccmax(zksav1,tur10(1))
              eta1=xis*(bnum/zksav1/re)**2
c   *********
              alpa1 = -al4/(al10-eta1*al1*cmuu)
              alpa2 = -2.*al3/(al10-eta1*cmuu*al1)
              t11 = 2.*rho*zksav2*re/3.
     +             -4.*alpa1*t(j+jiv,24)*bnum/
     +                 zksav1*(-s12*w12 - s13*w13)/re
     +           +2.*alpa2*t(j+jiv,24)*bnum/zksav1*
     +                 (s11*s11 + s12*s12 + s13*s13 - 0.33333*xis)/re
              xlnpt=2.*t(j+jiv,24)*(s11-(s11+s22+s33)/3.)
              t11=ccmax(t11,xlnpt)
              t22 = 2.*rho*zksav2*re/3.
     +             -4.*alpa1*t(j+jiv,24)*bnum/
     +                 zksav1*( s12*w12 - s23*w23)/re
     +           +2.*alpa2*t(j+jiv,24)*bnum/zksav1*
     +                 (s22*s22 + s12*s12 + s23*s23 - 0.33333*xis)/re
              xlnpt=2.*t(j+jiv,24)*(s22-(s11+s22+s33)/3.)
              t22=ccmax(t22,xlnpt)
              t33 = 2.*rho*zksav2*re/3.
     +             -4.*alpa1*t(j+jiv,24)*bnum/
     +                 zksav1*( s23*w23 + s13*w13)/re
     +           +2.*alpa2*t(j+jiv,24)*bnum/zksav1*
     +                 (s33*s33 + s23*s23 + s13*s13 - 0.33333*xis)/re
              xlnpt=2.*t(j+jiv,24)*(s33-(s11+s22+s33)/3.)
              t33=ccmax(t33,xlnpt)
              t12 =-2.*alpa1*t(j+jiv,24)*bnum/zksav1*
     +                 ( s11*w12 - s22*w12 - s13*w23 - s23*w13)/re
     +           +2.*alpa2*t(j+jiv,24)*bnum/zksav1*
     +                 (s11*s12 + s12*s22 + s13*s23)/re
              t13 =-2.*alpa1*t(j+jiv,24)*bnum/zksav1*
     +                 ( s11*w13 + s12*w23 - s23*w12 - s33*w13)/re
     +           +2.*alpa2*t(j+jiv,24)*bnum/zksav1*
     +                 (s11*s13 + s12*s23 + s13*s33)/re
              t23 =-2.*alpa1*t(j+jiv,24)*bnum/zksav1*
     +                 ( s12*w13 + s22*w23 + s13*w12 - s33*w23)/re
     +           +2.*alpa2*t(j+jiv,24)*bnum/zksav1*
     +                 (s12*s13 + s22*s23 + s23*s33)/re
c
              t(j+jiv,28)=xmre*(t(j+jiv,25)*t11+t(j+jiv,26)*t12+
     +                          t(j+jiv,27)*t13)/t(j+jiv,15)
              t(j+jiv,29)=xmre*(t(j+jiv,25)*t12+t(j+jiv,26)*t22+
     +                          t(j+jiv,27)*t23)/t(j+jiv,15)
              t(j+jiv,30)=xmre*(t(j+jiv,25)*t13+t(j+jiv,26)*t23+
     +                          t(j+jiv,27)*t33)/t(j+jiv,15)
              xmut=f7*t(j+jiv,24)+(1.-f7)*cmuc1*rho*zksav2**2/zksav1
              t(j+jiv,31)=xmre/t(j+jiv,15)*(
     +           t(j+jiv,25)*(0.5*(q(j,kk,i,2)+q(j,kk,i-1,2))*t11
     +                       +0.5*(q(j,kk,i,3)+q(j,kk,i-1,3))*t12
     +                       +0.5*(q(j,kk,i,4)+q(j,kk,i-1,4))*t13)
     +          +t(j+jiv,26)*(0.5*(q(j,kk,i,2)+q(j,kk,i-1,2))*t12
     +                       +0.5*(q(j,kk,i,3)+q(j,kk,i-1,3))*t22
     +                       +0.5*(q(j,kk,i,4)+q(j,kk,i-1,4))*t23)
     +          +t(j+jiv,27)*(0.5*(q(j,kk,i,2)+q(j,kk,i-1,2))*t13
     +                       +0.5*(q(j,kk,i,3)+q(j,kk,i-1,3))*t23
     +                       +0.5*(q(j,kk,i,4)+q(j,kk,i-1,4))*t33))
c   to get mulam, need mulam=t14*t15/xmre-t24
     +         -xmre/t(j+jiv,15)*(t(j+jiv,14)*
     +          t(j+jiv,15)/xmre-t(j+jiv,24)+sigk1*xmut)*
     +          (t(j+jiv,25)**2 + t(j+jiv,26)**2 + t(j+jiv,27)**2)*
     +          (zksav(j,kk,i,2)-zksav(j,kk,i-1,2))
 1282         continue
            end if
 1285     continue
 1286   continue
      else if(i_tauijs .ne. 1 .and. i_nonlin .ne. 0 .and. 
     +        level .ge. lglobal) then
        re=reue/xmach
c   Add nonlinear terms to RHS
        do 1296 kpl=1,npl
          kk=k+kpl-1
          jiv2=(kpl-1)*jdim1
          do 1295 i=1,idim
            jiv=(i-1)*npl*jdim1 + (kpl-1)*jdim1
            if (i.eq.1) then
              do 1294 j=1,jdim1
c            Determine Sij and Wij values (at cell interface)
              s11 = t(j+jiv,8) *t(j+jiv,25)
              s22 = t(j+jiv,9) *t(j+jiv,26)
              s33 = t(j+jiv,10)*t(j+jiv,27)
              s12 = 0.5*(t(j+jiv,8) *t(j+jiv,26) +
     +                   t(j+jiv,9) *t(j+jiv,25))
              s13 = 0.5*(t(j+jiv,8) *t(j+jiv,27) +
     +                   t(j+jiv,10)*t(j+jiv,25))
              s23 = 0.5*(t(j+jiv,9) *t(j+jiv,27) +
     +                   t(j+jiv,10)*t(j+jiv,26))
              w12 = 0.5*(t(j+jiv,8) *t(j+jiv,26) -
     +                   t(j+jiv,9) *t(j+jiv,25))
              w13 = 0.5*(t(j+jiv,8) *t(j+jiv,27) -
     +                   t(j+jiv,10)*t(j+jiv,25))
              w23 = 0.5*(t(j+jiv,9) *t(j+jiv,27) -
     +                   t(j+jiv,10)*t(j+jiv,26))
              tmp = (s11+s22+s33)/3.
              s11 = s11 - tmp
              s22 = s22 - tmp
              s33 = s33 - tmp
c             Note: denom term is taken at nearest cell center
              denom= ux(j,kk,i,1)**2 +ux(j,kk,i,2)**2 +ux(j,kk,i,3)**2
     +              +ux(j,kk,i,4)**2 +ux(j,kk,i,5)**2 +ux(j,kk,i,6)**2
     +              +ux(j,kk,i,7)**2 +ux(j,kk,i,8)**2 +ux(j,kk,i,9)**2
              denom=sqrt(denom)
              denom=ccmax(denom,snonlin_lim)
c            Find tauij values:
              rho=0.5*(q(j,kk,i,1)+qi0(j,kk,1,1))*(1.-wi0(j+jiv2,5))+
     +            qi0(j,kk,1,1)*wi0(j+jiv2,5)
              u=0.5*(q(j,kk,i,2)+qi0(j,kk,2,1))*(1.-wi0(j+jiv2,5))+
     +            qi0(j,kk,2,1)*wi0(j+jiv2,5)
              v=0.5*(q(j,kk,i,3)+qi0(j,kk,3,1))*(1.-wi0(j+jiv2,5))+
     +            qi0(j,kk,3,1)*wi0(j+jiv2,5)
              w=0.5*(q(j,kk,i,4)+qi0(j,kk,4,1))*(1.-wi0(j+jiv2,5))+
     +            qi0(j,kk,4,1)*wi0(j+jiv2,5)
              zk_sa=0.
              if (ivisc(1).eq.6 .or. ivisc(1).eq.7 .or.
     +            ivisc(1).eq.10.or. ivisc(1).eq.15) then
                zk_sa=0.5*(zksav(j,kk,i,2)+ti0(j,kk,2,1))
              end if
              t11 = 2.*rho*zk_sa*re/3.
     +             -8.*c_nonlin*t(j+jiv,24)/denom*(-s12*w12 - s13*w13)
c             xlnpt=2.*t(j+jiv,24)*s11
c             t11=ccmax(t11,xlnpt)
              t11 = t11*(1.-bci(j,kk,1))
              t22 = 2.*rho*zk_sa*re/3.
     +             -8.*c_nonlin*t(j+jiv,24)/denom*( s12*w12 - s23*w23)
c             xlnpt=2.*t(j+jiv,24)*s22
c             t22=ccmax(t22,xlnpt)
              t22 = t22*(1.-bci(j,kk,1))
              t33 = 2.*rho*zk_sa*re/3.
     +             -8.*c_nonlin*t(j+jiv,24)/denom*( s23*w23 + s13*w13)
c             xlnpt=2.*t(j+jiv,24)*s33
c             t33=ccmax(t33,xlnpt)
              t33 = t33*(1.-bci(j,kk,1))
              t12 =-4.*c_nonlin*t(j+jiv,24)/denom*
     +                 ( s11*w12 - s22*w12 - s13*w23 - s23*w13)
              t12 = t12*(1.-bci(j,kk,1))
              t13 =-4.*c_nonlin*t(j+jiv,24)/denom*
     +                 ( s11*w13 + s12*w23 - s23*w12 - s33*w13)
              t13 = t13*(1.-bci(j,kk,1))
              t23 =-4.*c_nonlin*t(j+jiv,24)/denom*
     +                 ( s12*w13 + s22*w23 + s13*w12 - s33*w23)
              t23 = t23*(1.-bci(j,kk,1))
c
              t(j+jiv,28)=xmre*(t(j+jiv,25)*t11+t(j+jiv,26)*t12+
     +                          t(j+jiv,27)*t13)/t(j+jiv,15)
              t(j+jiv,29)=xmre*(t(j+jiv,25)*t12+t(j+jiv,26)*t22+
     +                          t(j+jiv,27)*t23)/t(j+jiv,15)
              t(j+jiv,30)=xmre*(t(j+jiv,25)*t13+t(j+jiv,26)*t23+
     +                          t(j+jiv,27)*t33)/t(j+jiv,15)
              t(j+jiv,31)=xmre/t(j+jiv,15)*(
     +           t(j+jiv,25)*(u*t11 + v*t12 + w*t13)
     +          +t(j+jiv,26)*(u*t12 + v*t22 + w*t23)
     +          +t(j+jiv,27)*(u*t13 + v*t23 + w*t33))
 1294         continue
            else if(i.eq.idim) then
              do 1293 j=1,jdim1
c            Determine Sij and Wij values (at cell interface)
              s11 = t(j+jiv,8) *t(j+jiv,25)
              s22 = t(j+jiv,9) *t(j+jiv,26)
              s33 = t(j+jiv,10)*t(j+jiv,27)
              s12 = 0.5*(t(j+jiv,8) *t(j+jiv,26) +
     +                   t(j+jiv,9) *t(j+jiv,25))
              s13 = 0.5*(t(j+jiv,8) *t(j+jiv,27) +
     +                   t(j+jiv,10)*t(j+jiv,25))
              s23 = 0.5*(t(j+jiv,9) *t(j+jiv,27) +
     +                   t(j+jiv,10)*t(j+jiv,26))
              w12 = 0.5*(t(j+jiv,8) *t(j+jiv,26) -
     +                   t(j+jiv,9) *t(j+jiv,25))
              w13 = 0.5*(t(j+jiv,8) *t(j+jiv,27) -
     +                   t(j+jiv,10)*t(j+jiv,25))
              w23 = 0.5*(t(j+jiv,9) *t(j+jiv,27) -
     +                   t(j+jiv,10)*t(j+jiv,26))
              tmp = (s11+s22+s33)/3.
              s11 = s11 - tmp
              s22 = s22 - tmp
              s33 = s33 - tmp
c             Note: denom term is taken at nearest cell center
              denom= ux(j,kk,i-1,1)**2 +ux(j,kk,i-1,2)**2
     +              +ux(j,kk,i-1,3)**2
     +              +ux(j,kk,i-1,4)**2 +ux(j,kk,i-1,5)**2
     +              +ux(j,kk,i-1,6)**2
     +              +ux(j,kk,i-1,7)**2 +ux(j,kk,i-1,8)**2
     +              +ux(j,kk,i-1,9)**2
              denom=sqrt(denom)
              denom=ccmax(denom,snonlin_lim)
c            Find tauij values:
              rho=0.5*(q(j,kk,i-1,1)+qi0(j,kk,1,3))*(1.-wi0(j+jiv2,15))+
     +            qi0(j,kk,1,3)*wi0(j+jiv2,15)
              u=0.5*(q(j,kk,i-1,2)+qi0(j,kk,2,3))*(1.-wi0(j+jiv2,15))+
     +            qi0(j,kk,2,3)*wi0(j+jiv2,15)
              v=0.5*(q(j,kk,i-1,3)+qi0(j,kk,3,3))*(1.-wi0(j+jiv2,15))+
     +            qi0(j,kk,3,3)*wi0(j+jiv2,15)
              w=0.5*(q(j,kk,i-1,4)+qi0(j,kk,4,3))*(1.-wi0(j+jiv2,15))+
     +            qi0(j,kk,4,3)*wi0(j+jiv2,15)
              zk_sa=0.
              if (ivisc(1).eq.6 .or. ivisc(1).eq.7 .or.
     +            ivisc(1).eq.10.or. ivisc(1).eq.15) then
                zk_sa=0.5*(zksav(j,kk,i-1,2)+ti0(j,kk,2,3))
              end if
              t11 = 2.*rho*zk_sa*re/3.
     +             -8.*c_nonlin*t(j+jiv,24)/denom*(-s12*w12 - s13*w13)
c             xlnpt=2.*t(j+jiv,24)*s11
c             t11=ccmax(t11,xlnpt)
              t11 = t11*(1.-bci(j,kk,2))
              t22 = 2.*rho*zk_sa*re/3.
     +             -8.*c_nonlin*t(j+jiv,24)/denom*( s12*w12 - s23*w23)
c             xlnpt=2.*t(j+jiv,24)*s22
c             t22=ccmax(t22,xlnpt)
              t22 = t22*(1.-bci(j,kk,2))
              t33 = 2.*rho*zk_sa*re/3.
     +             -8.*c_nonlin*t(j+jiv,24)/denom*( s23*w23 + s13*w13)
c             xlnpt=2.*t(j+jiv,24)*s33
c             t33=ccmax(t33,xlnpt)
              t33 = t33*(1.-bci(j,kk,2))
              t12 =-4.*c_nonlin*t(j+jiv,24)/denom*
     +                 ( s11*w12 - s22*w12 - s13*w23 - s23*w13)
              t12 = t12*(1.-bci(j,kk,2))
              t13 =-4.*c_nonlin*t(j+jiv,24)/denom*
     +                 ( s11*w13 + s12*w23 - s23*w12 - s33*w13)
              t13 = t13*(1.-bci(j,kk,2))
              t23 =-4.*c_nonlin*t(j+jiv,24)/denom*
     +                 ( s12*w13 + s22*w23 + s13*w12 - s33*w23)
              t23 = t23*(1.-bci(j,kk,2))
c
              t(j+jiv,28)=xmre*(t(j+jiv,25)*t11+t(j+jiv,26)*t12+
     +                          t(j+jiv,27)*t13)/t(j+jiv,15)
              t(j+jiv,29)=xmre*(t(j+jiv,25)*t12+t(j+jiv,26)*t22+
     +                          t(j+jiv,27)*t23)/t(j+jiv,15)
              t(j+jiv,30)=xmre*(t(j+jiv,25)*t13+t(j+jiv,26)*t23+
     +                          t(j+jiv,27)*t33)/t(j+jiv,15)
              t(j+jiv,31)=xmre/t(j+jiv,15)*(
     +           t(j+jiv,25)*(u*t11 + v*t12 + w*t13)
     +          +t(j+jiv,26)*(u*t12 + v*t22 + w*t23)
     +          +t(j+jiv,27)*(u*t13 + v*t23 + w*t33))
 1293         continue
            else
              do 1292 j=1,jdim1
c            Determine Sij and Wij values (at cell interfaces):
              s11 = t(j+jiv,8) *t(j+jiv,25)
              s22 = t(j+jiv,9) *t(j+jiv,26)
              s33 = t(j+jiv,10)*t(j+jiv,27)
              s12 = 0.5*(t(j+jiv,8) *t(j+jiv,26) +
     +                   t(j+jiv,9) *t(j+jiv,25))
              s13 = 0.5*(t(j+jiv,8) *t(j+jiv,27) +
     +                   t(j+jiv,10)*t(j+jiv,25))
              s23 = 0.5*(t(j+jiv,9) *t(j+jiv,27) +
     +                   t(j+jiv,10)*t(j+jiv,26))
              w12 = 0.5*(t(j+jiv,8) *t(j+jiv,26) -
     +                   t(j+jiv,9) *t(j+jiv,25))
              w13 = 0.5*(t(j+jiv,8) *t(j+jiv,27) -
     +                   t(j+jiv,10)*t(j+jiv,25))
              w23 = 0.5*(t(j+jiv,9) *t(j+jiv,27) -
     +                   t(j+jiv,10)*t(j+jiv,26))
              tmp = (s11+s22+s33)/3.
              s11 = s11 - tmp
              s22 = s22 - tmp
              s33 = s33 - tmp
c             Note: denom term is averaged to cell face
              ux1 = 0.5*(ux(j,kk,i-1,1) + ux(j,kk,i,1))
              ux2 = 0.5*(ux(j,kk,i-1,2) + ux(j,kk,i,2))
              ux3 = 0.5*(ux(j,kk,i-1,3) + ux(j,kk,i,3))
              ux4 = 0.5*(ux(j,kk,i-1,4) + ux(j,kk,i,4))
              ux5 = 0.5*(ux(j,kk,i-1,5) + ux(j,kk,i,5))
              ux6 = 0.5*(ux(j,kk,i-1,6) + ux(j,kk,i,6))
              ux7 = 0.5*(ux(j,kk,i-1,7) + ux(j,kk,i,7))
              ux8 = 0.5*(ux(j,kk,i-1,8) + ux(j,kk,i,8))
              ux9 = 0.5*(ux(j,kk,i-1,9) + ux(j,kk,i,9))
              denom= ux1**2 +ux2**2 +ux3**2
     +              +ux4**2 +ux5**2 +ux6**2
     +              +ux7**2 +ux8**2 +ux9**2
              denom=sqrt(denom)
              denom=ccmax(denom,snonlin_lim)
c            Find tauij values:
              rho=0.5*(q(j,kk,i,1)+q(j,kk,i-1,1))
              zk_sa=0.
              if (ivisc(1).eq.6 .or. ivisc(1).eq.7 .or.
     +            ivisc(1).eq.10.or. ivisc(1).eq.15) then
                zk_sa=0.5*(zksav(j,kk,i,2)+zksav(j,kk,i-1,2))
              end if
              t11 = 2.*rho*zk_sa*re/3.
     +             -8.*c_nonlin*t(j+jiv,24)/denom*(-s12*w12 - s13*w13)
c             xlnpt=2.*t(j+jiv,24)*s11
c             t11=ccmax(t11,xlnpt)
              t22 = 2.*rho*zk_sa*re/3.
     +             -8.*c_nonlin*t(j+jiv,24)/denom*( s12*w12 - s23*w23)
c             xlnpt=2.*t(j+jiv,24)*s22
c             t22=ccmax(t22,xlnpt)
              t33 = 2.*rho*zk_sa*re/3.
     +             -8.*c_nonlin*t(j+jiv,24)/denom*( s23*w23 + s13*w13)
c             xlnpt=2.*t(j+jiv,24)*s33
c             t33=ccmax(t33,xlnpt)
              t12 =-4.*c_nonlin*t(j+jiv,24)/denom*
     +                 ( s11*w12 - s22*w12 - s13*w23 - s23*w13)
              t13 =-4.*c_nonlin*t(j+jiv,24)/denom*
     +                 ( s11*w13 + s12*w23 - s23*w12 - s33*w13)
              t23 =-4.*c_nonlin*t(j+jiv,24)/denom*
     +                 ( s12*w13 + s22*w23 + s13*w12 - s33*w23)
c
              t(j+jiv,28)=xmre*(t(j+jiv,25)*t11+t(j+jiv,26)*t12+
     +                          t(j+jiv,27)*t13)/t(j+jiv,15)
              t(j+jiv,29)=xmre*(t(j+jiv,25)*t12+t(j+jiv,26)*t22+
     +                          t(j+jiv,27)*t23)/t(j+jiv,15)
              t(j+jiv,30)=xmre*(t(j+jiv,25)*t13+t(j+jiv,26)*t23+
     +                          t(j+jiv,27)*t33)/t(j+jiv,15)
              t(j+jiv,31)=xmre/t(j+jiv,15)*(
     +           t(j+jiv,25)*(0.5*(q(j,kk,i,2)+q(j,kk,i-1,2))*t11
     +                       +0.5*(q(j,kk,i,3)+q(j,kk,i-1,3))*t12
     +                       +0.5*(q(j,kk,i,4)+q(j,kk,i-1,4))*t13)
     +          +t(j+jiv,26)*(0.5*(q(j,kk,i,2)+q(j,kk,i-1,2))*t12
     +                       +0.5*(q(j,kk,i,3)+q(j,kk,i-1,3))*t22
     +                       +0.5*(q(j,kk,i,4)+q(j,kk,i-1,4))*t23)
     +          +t(j+jiv,27)*(0.5*(q(j,kk,i,2)+q(j,kk,i-1,2))*t13
     +                       +0.5*(q(j,kk,i,3)+q(j,kk,i-1,3))*t23
     +                       +0.5*(q(j,kk,i,4)+q(j,kk,i-1,4))*t33))
 1292         continue
            end if
 1295     continue
 1296   continue
      end if
      if(i_tauijs .ne. 1 .and.
     +   (ivisc(1).eq.11.or.ivisc(1).eq.12.or.ivisc(1).eq.13.or.
     +    ivisc(1).eq.14.or.i_nonlin.ne.0)
     +    .and.level.ge.lglobal) then
        do 1271 l=2,5
        do 1270 izz=1,n
c   Take positive derivative, because want negative of all above for correct
c   sign (this term is going into -d/dxj(Hv) term)
          t(izz,l)=t(izz+js-1,l+26)-t(izz,l+26)
 1270   continue
 1271   continue
        do 1268 kpl=1,npl
          kk=k+kpl-1
          do 1268 i=1,idim1
            jiv=(i-1)*jdim1*npl + (kpl-1)*jdim1+1
            do 1268 l=2,5
              do 1267 izz=1,jdim1
                res(izz,kk,i,l)=res(izz,kk,i,l)+t(izz+jiv-1,l)
 1267         continue
 1268   continue
      end if
      return
      end
