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
      subroutine gfluxv(i,npl,jdim,kdim,idim,res,q,qj0,sj,vol,t,nvtq,
     .                  wj0,vist3d,vmuj,vj0,bcj,zksav,tj0,cmuv,volj0,
     .                  nou,bou,nbuf,ibufdim,iadv,nummem,ux)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Calculate right-hand residual contributions in the 
c     J-direction due to the viscous terms.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension sj(jdim,kdim,idim-1,5),vol(jdim,kdim,idim-1)
      dimension q(jdim,kdim,idim,5),qj0(kdim,idim-1,5,4)
      dimension res(jdim,kdim,idim-1,5),vist3d(jdim,kdim,idim),
     .          vj0(kdim,idim-1,1,4),volj0(kdim,idim-1,4)
      dimension t(nvtq,32),vmuj(kdim-1,idim-1,2),wj0(npl*(kdim-1),22)
      dimension bcj(kdim,idim-1,2)
      dimension zksav(jdim,kdim,idim,nummem),tj0(kdim,idim-1,nummem,4)
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
     .                 i_qcr2013,i_qcr2013v
      real :: coef_eddy
c
c******************************************************************************
c     right-hand-side residual contribution in J-direction (viscous terms)
c     written 3-35-88   J.L. Thomas and S.L. Taylor
c     Analytical Methods Branch, NASA Langley Research Center
c******************************************************************************
c
c     coefficient used to turn off eddy-viscosity based
c     turbulent stress calculation for ivmx==70 (full
c     Reynolds stress model)  - xxiao
c
      coef_eddy = 1.0
      if(ivisc(2)>=70) coef_eddy=0
      if(i_tauijs .eq. 1) coef_eddy=0.
      ccr2=2.5
c
      kdim1 = kdim-1
      jdim1 = jdim-1
c
c n  : number of cell centers for npl planes
c l0 : number of cell interfaces for npl planes
c kv : number of cell centers (and interfaces) on a j=constant plane
c
      n     = npl*jdim1*kdim1
      l0    = npl*kdim1*jdim
      kv    = npl*kdim1
      ks    =  kv+1
      nn    = n-kv
c
      xmre  = xmach/reue
      gpr   = gamma/pr
      gm1pr = gm1*pr
      prtr  = pr/prt
c     Durbin TCFD 1991 near-wall limiter (0=off)
c     (using it tends to delay separation - generally not desired!)
      idurbinlim=0
c
      if (isklton.gt.0 .and. i.eq.1 .and. iadv.ge.0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1442)
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
         if (ivisc(2).eq.11 .or. ivisc(2).eq.12 .or.
     +       ivisc(2).eq.13 .or. ivisc(2).eq.14 .or.
     +       i_nonlin.ne.0) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),1444)
         end if
         if (i_nonlin.ne.0) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''      snonlin_lim='',e19.8)')
     +       snonlin_lim
         end if
      end if
 1442 format(52h   computing viscous  fluxes, J-direction - central ,
     .12hdifferencing)
 1444 format(37h      nonlinear constitutive relation)
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
      l1  = jdim*kdim-1
      do 10 ipl=1,npl
      ii  = i+ipl-1
      do 50 j=1,jdim1
      jkv = (j-1)*npl*kdim1 + (ipl-1)*kdim1 
cdir$ ivdep
      do 1002 k=1,kdim1
      t(k+jkv,22) = q(j,k,ii,1)
      t(k+jkv,18) = q(j,k,ii,2)
      t(k+jkv,19) = q(j,k,ii,3)
      t(k+jkv,20) = q(j,k,ii,4)
      t(k+jkv,16) = q(j,k,ii,5) 
 1002 continue
   50 continue
   10 continue
c
c******************************************************************************
c
c  store dependent variables on j=0 and j=jdim boundaries in wj0
c
      do 1127 m=0,10,10
      if (m.eq.0) then
        mm=1
        mb=1
      else
        mm=3
        mb=2
      end if
      imin = i
      imax = i+npl-1
      kmin = 1
      kmax = kdim-1
c
      do 6120 ii=imin,imax
      kc         = (ii-i)*kdim1 + kmin-1
      do 6120 kk=kmin,kmax
      kc         = kc + 1
      wj0(kc,m+12) =  qj0(kk,ii,1,mm)
      wj0(kc,m+8)  =  qj0(kk,ii,2,mm)
      wj0(kc,m+9)  =  qj0(kk,ii,3,mm)
      wj0(kc,m+10) =  qj0(kk,ii,4,mm)
      wj0(kc,m+6)  =  qj0(kk,ii,5,mm)
c     m+5 flag indicates whether ghost cell or interface values stored in wj0
      wj0(kc,m+5)  =  bcj(kk,ii,mb)
 6120 continue
 6400 continue
 1127 continue
c
c******************************************************************************
c
c  t(17)/wj0(m+7)   : (u*u+v*v+w*w)/2
c  t(1)/wj0(m+1)    : 1/density
c  t(16)/wj0(m+6)   : c*c/(gm1*pr)
c  note: m=0 gives values at j=0, m=10 gives values at j=jdim
c
cdir$ ivdep
      do 7005 izz=1,n
      t(izz,17) = 0.5e0*(t(izz,18)*t(izz,18)
     .                  +t(izz,19)*t(izz,19)
     .                  +t(izz,20)*t(izz,20)) 
      t(izz,1)  = 1.e0/t(izz,22) 
      t(izz,16) = gpr*t(izz,16)*t(izz,1)/gm1
 7005 continue
      do 2121 m=0,10,10
cdir$ ivdep
      do 1005 izz=1,kv
      wj0(izz,m+7) = 0.5e0*(wj0(izz,m+8)*wj0(izz,m+8)
     .                     +wj0(izz,m+9)*wj0(izz,m+9)
     .                     +wj0(izz,m+10)*wj0(izz,m+10))
      wj0(izz,m+1)  = 1.e0/wj0(izz,m+12) 
      wj0(izz,m+6) = gpr*wj0(izz,m+6)*wj0(izz,m+1)/gm1
 1005 continue
 2121 continue
c
c******************************************************************************
c
c t(7) : laminar viscosity at cell centers (via sutherland relation)
c
      c2b  = cbar/tinf
      c2bp = c2b+1.e0
c
cdir$ ivdep
      do 1006 izz=1,n
      t5       = gm1pr*t(izz,16)
      t6       = sqrt(t5) 
      t(izz,7) = c2bp*t5*t6/(c2b+t5)
 1006 continue
c
c******************************************************************************
c
c t(14) : laminar viscosity values at cell interfaces
c
c     interior interfaces
cdir$ ivdep
      do 1007 izz=1,nn
      t(izz+ks-1,14) = (t(izz,7)+t(izz+ks-1,7))*.5e0
 1007 continue
c
c     j=0 and j=jdim interfaces
cdir$ ivdep
      do 1008 izz=1,kv
      ab=1.+wj0(izz,5)
      bb=1.-wj0(izz,5)
      wj05        = gm1pr*0.5*(ab*wj0(izz,6)+bb*t(izz,16))
      wj06        = sqrt(wj05) 
      t(izz,14)   = c2bp*wj05*wj06/(c2b+wj05) 
      ab=1.+wj0(izz,15)
      bb=1.-wj0(izz,15)
      wj05        = gm1pr*0.5*(ab*wj0(izz,16)+bb*t(izz+n-kv,16))
      wj06        = sqrt(wj05) 
      t(izz+n,14) = c2bp*wj05*wj06/(c2b+wj05) 
 1008 continue
c
c******************************************************************************
c
c t(15) : average jacobian (inverse volume) at cell interface
c
      n1  = kdim*jdim1 + 1
      l1  = kdim*jdim-1
      do 200 ipl=1,npl
      ii  = i+ipl-1
      do 200 j=1,jdim
      jkv  = (j-1)*npl*kdim1 + (ipl-1)*kdim1 
      if (j.eq.1) then 
c     inverse volume at j=0 interface
         do 7013 k=1,kdim1 
         t(k+jkv,15) = 2./(volj0(k,ii,1)+vol(1,k,ii))
 7013    continue
      else if (j.eq.jdim) then 
c     inverse volume at j=jdim interface
         do 7113 k=1,kdim1 
         t(k+jkv,15) = 2./(volj0(k,ii,3)+vol(jdim1,k,ii))
 7113    continue
      else
c     inverse volume at interior interfaces
         do 1013 k=1,kdim1 
         t(k+jkv,15) = 2./(vol(j,k,ii)+vol(j-1,k,ii))
 1013    continue
      end if
c
c******************************************************************************
c
c  t(25) - t(27) : components of grad(eta)
c  t1    : grad(eta)/J
c  t(25) : d(eta)/dx
c  t(26) : d(eta)/dy
c  t(27) : d(eta)/dz
c
      do 1014 k=1,kdim1 
      t1         =  sj(j,k,ii,4)*t(k+jkv,15) 
      t(k+jkv,25) =  sj(j,k,ii,1)*t1
      t(k+jkv,26) =  sj(j,k,ii,2)*t1
      t(k+jkv,27) =  sj(j,k,ii,3)*t1
 1014 continue
  200 continue
c
c******************************************************************************
c
c  t(6)-t(10) : gradients at cell interfaces
c  t(6)       : d( c*c/(gm1*pr) )/d(eta)
c  t(7)       : d( (u*u+v*v+w*w)/2 )/d(eta)
c  t(8)       : d(u)/d(eta)
c  t(9)       : d(v)/d(eta)
c  t(10)      : d(w)/d(eta)
c
      do 1000 k=1,5
      k1 = k+5
      k2 = k+15
c     interior interfaces
cdir$ ivdep
      do 1015 izz=1,nn
      t(izz+ks-1,k1) = t(izz+ks-1,k2)-t(izz,k2)
 1015 continue
c     interfaces at j=0/jdim
cdir$ ivdep
      do 1016 izz=1,kv
      z1=1.+wj0(izz,5)
      z2=1.+wj0(izz,15)
      t(izz,k1)   = z1*(t(izz,k2)-wj0(izz,k1))
      t(izz+n,k1) = z2*(wj0(izz,k2)-t(izz+n-kv,k2))
 1016 continue
 1000 continue
c
c******************************************************************************
c
c  t(24) : turbulent viscosity at interfaces (=0 for laminar flow)
c
      iviscc = ivisc(2)
      if (iviscc.gt.1) then
      do 1401 ipl=1,npl
      ii = i+ipl-1
      do 1401 j=1,jdim
      jkv=(j-1)*npl*kdim1 + (ipl-1)*kdim1
c     j=0 interface
      if (j.eq.1) then
        do 1017 k=1,kdim1
          t(k+jkv,24)=bcj(k,ii,1)*vj0(k,ii,1,1)+
     +     (1.-bcj(k,ii,1))*0.5*(vj0(k,ii,1,1)+vist3d(1,k,ii))
 1017   continue
c     j=jdim interface
      else if (j.eq.jdim) then
        do 1117 k=1,kdim1
          t(k+jkv,24)=bcj(k,ii,2)*vj0(k,ii,1,3)+
     +     (1.-bcj(k,ii,2))*0.5*(vj0(k,ii,1,3)+vist3d(jdim1,k,ii))
 1117   continue
      else
c     interior interfaces
        do 1217 k=1,kdim1
          t(k+jkv,24)=0.5*(vist3d(j,k,ii)+vist3d(j-1,k,ii))
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
cdir$ ivdep
      do 1022 izz=1,l0
c  ratio of turbulent to laminar viscosity
      t25       = t(izz,24)/t(izz,14) 
      t24       = 1.e0+t25*coef_eddy
c
c******************************************************************************
c
c  t(23) : ratio of laminar Prandtl number to total Prandtl number
c
      t(izz,23) = (1.e0+prtr*t25)/t24
c
c******************************************************************************
c
c  t(14): (mach/re)*viscosity/J
c
      t(izz,14) = xmre*t24*t(izz,14)/t(izz,15)
c
c******************************************************************************
c  t(5) : t(14)*(grad(zeta))**2
c
      t(izz,5)  = t(izz,14)*(t(izz,25)*t(izz,25)
     .                      +t(izz,26)*t(izz,26)
     .                      +t(izz,27)*t(izz,27))
 1022 continue
c
c******************************************************************************
c
c  t(11) : t(14)*((d(eta)/dx * d(u)/d(eta)) 
c               + (d(eta)/dy * d(v)/d(eta))
c               + (d(eta)/dz * d(w)/d(eta)))/3
c
      tr1        = 1.e0/3.e0
cdir$ ivdep
      do 1023 izz=1,l0
      t(izz,11)  = tr1*t(izz,14)*(t(izz,25)*t(izz,8)
     .                           +t(izz,26)*t(izz,9)
     .                           +t(izz,27)*t(izz,10))
 1023 continue
c
c******************************************************************************
c
c  store off vmuj....viscosity coefficient (laminar + turbulent) on j=0 and
c  j=jdim boundaries for later use in force integration
      do 1923 ipl=1,npl
      ii=i+ipl-1
cdir$ ivdep
      do 1922 k=1,kdim1
      izz = (ipl-1)*kdim1+k
      vmuj(k,ii,1) = t(izz,14)*t(izz,15)/xmre
      vmuj(k,ii,2) = t(izz+n,14)*t(izz+n,15)/xmre
 1922 continue
 1923 continue
c
      if (iadv .lt. 0) return
c
c******************************************************************************
c
c  t(12) : contravariant velocity V at cell interfaces
c
c     interior interfaces
cdir$ ivdep
      do 1024 izz=1,nn
      t(izz+ks-1,12) =0.5e0*(t(izz+ks-1,25)*(t(izz+ks-1,18)+t(izz,18))+
     .                       t(izz+ks-1,26)*(t(izz+ks-1,19)+t(izz,19))+
     .                       t(izz+ks-1,27)*(t(izz+ks-1,20)+t(izz,20)))
 1024 continue
c     interfaces at j=0 and j=jdim
cdir$ ivdep
      do 1025 izz=1,kv
      ab=1.+wj0(izz,15)
      bb=1.-wj0(izz,15)
      t(izz+n,12)=0.5*(t(izz+n,25)*(ab*wj0(izz,18)+bb*t(izz+n-kv,18))
     .                +t(izz+n,26)*(ab*wj0(izz,19)+bb*t(izz+n-kv,19))
     .                +t(izz+n,27)*(ab*wj0(izz,20)+bb*t(izz+n-kv,20)))
      ab=1.+wj0(izz,5)
      bb=1.-wj0(izz,5)
      t(izz,12) = 0.5*(t(izz,25)*(ab*wj0(izz,8) +bb*t(izz,18))
     .                +t(izz,26)*(ab*wj0(izz,9) +bb*t(izz,19))
     .                +t(izz,27)*(ab*wj0(izz,10)+bb*t(izz,20)))
 1025 continue
c
c******************************************************************************
c
c  calculate fluxes
c
      do 3000 l=2,4
c  viscous terms at interfaces (momentum eqs.)
cdir$ ivdep
      do 1026 izz=1,l0
      t(izz,16)      = t(izz,5)*t(izz,l+6)+t(izz,l+23)*t(izz,11)
 1026 continue
c  (-)viscous flux =  Gv(k-1/2)-Gv(k+1/2) (momentum eqs.)
cdir$ ivdep
      do 1027 izz=1,n
      t(izz,l)       = -t(izz+ks-1,16)+t(izz,16)
 1027 continue
 3000 continue
c
c  viscous terms at interfaces (energy eq.)
cdir$ ivdep
      do 1028 izz=1,l0
      t(izz,16)      =  t(izz,5)*(t(izz,23)*t(izz,6)+t(izz,7))
     .                + t(izz,12)*t(izz,11)
 1028 continue
c  (-)viscous flux = Gv(k-1/2)-Gv(k+1/2) (energy eq.)
cdir$ ivdep
      do 1029 izz=1,n
      t(izz,5)       = -t(izz+ks-1,16)+t(izz,16)
 1029 continue
c
c******************************************************************************
c
c  calculate residuals
c
      do 400 ipl=1,npl
      ii            = i+ipl-1
      do 400 j=1,jdim1
      jkv           = (j-1)*kdim1*npl + (ipl-1)*kdim1
      do 400 l=2,5
cdir$ ivdep
      do 1030 k=1,kdim1
      res(j,k,ii,l) = res(j,k,ii,l)+t(k+jkv,l)
 1030 continue
  400 continue
c  Add tau_ij's for general method:
      if (i_tauijs .eq. 1 .and. ivisc(2) .lt. 70) then
        re=reue/xmach
c   Add nonlinear terms to RHS
        do 1196 ipl=1,npl
          ii=i+ipl-1
          jkv2=(ipl-1)*kdim1
          do 1195 j=1,jdim
            jkv=(j-1)*npl*kdim1 + (ipl-1)*kdim1
            jp1=min(j+1,jdim-1)
            jm2=max(j-2,1)
            do 1194 k=1,kdim1
              if (j.eq.1) then
                rho=0.5*(q(j,k,ii,1)+qj0(k,ii,1,1))*(1.-wj0(k+jkv2,5))+
     +              qj0(k,ii,1,1)*wj0(k+jkv2,5)
                u=0.5*(q(j,k,ii,2)+qj0(k,ii,2,1))*(1.-wj0(k+jkv2,5))+
     +              qj0(k,ii,2,1)*wj0(k+jkv2,5)
                v=0.5*(q(j,k,ii,3)+qj0(k,ii,3,1))*(1.-wj0(k+jkv2,5))+
     +              qj0(k,ii,3,1)*wj0(k+jkv2,5)
                w=0.5*(q(j,k,ii,4)+qj0(k,ii,4,1))*(1.-wj0(k+jkv2,5))+
     +              qj0(k,ii,4,1)*wj0(k+jkv2,5)
             v3d=0.5*(vist3d(j,k,ii)+vj0(k,ii,1,1))*(1.-wj0(k+jkv2,5))+
     +              vj0(k,ii,1,1)*wj0(k+jkv2,5)
c               use 1st order interpolation at end
                ux1=1.5*ux(j,k,ii,1)-0.5*ux(jp1,k,ii,1)
                ux2=1.5*ux(j,k,ii,2)-0.5*ux(jp1,k,ii,2)
                ux3=1.5*ux(j,k,ii,3)-0.5*ux(jp1,k,ii,3)
                ux4=1.5*ux(j,k,ii,4)-0.5*ux(jp1,k,ii,4)
                ux5=1.5*ux(j,k,ii,5)-0.5*ux(jp1,k,ii,5)
                ux6=1.5*ux(j,k,ii,6)-0.5*ux(jp1,k,ii,6)
                ux7=1.5*ux(j,k,ii,7)-0.5*ux(jp1,k,ii,7)
                ux8=1.5*ux(j,k,ii,8)-0.5*ux(jp1,k,ii,8)
                ux9=1.5*ux(j,k,ii,9)-0.5*ux(jp1,k,ii,9)
              else if (j .eq. jdim) then
                rho=0.5*(q(j-1,k,ii,1)+qj0(k,ii,1,3))*
     +              (1.-wj0(k+jkv2,15))+
     +              qj0(k,ii,1,3)*wj0(k+jkv2,15)
                u=0.5*(q(j-1,k,ii,2)+qj0(k,ii,2,3))*(1.-wj0(k+jkv2,15))+
     +              qj0(k,ii,2,3)*wj0(k+jkv2,15)
                v=0.5*(q(j-1,k,ii,3)+qj0(k,ii,3,3))*(1.-wj0(k+jkv2,15))+
     +              qj0(k,ii,3,3)*wj0(k+jkv2,15)
                w=0.5*(q(j-1,k,ii,4)+qj0(k,ii,4,3))*(1.-wj0(k+jkv2,15))+
     +              qj0(k,ii,4,3)*wj0(k+jkv2,15)
           v3d=0.5*(vist3d(j-1,k,ii)+vj0(k,ii,1,3))*(1.-wj0(k+jkv2,15))+
     +              vj0(k,ii,1,3)*wj0(k+jkv2,15)
c               use 1st order interpolation at end
                ux1=1.5*ux(j-1,k,ii,1)-0.5*ux(jm2,k,ii,1)
                ux2=1.5*ux(j-1,k,ii,2)-0.5*ux(jm2,k,ii,2)
                ux3=1.5*ux(j-1,k,ii,3)-0.5*ux(jm2,k,ii,3)
                ux4=1.5*ux(j-1,k,ii,4)-0.5*ux(jm2,k,ii,4)
                ux5=1.5*ux(j-1,k,ii,5)-0.5*ux(jm2,k,ii,5)
                ux6=1.5*ux(j-1,k,ii,6)-0.5*ux(jm2,k,ii,6)
                ux7=1.5*ux(j-1,k,ii,7)-0.5*ux(jm2,k,ii,7)
                ux8=1.5*ux(j-1,k,ii,8)-0.5*ux(jm2,k,ii,8)
                ux9=1.5*ux(j-1,k,ii,9)-0.5*ux(jm2,k,ii,9)
              else
                rho=0.5*(q(j,k,ii,1)+q(j-1,k,ii,1))
                u=0.5*(q(j,k,ii,2)+q(j-1,k,ii,2))
                v=0.5*(q(j,k,ii,3)+q(j-1,k,ii,3))
                w=0.5*(q(j,k,ii,4)+q(j-1,k,ii,4))
                v3d=0.5*(vist3d(j,k,ii)+vist3d(j-1,k,ii))
                ux1=0.5*(ux(j,k,ii,1)+ux(j-1,k,ii,1))
                ux2=0.5*(ux(j,k,ii,2)+ux(j-1,k,ii,2))
                ux3=0.5*(ux(j,k,ii,3)+ux(j-1,k,ii,3))
                ux4=0.5*(ux(j,k,ii,4)+ux(j-1,k,ii,4))
                ux5=0.5*(ux(j,k,ii,5)+ux(j-1,k,ii,5))
                ux6=0.5*(ux(j,k,ii,6)+ux(j-1,k,ii,6))
                ux7=0.5*(ux(j,k,ii,7)+ux(j-1,k,ii,7))
                ux8=0.5*(ux(j,k,ii,8)+ux(j-1,k,ii,8))
                ux9=0.5*(ux(j,k,ii,9)+ux(j-1,k,ii,9))
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
              t(k+jkv,28)=xmre*(t(k+jkv,25)*t11+t(k+jkv,26)*t12+
     +                          t(k+jkv,27)*t13)/t(k+jkv,15)
              t(k+jkv,29)=xmre*(t(k+jkv,25)*t12+t(k+jkv,26)*t22+
     +                          t(k+jkv,27)*t23)/t(k+jkv,15)
              t(k+jkv,30)=xmre*(t(k+jkv,25)*t13+t(k+jkv,26)*t23+
     +                          t(k+jkv,27)*t33)/t(k+jkv,15)
              t(k+jkv,31)=xmre/t(k+jkv,15)*(
     +           t(k+jkv,25)*(u*t11 + v*t12 + w*t13)
     +          +t(k+jkv,26)*(u*t12 + v*t22 + w*t23)
     +          +t(k+jkv,27)*(u*t13 + v*t23 + w*t33))
 1194       continue
 1195     continue
 1196   continue
        do 1171 l=2,5
        do 1170 izz=1,n
c   Take positive derivative, because want negative of all above for correct
c   sign (this term is going into -d/dxj(Gv) term)
          t(izz,l)=t(izz+ks-1,l+26)-t(izz,l+26)
 1170   continue
 1171   continue
        do 1168 ipl=1,npl
          ii=i+ipl-1
          do 1168 j=1,jdim1
            jkv=(j-1)*kdim1*npl + (ipl-1)*kdim1+1
            do 1168 l=2,5
              do 1167 k=1,kdim1
                res(j,k,ii,l)=res(j,k,ii,l)+t(k+jkv-1,l)
 1167         continue
 1168   continue
      end if
c  Add tau_ij's for RSMs:
      if (ivisc(2) .ge. 70) then
        re=reue/xmach
c   Add nonlinear terms to RHS
        do 8296 ipl=1,npl
          ii=i+ipl-1
          jkv2=(ipl-1)*kdim1
          do 8295 j=1,jdim
            jkv=(j-1)*npl*kdim1 + (ipl-1)*kdim1
            if (j.eq.1) then
              do 8294 k=1,kdim1
              rho=0.5*(q(j,k,ii,1)+qj0(k,ii,1,1))*(1.-wj0(k+jkv2,5))+
     +            qj0(k,ii,1,1)*wj0(k+jkv2,5)
              u=0.5*(q(j,k,ii,2)+qj0(k,ii,2,1))*(1.-wj0(k+jkv2,5))+
     +            qj0(k,ii,2,1)*wj0(k+jkv2,5)
              v=0.5*(q(j,k,ii,3)+qj0(k,ii,3,1))*(1.-wj0(k+jkv2,5))+
     +            qj0(k,ii,3,1)*wj0(k+jkv2,5)
              w=0.5*(q(j,k,ii,4)+qj0(k,ii,4,1))*(1.-wj0(k+jkv2,5))+
     +            qj0(k,ii,4,1)*wj0(k+jkv2,5)
              t11=-rho*0.5*(zksav(j,k,ii,1)+tj0(k,ii,1,1))*re
              t22=-rho*0.5*(zksav(j,k,ii,2)+tj0(k,ii,2,1))*re
              t33=-rho*0.5*(zksav(j,k,ii,3)+tj0(k,ii,3,1))*re
              t12=-rho*0.5*(zksav(j,k,ii,4)+tj0(k,ii,4,1))*re
              t23=-rho*0.5*(zksav(j,k,ii,5)+tj0(k,ii,5,1))*re
              t13=-rho*0.5*(zksav(j,k,ii,6)+tj0(k,ii,6,1))*re
              t(k+jkv,28)=xmre*(t(k+jkv,25)*t11+t(k+jkv,26)*t12+
     +                          t(k+jkv,27)*t13)/t(k+jkv,15)
              t(k+jkv,29)=xmre*(t(k+jkv,25)*t12+t(k+jkv,26)*t22+
     +                          t(k+jkv,27)*t23)/t(k+jkv,15)
              t(k+jkv,30)=xmre*(t(k+jkv,25)*t13+t(k+jkv,26)*t23+
     +                          t(k+jkv,27)*t33)/t(k+jkv,15)
              t(k+jkv,31)=xmre/t(k+jkv,15)*(
     +           t(k+jkv,25)*(u*t11 + v*t12 + w*t13)
     +          +t(k+jkv,26)*(u*t12 + v*t22 + w*t23)
     +          +t(k+jkv,27)*(u*t13 + v*t23 + w*t33))
 8294         continue
            else if(j.eq.jdim) then
              do 8293 k=1,kdim1
              rho=0.5*(q(j-1,k,ii,1)+qj0(k,ii,1,3))*(1.-wj0(k+jkv2,15))+
     +            qj0(k,ii,1,3)*wj0(k+jkv2,15)
              u=0.5*(q(j-1,k,ii,2)+qj0(k,ii,2,3))*(1.-wj0(k+jkv2,15))+
     +            qj0(k,ii,2,3)*wj0(k+jkv2,15)
              v=0.5*(q(j-1,k,ii,3)+qj0(k,ii,3,3))*(1.-wj0(k+jkv2,15))+
     +            qj0(k,ii,3,3)*wj0(k+jkv2,15)
              w=0.5*(q(j-1,k,ii,4)+qj0(k,ii,4,3))*(1.-wj0(k+jkv2,15))+
     +            qj0(k,ii,4,3)*wj0(k+jkv2,15)
              t11=-rho*0.5*(zksav(j-1,k,ii,1)+tj0(k,ii,1,3))*re
              t22=-rho*0.5*(zksav(j-1,k,ii,2)+tj0(k,ii,2,3))*re
              t33=-rho*0.5*(zksav(j-1,k,ii,3)+tj0(k,ii,3,3))*re
              t12=-rho*0.5*(zksav(j-1,k,ii,4)+tj0(k,ii,4,3))*re
              t23=-rho*0.5*(zksav(j-1,k,ii,5)+tj0(k,ii,5,3))*re
              t13=-rho*0.5*(zksav(j-1,k,ii,6)+tj0(k,ii,6,3))*re
              t(k+jkv,28)=xmre*(t(k+jkv,25)*t11+t(k+jkv,26)*t12+
     +                          t(k+jkv,27)*t13)/t(k+jkv,15)
              t(k+jkv,29)=xmre*(t(k+jkv,25)*t12+t(k+jkv,26)*t22+
     +                          t(k+jkv,27)*t23)/t(k+jkv,15)
              t(k+jkv,30)=xmre*(t(k+jkv,25)*t13+t(k+jkv,26)*t23+
     +                          t(k+jkv,27)*t33)/t(k+jkv,15)
              t(k+jkv,31)=xmre/t(k+jkv,15)*(
     +           t(k+jkv,25)*(u*t11 + v*t12 + w*t13)
     +          +t(k+jkv,26)*(u*t12 + v*t22 + w*t23)
     +          +t(k+jkv,27)*(u*t13 + v*t23 + w*t33))
 8293         continue
            else
              do 8292 k=1,kdim1
              rho=0.5*(q(j,k,ii,1)+q(j-1,k,ii,1))
              t11=-rho*0.5*(zksav(j,k,ii,1)+zksav(j-1,k,ii,1))*re
              t22=-rho*0.5*(zksav(j,k,ii,2)+zksav(j-1,k,ii,2))*re
              t33=-rho*0.5*(zksav(j,k,ii,3)+zksav(j-1,k,ii,3))*re
              t12=-rho*0.5*(zksav(j,k,ii,4)+zksav(j-1,k,ii,4))*re
              t23=-rho*0.5*(zksav(j,k,ii,5)+zksav(j-1,k,ii,5))*re
              t13=-rho*0.5*(zksav(j,k,ii,6)+zksav(j-1,k,ii,6))*re
              t(k+jkv,28)=xmre*(t(k+jkv,25)*t11+t(k+jkv,26)*t12+
     +                          t(k+jkv,27)*t13)/t(k+jkv,15)
              t(k+jkv,29)=xmre*(t(k+jkv,25)*t12+t(k+jkv,26)*t22+
     +                          t(k+jkv,27)*t23)/t(k+jkv,15)
              t(k+jkv,30)=xmre*(t(k+jkv,25)*t13+t(k+jkv,26)*t23+
     +                          t(k+jkv,27)*t33)/t(k+jkv,15)
              t(k+jkv,31)=xmre/t(k+jkv,15)*(
     +           t(k+jkv,25)*(0.5*(q(j,k,ii,2)+q(j-1,k,ii,2))*t11
     +                       +0.5*(q(j,k,ii,3)+q(j-1,k,ii,3))*t12
     +                       +0.5*(q(j,k,ii,4)+q(j-1,k,ii,4))*t13)
     +          +t(k+jkv,26)*(0.5*(q(j,k,ii,2)+q(j-1,k,ii,2))*t12
     +                       +0.5*(q(j,k,ii,3)+q(j-1,k,ii,3))*t22
     +                       +0.5*(q(j,k,ii,4)+q(j-1,k,ii,4))*t23)
     +          +t(k+jkv,27)*(0.5*(q(j,k,ii,2)+q(j-1,k,ii,2))*t13
     +                       +0.5*(q(j,k,ii,3)+q(j-1,k,ii,3))*t23
     +                       +0.5*(q(j,k,ii,4)+q(j-1,k,ii,4))*t33))
 8292         continue
            end if
 8295     continue
 8296   continue
        do 8271 l=2,5
        do 8270 izz=1,n
c   Take positive derivative, because want negative of all above for correct
c   sign (this term is going into -d/dxj(Gv) term)
          t(izz,l)=t(izz+ks-1,l+26)-t(izz,l+26)
 8270   continue
 8271   continue
        do 8268 ipl=1,npl
          ii=i+ipl-1
          do 8268 j=1,jdim1
            jkv=(j-1)*kdim1*npl + (ipl-1)*kdim1+1
            do 8268 l=2,5
              do 8267 k=1,kdim1
                res(j,k,ii,l)=res(j,k,ii,l)+t(k+jkv-1,l)
 8267         continue
 8268   continue
      end if
c
      if(i_tauijs .ne. 1 .and. (ivisc(2) .eq. 11 .or. ivisc(2) .eq. 12)
     .    .and. level .ge. lglobal) then
c f7=factor for determining whether it's a k-omega or k-epsilon formulation
c   f7=0 for k-epsilon, 1 for k-omega
        f7=0.
        if (ivisc(2) .eq. 12) f7=1.
        re=reue/xmach
        gg=1./(c10+c5-1.)
c   Add nonlinear terms to RHS
        do 1276 ipl=1,npl
          ii=i+ipl-1
          jkv2=(ipl-1)*kdim1
          do 1275 j=1,jdim
            jkv=(j-1)*npl*kdim1 + (ipl-1)*kdim1
            if (j.eq.1) then
              do 1274 k=1,kdim1
c            Determine Sij and Wij values (at cell interface)
              s11 = t(k+jkv,8) *t(k+jkv,25)
              s22 = t(k+jkv,9) *t(k+jkv,26)
              s33 = t(k+jkv,10)*t(k+jkv,27)
              s12 = 0.5*(t(k+jkv,8) *t(k+jkv,26) +
     +                   t(k+jkv,9) *t(k+jkv,25))
              s13 = 0.5*(t(k+jkv,8) *t(k+jkv,27) +
     +                   t(k+jkv,10)*t(k+jkv,25))
              s23 = 0.5*(t(k+jkv,9) *t(k+jkv,27) +
     +                   t(k+jkv,10)*t(k+jkv,26))
              w12 = 0.5*(t(k+jkv,8) *t(k+jkv,26) -
     +                   t(k+jkv,9) *t(k+jkv,25))
              w13 = 0.5*(t(k+jkv,8) *t(k+jkv,27) -
     +                   t(k+jkv,10)*t(k+jkv,25))
              w23 = 0.5*(t(k+jkv,9) *t(k+jkv,27) -
     +                   t(k+jkv,10)*t(k+jkv,26))
              xis = s11*s11 + s22*s22 + s33*s33 +
     +              2.*s12*s12 + 2.*s13*s13 + 2.*s23*s23
              wis = 2.*w12*w12 + 2.*w13*w13 + 2.*w23*w23
c            "Constants" are function of c3, c4, and gg:
              alpa1=(2.-c4)/2.*gg
              alpa2=(2.-c3)*gg
c            Find tauij values:
              rho=0.5*(q(j,k,ii,1)+qj0(k,ii,1,1))*(1.-wj0(k+jkv2,5))+
     +            qj0(k,ii,1,1)*wj0(k+jkv2,5)
              u=0.5*(q(j,k,ii,2)+qj0(k,ii,2,1))*(1.-wj0(k+jkv2,5))+
     +            qj0(k,ii,2,1)*wj0(k+jkv2,5)
              v=0.5*(q(j,k,ii,3)+qj0(k,ii,3,1))*(1.-wj0(k+jkv2,5))+
     +            qj0(k,ii,3,1)*wj0(k+jkv2,5)
              w=0.5*(q(j,k,ii,4)+qj0(k,ii,4,1))*(1.-wj0(k+jkv2,5))+
     +            qj0(k,ii,4,1)*wj0(k+jkv2,5)
              zksav1=0.5*(zksav(j,k,ii,1)+tj0(k,ii,1,1))
              zksav1=ccmax(zksav1,tur10(1))
              zksav2=0.5*(zksav(j,k,ii,2)+tj0(k,ii,2,1))
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
     +             -4.*alpa1*t(k+jkv,24)*factre*bnum/
     +                 zksav1*(-s12*w12 - s13*w13)/re
     +           +2.*alpa2*t(k+jkv,24)*factre*bnum/zksav1*
     +                 (s11*s11 + s12*s12 + s13*s13 - 0.33333*xis)/re
              xlnpt=2.*t(k+jkv,24)*(s11-(s11+s22+s33)/3.)
              t11=ccmax(t11,xlnpt)
              t11 = t11*(1.-bcj(k,ii,1))
              t22 = 2.*rho*zksav2*re/3.
     +             -4.*alpa1*t(k+jkv,24)*factre*bnum/
     +                 zksav1*( s12*w12 - s23*w23)/re
     +           +2.*alpa2*t(k+jkv,24)*factre*bnum/zksav1*
     +                 (s22*s22 + s12*s12 + s23*s23 - 0.33333*xis)/re
              xlnpt=2.*t(k+jkv,24)*(s22-(s11+s22+s33)/3.)
              t22=ccmax(t22,xlnpt)
              t22 = t22*(1.-bcj(k,ii,1))
              t33 = 2.*rho*zksav2*re/3.
     +             -4.*alpa1*t(k+jkv,24)*factre*bnum/
     +                 zksav1*( s23*w23 + s13*w13)/re
     +           +2.*alpa2*t(k+jkv,24)*factre*bnum/zksav1*
     +                 (s33*s33 + s23*s23 + s13*s13 - 0.33333*xis)/re
              xlnpt=2.*t(k+jkv,24)*(s33-(s11+s22+s33)/3.)
              t33=ccmax(t33,xlnpt)
              t33 = t33*(1.-bcj(k,ii,1))
              t12 =-2.*alpa1*t(k+jkv,24)*factre*bnum/zksav1*
     +                 ( s11*w12 - s22*w12 - s13*w23 - s23*w13)/re
     +           +2.*alpa2*t(k+jkv,24)*factre*bnum/zksav1*
     +                 (s11*s12 + s12*s22 + s13*s23)/re
              t12 = t12*(1.-bcj(k,ii,1))
              t13 =-2.*alpa1*t(k+jkv,24)*factre*bnum/zksav1*
     +                 ( s11*w13 + s12*w23 - s23*w12 - s33*w13)/re
     +           +2.*alpa2*t(k+jkv,24)*factre*bnum/zksav1*
     +                 (s11*s13 + s12*s23 + s13*s33)/re
              t13 = t13*(1.-bcj(k,ii,1))
              t23 =-2.*alpa1*t(k+jkv,24)*factre*bnum/zksav1*
     +                 ( s12*w13 + s22*w23 + s13*w12 - s33*w23)/re
     +           +2.*alpa2*t(k+jkv,24)*factre*bnum/zksav1*
     +                 (s12*s13 + s22*s23 + s23*s33)/re
              t23 = t23*(1.-bcj(k,ii,1))
c
              t(k+jkv,28)=xmre*(t(k+jkv,25)*t11+t(k+jkv,26)*t12+
     +                          t(k+jkv,27)*t13)/t(k+jkv,15)
              t(k+jkv,29)=xmre*(t(k+jkv,25)*t12+t(k+jkv,26)*t22+
     +                          t(k+jkv,27)*t23)/t(k+jkv,15)
              t(k+jkv,30)=xmre*(t(k+jkv,25)*t13+t(k+jkv,26)*t23+
     +                          t(k+jkv,27)*t33)/t(k+jkv,15)
              t(k+jkv,31)=xmre/t(k+jkv,15)*(
     +           t(k+jkv,25)*(u*t11 + v*t12 + w*t13)
     +          +t(k+jkv,26)*(u*t12 + v*t22 + w*t23)
     +          +t(k+jkv,27)*(u*t13 + v*t23 + w*t33))
c   to get mulam, need mulam=t14*t15/xmre-t24
c   (currently uses eddy viscosity mu_t (t(24)), even though diffusion term 
c   in k-eqn uses cmuc1*rho*k/omega or cmuc1*rho*k^2/epsilon)
     +         -xmre/t(k+jkv,15)*(t(k+jkv,14)*
     +          t(k+jkv,15)/xmre-t(k+jkv,24)+sigk1*t(k+jkv,24))*
     +          (t(k+jkv,25)**2 + t(k+jkv,26)**2 + t(k+jkv,27)**2)*
     +          (zksav(j,k,ii,2)-tj0(k,ii,2,1))
 1274         continue
            else if(j.eq.jdim) then
              do 1273 k=1,kdim1
c            Determine Sij and Wij values (at cell interface)
              s11 = t(k+jkv,8) *t(k+jkv,25)
              s22 = t(k+jkv,9) *t(k+jkv,26)
              s33 = t(k+jkv,10)*t(k+jkv,27)
              s12 = 0.5*(t(k+jkv,8) *t(k+jkv,26) +
     +                   t(k+jkv,9) *t(k+jkv,25))
              s13 = 0.5*(t(k+jkv,8) *t(k+jkv,27) +
     +                   t(k+jkv,10)*t(k+jkv,25))
              s23 = 0.5*(t(k+jkv,9) *t(k+jkv,27) +
     +                   t(k+jkv,10)*t(k+jkv,26))
              w12 = 0.5*(t(k+jkv,8) *t(k+jkv,26) -
     +                   t(k+jkv,9) *t(k+jkv,25))
              w13 = 0.5*(t(k+jkv,8) *t(k+jkv,27) -
     +                   t(k+jkv,10)*t(k+jkv,25))
              w23 = 0.5*(t(k+jkv,9) *t(k+jkv,27) -
     +                   t(k+jkv,10)*t(k+jkv,26))
              xis = s11*s11 + s22*s22 + s33*s33 +
     +              2.*s12*s12 + 2.*s13*s13 + 2.*s23*s23
              wis = 2.*w12*w12 + 2.*w13*w13 + 2.*w23*w23
c            "Constants" are function of c3, c4, and gg:
              alpa1=(2.-c4)/2.*gg
              alpa2=(2.-c3)*gg
c            Find tauij values:
              rho=0.5*(q(j-1,k,ii,1)+qj0(k,ii,1,3))*(1.-wj0(k+jkv2,15))+
     +            qj0(k,ii,1,3)*wj0(k+jkv2,15)
              u=0.5*(q(j-1,k,ii,2)+qj0(k,ii,2,3))*(1.-wj0(k+jkv2,15))+
     +            qj0(k,ii,2,3)*wj0(k+jkv2,15)
              v=0.5*(q(j-1,k,ii,3)+qj0(k,ii,3,3))*(1.-wj0(k+jkv2,15))+
     +            qj0(k,ii,3,3)*wj0(k+jkv2,15)
              w=0.5*(q(j-1,k,ii,4)+qj0(k,ii,4,3))*(1.-wj0(k+jkv2,15))+
     +            qj0(k,ii,4,3)*wj0(k+jkv2,15)
              zksav1=0.5*(zksav(j-1,k,ii,1)+tj0(k,ii,1,3))
              zksav1=ccmax(zksav1,tur10(1))
              zksav2=0.5*(zksav(j-1,k,ii,2)+tj0(k,ii,2,3))
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
     +             -4.*alpa1*t(k+jkv,24)*factre*bnum/
     +                 zksav1*(-s12*w12 - s13*w13)/re
     +           +2.*alpa2*t(k+jkv,24)*factre*bnum/zksav1*
     +                 (s11*s11 + s12*s12 + s13*s13 - 0.33333*xis)/re
              xlnpt=2.*t(k+jkv,24)*(s11-(s11+s22+s33)/3.)
              t11=ccmax(t11,xlnpt)
              t11 = t11*(1.-bcj(k,ii,2))
              t22 = 2.*rho*zksav2*re/3.
     +             -4.*alpa1*t(k+jkv,24)*factre*bnum/
     +                 zksav1*( s12*w12 - s23*w23)/re
     +           +2.*alpa2*t(k+jkv,24)*factre*bnum/zksav1*
     +                 (s22*s22 + s12*s12 + s23*s23 - 0.33333*xis)/re
              xlnpt=2.*t(k+jkv,24)*(s22-(s11+s22+s33)/3.)
              t22=ccmax(t22,xlnpt)
              t22 = t22*(1.-bcj(k,ii,2))
              t33 = 2.*rho*zksav2*re/3.
     +             -4.*alpa1*t(k+jkv,24)*factre*bnum/
     +                 zksav1*( s23*w23 + s13*w13)/re
     +           +2.*alpa2*t(k+jkv,24)*factre*bnum/zksav1*
     +                 (s33*s33 + s23*s23 + s13*s13 - 0.33333*xis)/re
              xlnpt=2.*t(k+jkv,24)*(s33-(s11+s22+s33)/3.)
              t33=ccmax(t33,xlnpt)
              t33 = t33*(1.-bcj(k,ii,2))
              t12 =-2.*alpa1*t(k+jkv,24)*factre*bnum/zksav1*
     +                 ( s11*w12 - s22*w12 - s13*w23 - s23*w13)/re
     +           +2.*alpa2*t(k+jkv,24)*factre*bnum/zksav1*
     +                 (s11*s12 + s12*s22 + s13*s23)/re
              t12 = t12*(1.-bcj(k,ii,2))
              t13 =-2.*alpa1*t(k+jkv,24)*factre*bnum/zksav1*
     +                 ( s11*w13 + s12*w23 - s23*w12 - s33*w13)/re
     +           +2.*alpa2*t(k+jkv,24)*factre*bnum/zksav1*
     +                 (s11*s13 + s12*s23 + s13*s33)/re
              t13 = t13*(1.-bcj(k,ii,2))
              t23 =-2.*alpa1*t(k+jkv,24)*factre*bnum/zksav1*
     +                 ( s12*w13 + s22*w23 + s13*w12 - s33*w23)/re
     +           +2.*alpa2*t(k+jkv,24)*factre*bnum/zksav1*
     +                 (s12*s13 + s22*s23 + s23*s33)/re
              t23 = t23*(1.-bcj(k,ii,2))
c
              t(k+jkv,28)=xmre*(t(k+jkv,25)*t11+t(k+jkv,26)*t12+
     +                          t(k+jkv,27)*t13)/t(k+jkv,15)
              t(k+jkv,29)=xmre*(t(k+jkv,25)*t12+t(k+jkv,26)*t22+
     +                          t(k+jkv,27)*t23)/t(k+jkv,15)
              t(k+jkv,30)=xmre*(t(k+jkv,25)*t13+t(k+jkv,26)*t23+
     +                          t(k+jkv,27)*t33)/t(k+jkv,15)
              t(k+jkv,31)=xmre/t(k+jkv,15)*(
     +           t(k+jkv,25)*(u*t11 + v*t12 + w*t13)
     +          +t(k+jkv,26)*(u*t12 + v*t22 + w*t23)
     +          +t(k+jkv,27)*(u*t13 + v*t23 + w*t33))
c   to get mulam, need mulam=t14*t15/xmre-t24
c   (currently uses eddy viscosity mu_t (t(24)), even though diffusion term 
c   in k-eqn uses cmuc1*rho*k/omega or cmuc1*rho*k^2/epsilon)
     +         -xmre/t(k+jkv,15)*(t(k+jkv,14)*
     +          t(k+jkv,15)/xmre-t(k+jkv,24)+sigk1*t(k+jkv,24))*
     +          (t(k+jkv,25)**2 + t(k+jkv,26)**2 + t(k+jkv,27)**2)*
     +          (tj0(k,ii,2,3)-zksav(j-1,k,ii,2))
 1273         continue
            else
              do 1272 k=1,kdim1
c            Determine Sij and Wij values (at cell interfaces):
              s11 = t(k+jkv,8) *t(k+jkv,25)
              s22 = t(k+jkv,9) *t(k+jkv,26)
              s33 = t(k+jkv,10)*t(k+jkv,27)
              s12 = 0.5*(t(k+jkv,8) *t(k+jkv,26) +
     +                   t(k+jkv,9) *t(k+jkv,25))
              s13 = 0.5*(t(k+jkv,8) *t(k+jkv,27) +
     +                   t(k+jkv,10)*t(k+jkv,25))
              s23 = 0.5*(t(k+jkv,9) *t(k+jkv,27) +
     +                   t(k+jkv,10)*t(k+jkv,26))
              w12 = 0.5*(t(k+jkv,8) *t(k+jkv,26) -
     +                   t(k+jkv,9) *t(k+jkv,25))
              w13 = 0.5*(t(k+jkv,8) *t(k+jkv,27) -
     +                   t(k+jkv,10)*t(k+jkv,25))
              w23 = 0.5*(t(k+jkv,9) *t(k+jkv,27) -
     +                   t(k+jkv,10)*t(k+jkv,26))
              xis = s11*s11 + s22*s22 + s33*s33 +
     +              2.*s12*s12 + 2.*s13*s13 + 2.*s23*s23
              wis = 2.*w12*w12 + 2.*w13*w13 + 2.*w23*w23
c            "Constants" are function of c3, c4, and gg:
              alpa1=(2.-c4)/2.*gg
              alpa2=(2.-c3)*gg
c            Find tauij values:
              rho=0.5*(q(j,k,ii,1)+q(j-1,k,ii,1))
              zksav1=0.5*(zksav(j,k,ii,1)+zksav(j-1,k,ii,1))
              zksav1=ccmax(zksav1,tur10(1))
              zksav2=0.5*(zksav(j,k,ii,2)+zksav(j-1,k,ii,2))
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
     +             -4.*alpa1*t(k+jkv,24)*factre*bnum/
     +                 zksav1*(-s12*w12 - s13*w13)/re
     +           +2.*alpa2*t(k+jkv,24)*factre*bnum/zksav1*
     +                 (s11*s11 + s12*s12 + s13*s13 - 0.33333*xis)/re
              xlnpt=2.*t(k+jkv,24)*(s11-(s11+s22+s33)/3.)
              t11=ccmax(t11,xlnpt)
              t22 = 2.*rho*zksav2*re/3.
     +             -4.*alpa1*t(k+jkv,24)*factre*bnum/
     +                 zksav1*( s12*w12 - s23*w23)/re
     +           +2.*alpa2*t(k+jkv,24)*factre*bnum/zksav1*
     +                 (s22*s22 + s12*s12 + s23*s23 - 0.33333*xis)/re
              xlnpt=2.*t(k+jkv,24)*(s22-(s11+s22+s33)/3.)
              t22=ccmax(t22,xlnpt)
              t33 = 2.*rho*zksav2*re/3.
     +             -4.*alpa1*t(k+jkv,24)*factre*bnum/
     +                 zksav1*( s23*w23 + s13*w13)/re
     +           +2.*alpa2*t(k+jkv,24)*factre*bnum/zksav1*
     +                 (s33*s33 + s23*s23 + s13*s13 - 0.33333*xis)/re
              xlnpt=2.*t(k+jkv,24)*(s33-(s11+s22+s33)/3.)
              t33=ccmax(t33,xlnpt)
              t12 =-2.*alpa1*t(k+jkv,24)*factre*bnum/zksav1*
     +                 ( s11*w12 - s22*w12 - s13*w23 - s23*w13)/re
     +           +2.*alpa2*t(k+jkv,24)*factre*bnum/zksav1*
     +                 (s11*s12 + s12*s22 + s13*s23)/re
              t13 =-2.*alpa1*t(k+jkv,24)*factre*bnum/zksav1*
     +                 ( s11*w13 + s12*w23 - s23*w12 - s33*w13)/re
     +           +2.*alpa2*t(k+jkv,24)*factre*bnum/zksav1*
     +                 (s11*s13 + s12*s23 + s13*s33)/re
              t23 =-2.*alpa1*t(k+jkv,24)*factre*bnum/zksav1*
     +                 ( s12*w13 + s22*w23 + s13*w12 - s33*w23)/re
     +           +2.*alpa2*t(k+jkv,24)*factre*bnum/zksav1*
     +                 (s12*s13 + s22*s23 + s23*s33)/re
c
              t(k+jkv,28)=xmre*(t(k+jkv,25)*t11+t(k+jkv,26)*t12+
     +                          t(k+jkv,27)*t13)/t(k+jkv,15)
              t(k+jkv,29)=xmre*(t(k+jkv,25)*t12+t(k+jkv,26)*t22+
     +                          t(k+jkv,27)*t23)/t(k+jkv,15)
              t(k+jkv,30)=xmre*(t(k+jkv,25)*t13+t(k+jkv,26)*t23+
     +                          t(k+jkv,27)*t33)/t(k+jkv,15)
              t(k+jkv,31)=xmre/t(k+jkv,15)*(
     +           t(k+jkv,25)*(0.5*(q(j,k,ii,2)+q(j-1,k,ii,2))*t11
     +                       +0.5*(q(j,k,ii,3)+q(j-1,k,ii,3))*t12
     +                       +0.5*(q(j,k,ii,4)+q(j-1,k,ii,4))*t13)
     +          +t(k+jkv,26)*(0.5*(q(j,k,ii,2)+q(j-1,k,ii,2))*t12
     +                       +0.5*(q(j,k,ii,3)+q(j-1,k,ii,3))*t22
     +                       +0.5*(q(j,k,ii,4)+q(j-1,k,ii,4))*t23)
     +          +t(k+jkv,27)*(0.5*(q(j,k,ii,2)+q(j-1,k,ii,2))*t13
     +                       +0.5*(q(j,k,ii,3)+q(j-1,k,ii,3))*t23
     +                       +0.5*(q(j,k,ii,4)+q(j-1,k,ii,4))*t33))
c   to get mulam, need mulam=t14*t15/xmre-t24
c   (currently uses eddy viscosity mu_t (t(24)), even though diffusion term 
c   in k-eqn uses cmuc1*rho*k/omega or cmuc1*rho*k^2/epsilon)
     +         -xmre/t(k+jkv,15)*(t(k+jkv,14)*
     +          t(k+jkv,15)/xmre-t(k+jkv,24)+sigk1*t(k+jkv,24))*
     +          (t(k+jkv,25)**2 + t(k+jkv,26)**2 + t(k+jkv,27)**2)*
     +          (zksav(j,k,ii,2)-zksav(j-1,k,ii,2))
 1272         continue
            end if
 1275     continue
 1276   continue
      else if(i_tauijs .ne. 1 .and. (ivisc(2) .eq. 13 .or.
     +  ivisc(2) .eq. 14) .and. level .ge. lglobal) then
c f7=factor for determining whether it's a k-omega or k-epsilon formulation
c   f7=0 for k-epsilon, 1 for k-omega
        f7=0.
        if (ivisc(2) .eq. 14) f7=1.
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
        do 1286 ipl=1,npl
          ii=i+ipl-1
          jkv2=(ipl-1)*kdim1
          do 1285 j=1,jdim
            jkv=(j-1)*npl*kdim1 + (ipl-1)*kdim1
            if (j.eq.1) then
              do 1284 k=1,kdim1
c            Determine Sij and Wij values (at cell interface)
              s11 = t(k+jkv,8) *t(k+jkv,25)
              s22 = t(k+jkv,9) *t(k+jkv,26)
              s33 = t(k+jkv,10)*t(k+jkv,27)
              s12 = 0.5*(t(k+jkv,8) *t(k+jkv,26) +
     +                   t(k+jkv,9) *t(k+jkv,25))
              s13 = 0.5*(t(k+jkv,8) *t(k+jkv,27) +
     +                   t(k+jkv,10)*t(k+jkv,25))
              s23 = 0.5*(t(k+jkv,9) *t(k+jkv,27) +
     +                   t(k+jkv,10)*t(k+jkv,26))
              w12 = 0.5*(t(k+jkv,8) *t(k+jkv,26) -
     +                   t(k+jkv,9) *t(k+jkv,25))
              w13 = 0.5*(t(k+jkv,8) *t(k+jkv,27) -
     +                   t(k+jkv,10)*t(k+jkv,25))
              w23 = 0.5*(t(k+jkv,9) *t(k+jkv,27) -
     +                   t(k+jkv,10)*t(k+jkv,26))
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
              rho=0.5*(q(j,k,ii,1)+qj0(k,ii,1,1))*(1.-wj0(k+jkv2,5))+
     +            qj0(k,ii,1,1)*wj0(k+jkv2,5)
              u=0.5*(q(j,k,ii,2)+qj0(k,ii,2,1))*(1.-wj0(k+jkv2,5))+
     +            qj0(k,ii,2,1)*wj0(k+jkv2,5)
              v=0.5*(q(j,k,ii,3)+qj0(k,ii,3,1))*(1.-wj0(k+jkv2,5))+
     +            qj0(k,ii,3,1)*wj0(k+jkv2,5)
              w=0.5*(q(j,k,ii,4)+qj0(k,ii,4,1))*(1.-wj0(k+jkv2,5))+
     +            qj0(k,ii,4,1)*wj0(k+jkv2,5)
              zksav1=0.5*(zksav(j,k,ii,1)+tj0(k,ii,1,1))
              zksav2=0.5*(zksav(j,k,ii,2)+tj0(k,ii,2,1))
              cmuu = cmuv(j,k,ii)
              bnum=zksav2*(1.-f7)+f7
c             Durbin TCFD 1991 near-wall limiter
              if (idurbinlim .ne. 0 .and. (ieasm_type .eq. 3 .or.
     +            ieasm_type .eq. 4)) then
                tau=bnum/zksav1
                fmu=t(k+jkv,14)*t(k+jkv,15)/xmre-t(k+jkv,24)
                zkolmog=6.*sqrt(fmu*bnum/(rho*zksav2*zksav1))
                tau=ccmax(tau,zkolmog)
                zksav1=bnum/tau
              end if
              zksav1=ccmax(zksav1,tur10(1))
              eta1=xis*(bnum/zksav1/re)**2
              alpa1 = -al4/(al10-eta1*al1*cmuu)
              alpa2 = -2.*al3/(al10-eta1*cmuu*al1)
              t11 = 2.*rho*zksav2*re/3.
     +             -4.*alpa1*t(k+jkv,24)*bnum/
     +                 zksav1*(-s12*w12 - s13*w13)/re
     +           +2.*alpa2*t(k+jkv,24)*bnum/zksav1*
     +                 (s11*s11 + s12*s12 + s13*s13 - 0.33333*xis)/re
              xlnpt=2.*t(k+jkv,24)*(s11-(s11+s22+s33)/3.)
              t11=ccmax(t11,xlnpt)
              t11 = t11*(1.-bcj(k,ii,1))
              t22 = 2.*rho*zksav2*re/3.
     +             -4.*alpa1*t(k+jkv,24)*bnum/
     +                 zksav1*( s12*w12 - s23*w23)/re
     +           +2.*alpa2*t(k+jkv,24)*bnum/zksav1*
     +                 (s22*s22 + s12*s12 + s23*s23 - 0.33333*xis)/re
              xlnpt=2.*t(k+jkv,24)*(s22-(s11+s22+s33)/3.)
              t22=ccmax(t22,xlnpt)
              t22 = t22*(1.-bcj(k,ii,1))
              t33 = 2.*rho*zksav2*re/3.
     +             -4.*alpa1*t(k+jkv,24)*bnum/
     +                 zksav1*( s23*w23 + s13*w13)/re
     +           +2.*alpa2*t(k+jkv,24)*bnum/zksav1*
     +                 (s33*s33 + s23*s23 + s13*s13 - 0.33333*xis)/re
              xlnpt=2.*t(k+jkv,24)*(s33-(s11+s22+s33)/3.)
              t33=ccmax(t33,xlnpt)
              t33 = t33*(1.-bcj(k,ii,1))
              t12 =-2.*alpa1*t(k+jkv,24)*bnum/zksav1*
     +                 ( s11*w12 - s22*w12 - s13*w23 - s23*w13)/re
     +           +2.*alpa2*t(k+jkv,24)*bnum/zksav1*
     +                 (s11*s12 + s12*s22 + s13*s23)/re
              t12 = t12*(1.-bcj(k,ii,1))
              t13 =-2.*alpa1*t(k+jkv,24)*bnum/zksav1*
     +                 ( s11*w13 + s12*w23 - s23*w12 - s33*w13)/re
     +           +2.*alpa2*t(k+jkv,24)*bnum/zksav1*
     +                 (s11*s13 + s12*s23 + s13*s33)/re
              t13 = t13*(1.-bcj(k,ii,1))
              t23 =-2.*alpa1*t(k+jkv,24)*bnum/zksav1*
     +                 ( s12*w13 + s22*w23 + s13*w12 - s33*w23)/re
     +           +2.*alpa2*t(k+jkv,24)*bnum/zksav1*
     +                 (s12*s13 + s22*s23 + s23*s33)/re
              t23 = t23*(1.-bcj(k,ii,1))
c
              t(k+jkv,28)=xmre*(t(k+jkv,25)*t11+t(k+jkv,26)*t12+
     +                          t(k+jkv,27)*t13)/t(k+jkv,15)
              t(k+jkv,29)=xmre*(t(k+jkv,25)*t12+t(k+jkv,26)*t22+
     +                          t(k+jkv,27)*t23)/t(k+jkv,15)
              t(k+jkv,30)=xmre*(t(k+jkv,25)*t13+t(k+jkv,26)*t23+
     +                          t(k+jkv,27)*t33)/t(k+jkv,15)
              xmut=f7*t(k+jkv,24)+(1.-f7)*cmuc1*rho*zksav2**2/zksav1
              t(k+jkv,31)=xmre/t(k+jkv,15)*(
     +           t(k+jkv,25)*(u*t11 + v*t12 + w*t13)
     +          +t(k+jkv,26)*(u*t12 + v*t22 + w*t23)
     +          +t(k+jkv,27)*(u*t13 + v*t23 + w*t33))
c   to get mulam, need mulam=t14*t15/xmre-t24
     +         -xmre/t(k+jkv,15)*(t(k+jkv,14)*
     +          t(k+jkv,15)/xmre-t(k+jkv,24)+sigk1*xmut)*
     +          (t(k+jkv,25)**2 + t(k+jkv,26)**2 + t(k+jkv,27)**2)*
     +          (zksav(j,k,ii,2)-tj0(k,ii,2,1))
 1284         continue
            else if(j.eq.jdim) then
              do 1283 k=1,kdim1
c            Determine Sij and Wij values (at cell interface)
              s11 = t(k+jkv,8) *t(k+jkv,25)
              s22 = t(k+jkv,9) *t(k+jkv,26)
              s33 = t(k+jkv,10)*t(k+jkv,27)
              s12 = 0.5*(t(k+jkv,8) *t(k+jkv,26) +
     +                   t(k+jkv,9) *t(k+jkv,25))
              s13 = 0.5*(t(k+jkv,8) *t(k+jkv,27) +
     +                   t(k+jkv,10)*t(k+jkv,25))
              s23 = 0.5*(t(k+jkv,9) *t(k+jkv,27) +
     +                   t(k+jkv,10)*t(k+jkv,26))
              w12 = 0.5*(t(k+jkv,8) *t(k+jkv,26) -
     +                   t(k+jkv,9) *t(k+jkv,25))
              w13 = 0.5*(t(k+jkv,8) *t(k+jkv,27) -
     +                   t(k+jkv,10)*t(k+jkv,25))
              w23 = 0.5*(t(k+jkv,9) *t(k+jkv,27) -
     +                   t(k+jkv,10)*t(k+jkv,26))
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
              rho=0.5*(q(j-1,k,ii,1)+qj0(k,ii,1,3))*(1.-wj0(k+jkv2,15))+
     +            qj0(k,ii,1,3)*wj0(k+jkv2,15)
              u=0.5*(q(j-1,k,ii,2)+qj0(k,ii,2,3))*(1.-wj0(k+jkv2,15))+
     +            qj0(k,ii,2,3)*wj0(k+jkv2,15)
              v=0.5*(q(j-1,k,ii,3)+qj0(k,ii,3,3))*(1.-wj0(k+jkv2,15))+
     +            qj0(k,ii,3,3)*wj0(k+jkv2,15)
              w=0.5*(q(j-1,k,ii,4)+qj0(k,ii,4,3))*(1.-wj0(k+jkv2,15))+
     +            qj0(k,ii,4,3)*wj0(k+jkv2,15)
              zksav1=0.5*(zksav(j-1,k,ii,1)+tj0(k,ii,1,3))
              zksav2=0.5*(zksav(j-1,k,ii,2)+tj0(k,ii,2,3))
              cmuu=cmuv(j-1,k,ii)
              bnum=zksav2*(1.-f7)+f7
c             Durbin TCFD 1991 near-wall limiter
              if (idurbinlim .ne. 0 .and. (ieasm_type .eq. 3 .or.
     +            ieasm_type .eq. 4)) then
                tau=bnum/zksav1
                fmu=t(k+jkv,14)*t(k+jkv,15)/xmre-t(k+jkv,24)
                zkolmog=6.*sqrt(fmu*bnum/(rho*zksav2*zksav1))
                tau=ccmax(tau,zkolmog)
                zksav1=bnum/tau
              end if
              zksav1=ccmax(zksav1,tur10(1))
              eta1=xis*(bnum/zksav1/re)**2
              alpa1 = -al4/(al10-eta1*al1*cmuu)
              alpa2 = -2.*al3/(al10-eta1*cmuu*al1)
              t11 = 2.*rho*zksav2*re/3.
     +             -4.*alpa1*t(k+jkv,24)*bnum/
     +                 zksav1*(-s12*w12 - s13*w13)/re
     +           +2.*alpa2*t(k+jkv,24)*bnum/zksav1*
     +                 (s11*s11 + s12*s12 + s13*s13 - 0.33333*xis)/re
              xlnpt=2.*t(k+jkv,24)*(s11-(s11+s22+s33)/3.)
              t11=ccmax(t11,xlnpt)
              t11 = t11*(1.-bcj(k,ii,2))
              t22 = 2.*rho*zksav2*re/3.
     +             -4.*alpa1*t(k+jkv,24)*bnum/
     +                 zksav1*( s12*w12 - s23*w23)/re
     +           +2.*alpa2*t(k+jkv,24)*bnum/zksav1*
     +                 (s22*s22 + s12*s12 + s23*s23 - 0.33333*xis)/re
              xlnpt=2.*t(k+jkv,24)*(s22-(s11+s22+s33)/3.)
              t22=ccmax(t22,xlnpt)
              t22 = t22*(1.-bcj(k,ii,2))
              t33 = 2.*rho*zksav2*re/3.
     +             -4.*alpa1*t(k+jkv,24)*bnum/
     +                 zksav1*( s23*w23 + s13*w13)/re
     +           +2.*alpa2*t(k+jkv,24)*bnum/zksav1*
     +                 (s33*s33 + s23*s23 + s13*s13 - 0.33333*xis)/re
              xlnpt=2.*t(k+jkv,24)*(s33-(s11+s22+s33)/3.)
              t33=ccmax(t33,xlnpt)
              t33 = t33*(1.-bcj(k,ii,2))
              t12 =-2.*alpa1*t(k+jkv,24)*bnum/zksav1*
     +                 ( s11*w12 - s22*w12 - s13*w23 - s23*w13)/re
     +           +2.*alpa2*t(k+jkv,24)*bnum/zksav1*
     +                 (s11*s12 + s12*s22 + s13*s23)/re
              t12 = t12*(1.-bcj(k,ii,2))
              t13 =-2.*alpa1*t(k+jkv,24)*bnum/zksav1*
     +                 ( s11*w13 + s12*w23 - s23*w12 - s33*w13)/re
     +           +2.*alpa2*t(k+jkv,24)*bnum/zksav1*
     +                 (s11*s13 + s12*s23 + s13*s33)/re
              t13 = t13*(1.-bcj(k,ii,2))
              t23 =-2.*alpa1*t(k+jkv,24)*bnum/zksav1*
     +                 ( s12*w13 + s22*w23 + s13*w12 - s33*w23)/re
     +           +2.*alpa2*t(k+jkv,24)*bnum/zksav1*
     +                 (s12*s13 + s22*s23 + s23*s33)/re
              t23 = t23*(1.-bcj(k,ii,2))
c
              t(k+jkv,28)=xmre*(t(k+jkv,25)*t11+t(k+jkv,26)*t12+
     +                          t(k+jkv,27)*t13)/t(k+jkv,15)
              t(k+jkv,29)=xmre*(t(k+jkv,25)*t12+t(k+jkv,26)*t22+
     +                          t(k+jkv,27)*t23)/t(k+jkv,15)
              t(k+jkv,30)=xmre*(t(k+jkv,25)*t13+t(k+jkv,26)*t23+
     +                          t(k+jkv,27)*t33)/t(k+jkv,15)
              xmut=f7*t(k+jkv,24)+(1.-f7)*cmuc1*rho*zksav2**2/zksav1
              t(k+jkv,31)=xmre/t(k+jkv,15)*(
     +           t(k+jkv,25)*(u*t11 + v*t12 + w*t13)
     +          +t(k+jkv,26)*(u*t12 + v*t22 + w*t23)
     +          +t(k+jkv,27)*(u*t13 + v*t23 + w*t33))
c   to get mulam, need mulam=t14*t15/xmre-t24
     +         -xmre/t(k+jkv,15)*(t(k+jkv,14)*
     +          t(k+jkv,15)/xmre-t(k+jkv,24)+sigk1*xmut)*
     +          (t(k+jkv,25)**2 + t(k+jkv,26)**2 + t(k+jkv,27)**2)*
     +          (tj0(k,ii,2,3)-zksav(j-1,k,ii,2))
 1283         continue
            else
              do 1282 k=1,kdim1
c            Determine Sij and Wij values (at cell interfaces):
              s11 = t(k+jkv,8) *t(k+jkv,25)
              s22 = t(k+jkv,9) *t(k+jkv,26)
              s33 = t(k+jkv,10)*t(k+jkv,27)
              s12 = 0.5*(t(k+jkv,8) *t(k+jkv,26) +
     +                   t(k+jkv,9) *t(k+jkv,25))
              s13 = 0.5*(t(k+jkv,8) *t(k+jkv,27) +
     +                   t(k+jkv,10)*t(k+jkv,25))
              s23 = 0.5*(t(k+jkv,9) *t(k+jkv,27) +
     +                   t(k+jkv,10)*t(k+jkv,26))
              w12 = 0.5*(t(k+jkv,8) *t(k+jkv,26) -
     +                   t(k+jkv,9) *t(k+jkv,25))
              w13 = 0.5*(t(k+jkv,8) *t(k+jkv,27) -
     +                   t(k+jkv,10)*t(k+jkv,25))
              w23 = 0.5*(t(k+jkv,9) *t(k+jkv,27) -
     +                   t(k+jkv,10)*t(k+jkv,26))
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
              rho=0.5*(q(j,k,ii,1)+q(j-1,k,ii,1))
              zksav1=0.5*(zksav(j,k,ii,1)+zksav(j-1,k,ii,1))
              zksav2=0.5*(zksav(j,k,ii,2)+zksav(j-1,k,ii,2))
              cmuu=0.5*(cmuv(j,k,ii)+cmuv(j-1,k,ii))
              bnum=zksav2*(1.-f7)+f7
c             Durbin TCFD 1991 near-wall limiter
              if (idurbinlim .ne. 0 .and. (ieasm_type .eq. 3 .or.
     +            ieasm_type .eq. 4)) then
                tau=bnum/zksav1
                fmu=t(k+jkv,14)*t(k+jkv,15)/xmre-t(k+jkv,24)
                zkolmog=6.*sqrt(fmu*bnum/(rho*zksav2*zksav1))
                tau=ccmax(tau,zkolmog)
                zksav1=bnum/tau
              end if
              zksav1=ccmax(zksav1,tur10(1))
              eta1=xis*(bnum/zksav1/re)**2
              alpa1 = -al4/(al10-eta1*al1*cmuu)
              alpa2 = -2.*al3/(al10-eta1*cmuu*al1)
              t11 = 2.*rho*zksav2*re/3.
     +             -4.*alpa1*t(k+jkv,24)*bnum/
     +                 zksav1*(-s12*w12 - s13*w13)/re
     +           +2.*alpa2*t(k+jkv,24)*bnum/zksav1*
     +                 (s11*s11 + s12*s12 + s13*s13 - 0.33333*xis)/re
              xlnpt=2.*t(k+jkv,24)*(s11-(s11+s22+s33)/3.)
              t11=ccmax(t11,xlnpt)
              t22 = 2.*rho*zksav2*re/3.
     +             -4.*alpa1*t(k+jkv,24)*bnum/
     +                 zksav1*( s12*w12 - s23*w23)/re
     +           +2.*alpa2*t(k+jkv,24)*bnum/zksav1*
     +                 (s22*s22 + s12*s12 + s23*s23 - 0.33333*xis)/re
              xlnpt=2.*t(k+jkv,24)*(s22-(s11+s22+s33)/3.)
              t22=ccmax(t22,xlnpt)
              t33 = 2.*rho*zksav2*re/3.
     +             -4.*alpa1*t(k+jkv,24)*bnum/
     +                 zksav1*( s23*w23 + s13*w13)/re
     +           +2.*alpa2*t(k+jkv,24)*bnum/zksav1*
     +                 (s33*s33 + s23*s23 + s13*s13 - 0.33333*xis)/re
              xlnpt=2.*t(k+jkv,24)*(s33-(s11+s22+s33)/3.)
              t33=ccmax(t33,xlnpt)
              t12 =-2.*alpa1*t(k+jkv,24)*bnum/zksav1*
     +                 ( s11*w12 - s22*w12 - s13*w23 - s23*w13)/re
     +           +2.*alpa2*t(k+jkv,24)*bnum/zksav1*
     +                 (s11*s12 + s12*s22 + s13*s23)/re
              t13 =-2.*alpa1*t(k+jkv,24)*bnum/zksav1*
     +                 ( s11*w13 + s12*w23 - s23*w12 - s33*w13)/re
     +           +2.*alpa2*t(k+jkv,24)*bnum/zksav1*
     +                 (s11*s13 + s12*s23 + s13*s33)/re
              t23 =-2.*alpa1*t(k+jkv,24)*bnum/zksav1*
     +                 ( s12*w13 + s22*w23 + s13*w12 - s33*w23)/re
     +           +2.*alpa2*t(k+jkv,24)*bnum/zksav1*
     +                 (s12*s13 + s22*s23 + s23*s33)/re
c
              t(k+jkv,28)=xmre*(t(k+jkv,25)*t11+t(k+jkv,26)*t12+
     +                          t(k+jkv,27)*t13)/t(k+jkv,15)
              t(k+jkv,29)=xmre*(t(k+jkv,25)*t12+t(k+jkv,26)*t22+
     +                          t(k+jkv,27)*t23)/t(k+jkv,15)
              t(k+jkv,30)=xmre*(t(k+jkv,25)*t13+t(k+jkv,26)*t23+
     +                          t(k+jkv,27)*t33)/t(k+jkv,15)
              xmut=f7*t(k+jkv,24)+(1.-f7)*cmuc1*rho*zksav2**2/zksav1
              t(k+jkv,31)=xmre/t(k+jkv,15)*(
     +           t(k+jkv,25)*(0.5*(q(j,k,ii,2)+q(j-1,k,ii,2))*t11
     +                       +0.5*(q(j,k,ii,3)+q(j-1,k,ii,3))*t12
     +                       +0.5*(q(j,k,ii,4)+q(j-1,k,ii,4))*t13)
     +          +t(k+jkv,26)*(0.5*(q(j,k,ii,2)+q(j-1,k,ii,2))*t12
     +                       +0.5*(q(j,k,ii,3)+q(j-1,k,ii,3))*t22
     +                       +0.5*(q(j,k,ii,4)+q(j-1,k,ii,4))*t23)
     +          +t(k+jkv,27)*(0.5*(q(j,k,ii,2)+q(j-1,k,ii,2))*t13
     +                       +0.5*(q(j,k,ii,3)+q(j-1,k,ii,3))*t23
     +                       +0.5*(q(j,k,ii,4)+q(j-1,k,ii,4))*t33))
c   to get mulam, need mulam=t14*t15/xmre-t24
     +         -xmre/t(k+jkv,15)*(t(k+jkv,14)*
     +          t(k+jkv,15)/xmre-t(k+jkv,24)+sigk1*xmut)*
     +          (t(k+jkv,25)**2 + t(k+jkv,26)**2 + t(k+jkv,27)**2)*
     +          (zksav(j,k,ii,2)-zksav(j-1,k,ii,2))
 1282         continue
            end if
 1285     continue
 1286   continue
      else if(i_tauijs .ne. 1 .and. i_nonlin .ne. 0 .and. 
     +        level .ge. lglobal) then
        re=reue/xmach
c   Add nonlinear terms to RHS
        do 1296 ipl=1,npl
          ii=i+ipl-1
          jkv2=(ipl-1)*kdim1
          do 1295 j=1,jdim
            jkv=(j-1)*npl*kdim1 + (ipl-1)*kdim1
            if (j.eq.1) then
              do 1294 k=1,kdim1
c            Determine Sij and Wij values (at cell interface)
              s11 = t(k+jkv,8) *t(k+jkv,25)
              s22 = t(k+jkv,9) *t(k+jkv,26)
              s33 = t(k+jkv,10)*t(k+jkv,27)
              s12 = 0.5*(t(k+jkv,8) *t(k+jkv,26) +
     +                   t(k+jkv,9) *t(k+jkv,25))
              s13 = 0.5*(t(k+jkv,8) *t(k+jkv,27) +
     +                   t(k+jkv,10)*t(k+jkv,25))
              s23 = 0.5*(t(k+jkv,9) *t(k+jkv,27) +
     +                   t(k+jkv,10)*t(k+jkv,26))
              w12 = 0.5*(t(k+jkv,8) *t(k+jkv,26) -
     +                   t(k+jkv,9) *t(k+jkv,25))
              w13 = 0.5*(t(k+jkv,8) *t(k+jkv,27) -
     +                   t(k+jkv,10)*t(k+jkv,25))
              w23 = 0.5*(t(k+jkv,9) *t(k+jkv,27) -
     +                   t(k+jkv,10)*t(k+jkv,26))
              tmp = (s11+s22+s33)/3.
              s11 = s11 - tmp
              s22 = s22 - tmp
              s33 = s33 - tmp
c             Note: denom term is taken at nearest cell center
              denom= ux(j,k,ii,1)**2 +ux(j,k,ii,2)**2 +ux(j,k,ii,3)**2
     +              +ux(j,k,ii,4)**2 +ux(j,k,ii,5)**2 +ux(j,k,ii,6)**2
     +              +ux(j,k,ii,7)**2 +ux(j,k,ii,8)**2 +ux(j,k,ii,9)**2
              denom=sqrt(denom)
              denom=ccmax(denom,snonlin_lim)
c            Find tauij values:
              rho=0.5*(q(j,k,ii,1)+qj0(k,ii,1,1))*(1.-wj0(k+jkv2,5))+
     +            qj0(k,ii,1,1)*wj0(k+jkv2,5)
              u=0.5*(q(j,k,ii,2)+qj0(k,ii,2,1))*(1.-wj0(k+jkv2,5))+
     +            qj0(k,ii,2,1)*wj0(k+jkv2,5)
              v=0.5*(q(j,k,ii,3)+qj0(k,ii,3,1))*(1.-wj0(k+jkv2,5))+
     +            qj0(k,ii,3,1)*wj0(k+jkv2,5)
              w=0.5*(q(j,k,ii,4)+qj0(k,ii,4,1))*(1.-wj0(k+jkv2,5))+
     +            qj0(k,ii,4,1)*wj0(k+jkv2,5)
              zk_sa=0.
              if (ivisc(2).eq.6 .or. ivisc(2).eq.7 .or.
     +            ivisc(2).eq.10.or. ivisc(2).eq.15) then
                zk_sa=0.5*(zksav(j,k,ii,2)+tj0(k,ii,2,1))
              end if
              t11 = 2.*rho*zk_sa*re/3.
     +             -8.*c_nonlin*t(k+jkv,24)/denom*(-s12*w12 - s13*w13)
c             xlnpt=2.*t(k+jkv,24)*s11
c             t11=ccmax(t11,xlnpt)
              t11 = t11*(1.-bcj(k,ii,1))
              t22 = 2.*rho*zk_sa*re/3.
     +             -8.*c_nonlin*t(k+jkv,24)/denom*( s12*w12 - s23*w23)
c             xlnpt=2.*t(k+jkv,24)*s22
c             t22=ccmax(t22,xlnpt)
              t22 = t22*(1.-bcj(k,ii,1))
              t33 = 2.*rho*zk_sa*re/3.
     +             -8.*c_nonlin*t(k+jkv,24)/denom*( s23*w23 + s13*w13)
c             xlnpt=2.*t(k+jkv,24)*s33
c             t33=ccmax(t33,xlnpt)
              t33 = t33*(1.-bcj(k,ii,1))
              t12 =-4.*c_nonlin*t(k+jkv,24)/denom*
     +                 ( s11*w12 - s22*w12 - s13*w23 - s23*w13)
              t12 = t12*(1.-bcj(k,ii,1))
              t13 =-4.*c_nonlin*t(k+jkv,24)/denom*
     +                 ( s11*w13 + s12*w23 - s23*w12 - s33*w13)
              t13 = t13*(1.-bcj(k,ii,1))
              t23 =-4.*c_nonlin*t(k+jkv,24)/denom*
     +                 ( s12*w13 + s22*w23 + s13*w12 - s33*w23)
              t23 = t23*(1.-bcj(k,ii,1))
c
              t(k+jkv,28)=xmre*(t(k+jkv,25)*t11+t(k+jkv,26)*t12+
     +                          t(k+jkv,27)*t13)/t(k+jkv,15)
              t(k+jkv,29)=xmre*(t(k+jkv,25)*t12+t(k+jkv,26)*t22+
     +                          t(k+jkv,27)*t23)/t(k+jkv,15)
              t(k+jkv,30)=xmre*(t(k+jkv,25)*t13+t(k+jkv,26)*t23+
     +                          t(k+jkv,27)*t33)/t(k+jkv,15)
              t(k+jkv,31)=xmre/t(k+jkv,15)*(
     +           t(k+jkv,25)*(u*t11 + v*t12 + w*t13)
     +          +t(k+jkv,26)*(u*t12 + v*t22 + w*t23)
     +          +t(k+jkv,27)*(u*t13 + v*t23 + w*t33))
 1294         continue
            else if(j.eq.jdim) then
              do 1293 k=1,kdim1
c            Determine Sij and Wij values (at cell interface)
              s11 = t(k+jkv,8) *t(k+jkv,25)
              s22 = t(k+jkv,9) *t(k+jkv,26)
              s33 = t(k+jkv,10)*t(k+jkv,27)
              s12 = 0.5*(t(k+jkv,8) *t(k+jkv,26) +
     +                   t(k+jkv,9) *t(k+jkv,25))
              s13 = 0.5*(t(k+jkv,8) *t(k+jkv,27) +
     +                   t(k+jkv,10)*t(k+jkv,25))
              s23 = 0.5*(t(k+jkv,9) *t(k+jkv,27) +
     +                   t(k+jkv,10)*t(k+jkv,26))
              w12 = 0.5*(t(k+jkv,8) *t(k+jkv,26) -
     +                   t(k+jkv,9) *t(k+jkv,25))
              w13 = 0.5*(t(k+jkv,8) *t(k+jkv,27) -
     +                   t(k+jkv,10)*t(k+jkv,25))
              w23 = 0.5*(t(k+jkv,9) *t(k+jkv,27) -
     +                   t(k+jkv,10)*t(k+jkv,26))
              tmp = (s11+s22+s33)/3.
              s11 = s11 - tmp
              s22 = s22 - tmp
              s33 = s33 - tmp
c             Note: denom term is taken at nearest cell center
              denom= ux(j-1,k,ii,1)**2 +ux(j-1,k,ii,2)**2
     +              +ux(j-1,k,ii,3)**2
     +              +ux(j-1,k,ii,4)**2 +ux(j-1,k,ii,5)**2
     +              +ux(j-1,k,ii,6)**2
     +              +ux(j-1,k,ii,7)**2 +ux(j-1,k,ii,8)**2
     +              +ux(j-1,k,ii,9)**2
              denom=sqrt(denom)
              denom=ccmax(denom,snonlin_lim)
c            Find tauij values:
              rho=0.5*(q(j-1,k,ii,1)+qj0(k,ii,1,3))*(1.-wj0(k+jkv2,15))+
     +            qj0(k,ii,1,3)*wj0(k+jkv2,15)
              u=0.5*(q(j-1,k,ii,2)+qj0(k,ii,2,3))*(1.-wj0(k+jkv2,15))+
     +            qj0(k,ii,2,3)*wj0(k+jkv2,15)
              v=0.5*(q(j-1,k,ii,3)+qj0(k,ii,3,3))*(1.-wj0(k+jkv2,15))+
     +            qj0(k,ii,3,3)*wj0(k+jkv2,15)
              w=0.5*(q(j-1,k,ii,4)+qj0(k,ii,4,3))*(1.-wj0(k+jkv2,15))+
     +            qj0(k,ii,4,3)*wj0(k+jkv2,15)
              zk_sa=0.
              if (ivisc(2).eq.6 .or. ivisc(2).eq.7 .or.
     +            ivisc(2).eq.10.or. ivisc(2).eq.15) then
                zk_sa=0.5*(zksav(j-1,k,ii,2)+tj0(k,ii,2,3))
              end if
              t11 = 2.*rho*zk_sa*re/3.
     +             -8.*c_nonlin*t(k+jkv,24)/denom*(-s12*w12 - s13*w13)
c             xlnpt=2.*t(k+jkv,24)*s11
c             t11=ccmax(t11,xlnpt)
              t11 = t11*(1.-bcj(k,ii,2))
              t22 = 2.*rho*zk_sa*re/3.
     +             -8.*c_nonlin*t(k+jkv,24)/denom*( s12*w12 - s23*w23)
c             xlnpt=2.*t(k+jkv,24)*s22
c             t22=ccmax(t22,xlnpt)
              t22 = t22*(1.-bcj(k,ii,2))
              t33 = 2.*rho*zk_sa*re/3.
     +             -8.*c_nonlin*t(k+jkv,24)/denom*( s23*w23 + s13*w13)
c             xlnpt=2.*t(k+jkv,24)*s33
c             t33=ccmax(t33,xlnpt)
              t33 = t33*(1.-bcj(k,ii,2))
              t12 =-4.*c_nonlin*t(k+jkv,24)/denom*
     +                 ( s11*w12 - s22*w12 - s13*w23 - s23*w13)
              t12 = t12*(1.-bcj(k,ii,2))
              t13 =-4.*c_nonlin*t(k+jkv,24)/denom*
     +                 ( s11*w13 + s12*w23 - s23*w12 - s33*w13)
              t13 = t13*(1.-bcj(k,ii,2))
              t23 =-4.*c_nonlin*t(k+jkv,24)/denom*
     +                 ( s12*w13 + s22*w23 + s13*w12 - s33*w23)
              t23 = t23*(1.-bcj(k,ii,2))
c
              t(k+jkv,28)=xmre*(t(k+jkv,25)*t11+t(k+jkv,26)*t12+
     +                          t(k+jkv,27)*t13)/t(k+jkv,15)
              t(k+jkv,29)=xmre*(t(k+jkv,25)*t12+t(k+jkv,26)*t22+
     +                          t(k+jkv,27)*t23)/t(k+jkv,15)
              t(k+jkv,30)=xmre*(t(k+jkv,25)*t13+t(k+jkv,26)*t23+
     +                          t(k+jkv,27)*t33)/t(k+jkv,15)
              t(k+jkv,31)=xmre/t(k+jkv,15)*(
     +           t(k+jkv,25)*(u*t11 + v*t12 + w*t13)
     +          +t(k+jkv,26)*(u*t12 + v*t22 + w*t23)
     +          +t(k+jkv,27)*(u*t13 + v*t23 + w*t33))
 1293         continue
            else
              do 1292 k=1,kdim1
c            Determine Sij and Wij values (at cell interfaces):
              s11 = t(k+jkv,8) *t(k+jkv,25)
              s22 = t(k+jkv,9) *t(k+jkv,26)
              s33 = t(k+jkv,10)*t(k+jkv,27)
              s12 = 0.5*(t(k+jkv,8) *t(k+jkv,26) +
     +                   t(k+jkv,9) *t(k+jkv,25))
              s13 = 0.5*(t(k+jkv,8) *t(k+jkv,27) +
     +                   t(k+jkv,10)*t(k+jkv,25))
              s23 = 0.5*(t(k+jkv,9) *t(k+jkv,27) +
     +                   t(k+jkv,10)*t(k+jkv,26))
              w12 = 0.5*(t(k+jkv,8) *t(k+jkv,26) -
     +                   t(k+jkv,9) *t(k+jkv,25))
              w13 = 0.5*(t(k+jkv,8) *t(k+jkv,27) -
     +                   t(k+jkv,10)*t(k+jkv,25))
              w23 = 0.5*(t(k+jkv,9) *t(k+jkv,27) -
     +                   t(k+jkv,10)*t(k+jkv,26))
              tmp = (s11+s22+s33)/3.
              s11 = s11 - tmp
              s22 = s22 - tmp
              s33 = s33 - tmp
c             Note: denom term is averaged to cell face
              ux1 = 0.5*(ux(j-1,k,ii,1) + ux(j,k,ii,1))
              ux2 = 0.5*(ux(j-1,k,ii,2) + ux(j,k,ii,2))
              ux3 = 0.5*(ux(j-1,k,ii,3) + ux(j,k,ii,3))
              ux4 = 0.5*(ux(j-1,k,ii,4) + ux(j,k,ii,4))
              ux5 = 0.5*(ux(j-1,k,ii,5) + ux(j,k,ii,5))
              ux6 = 0.5*(ux(j-1,k,ii,6) + ux(j,k,ii,6))
              ux7 = 0.5*(ux(j-1,k,ii,7) + ux(j,k,ii,7))
              ux8 = 0.5*(ux(j-1,k,ii,8) + ux(j,k,ii,8))
              ux9 = 0.5*(ux(j-1,k,ii,9) + ux(j,k,ii,9))
              denom= ux1**2 +ux2**2 +ux3**2
     +              +ux4**2 +ux5**2 +ux6**2
     +              +ux7**2 +ux8**2 +ux9**2
              denom=sqrt(denom)
              denom=ccmax(denom,snonlin_lim)
c            Find tauij values:
              rho=0.5*(q(j,k,ii,1)+q(j-1,k,ii,1))
              zk_sa=0.
              if (ivisc(2).eq.6 .or. ivisc(2).eq.7 .or.
     +            ivisc(2).eq.10.or. ivisc(2).eq.15) then
                zk_sa=0.5*(zksav(j,k,ii,2)+zksav(j-1,k,ii,2))
              end if
              t11 = 2.*rho*zk_sa*re/3.
     +             -8.*c_nonlin*t(k+jkv,24)/denom*(-s12*w12 - s13*w13)
c             xlnpt=2.*t(k+jkv,24)*s11
c             t11=ccmax(t11,xlnpt)
              t22 = 2.*rho*zk_sa*re/3.
     +             -8.*c_nonlin*t(k+jkv,24)/denom*( s12*w12 - s23*w23)
c             xlnpt=2.*t(k+jkv,24)*s22
c             t22=ccmax(t22,xlnpt)
              t33 = 2.*rho*zk_sa*re/3.
     +             -8.*c_nonlin*t(k+jkv,24)/denom*( s23*w23 + s13*w13)
c             xlnpt=2.*t(k+jkv,24)*s33
c             t33=ccmax(t33,xlnpt)
              t12 =-4.*c_nonlin*t(k+jkv,24)/denom*
     +                 ( s11*w12 - s22*w12 - s13*w23 - s23*w13)
              t13 =-4.*c_nonlin*t(k+jkv,24)/denom*
     +                 ( s11*w13 + s12*w23 - s23*w12 - s33*w13)
              t23 =-4.*c_nonlin*t(k+jkv,24)/denom*
     +                 ( s12*w13 + s22*w23 + s13*w12 - s33*w23)
c
              t(k+jkv,28)=xmre*(t(k+jkv,25)*t11+t(k+jkv,26)*t12+
     +                          t(k+jkv,27)*t13)/t(k+jkv,15)
              t(k+jkv,29)=xmre*(t(k+jkv,25)*t12+t(k+jkv,26)*t22+
     +                          t(k+jkv,27)*t23)/t(k+jkv,15)
              t(k+jkv,30)=xmre*(t(k+jkv,25)*t13+t(k+jkv,26)*t23+
     +                          t(k+jkv,27)*t33)/t(k+jkv,15)
              t(k+jkv,31)=xmre/t(k+jkv,15)*(
     +           t(k+jkv,25)*(0.5*(q(j,k,ii,2)+q(j-1,k,ii,2))*t11
     +                       +0.5*(q(j,k,ii,3)+q(j-1,k,ii,3))*t12
     +                       +0.5*(q(j,k,ii,4)+q(j-1,k,ii,4))*t13)
     +          +t(k+jkv,26)*(0.5*(q(j,k,ii,2)+q(j-1,k,ii,2))*t12
     +                       +0.5*(q(j,k,ii,3)+q(j-1,k,ii,3))*t22
     +                       +0.5*(q(j,k,ii,4)+q(j-1,k,ii,4))*t23)
     +          +t(k+jkv,27)*(0.5*(q(j,k,ii,2)+q(j-1,k,ii,2))*t13
     +                       +0.5*(q(j,k,ii,3)+q(j-1,k,ii,3))*t23
     +                       +0.5*(q(j,k,ii,4)+q(j-1,k,ii,4))*t33))
 1292         continue
            end if
 1295     continue
 1296   continue
      end if
      if(i_tauijs .ne. 1 .and.
     +   (ivisc(2).eq.11.or.ivisc(2).eq.12.or.ivisc(2).eq.13.or.
     +    ivisc(2).eq.14.or.i_nonlin.ne.0)
     +    .and.level.ge.lglobal) then
        do 1271 l=2,5
        do 1270 izz=1,n
c   Take positive derivative, because want negative of all above for correct
c   sign (this term is going into -d/dxj(Gv) term)
          t(izz,l)=t(izz+ks-1,l+26)-t(izz,l+26)
 1270   continue
 1271   continue
        do 1268 ipl=1,npl
          ii=i+ipl-1
          do 1268 j=1,jdim1
            jkv=(j-1)*kdim1*npl + (ipl-1)*kdim1+1
            do 1268 l=2,5
              do 1267 k=1,kdim1
                res(j,k,ii,l)=res(j,k,ii,l)+t(k+jkv-1,l)
 1267         continue
 1268   continue
      end if
      return
      end 
