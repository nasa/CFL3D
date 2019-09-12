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
      subroutine hfluxv1(i,npl,jdim,kdim,idim,res,q,qj0,qk0,qi0,
     .                  sj,sk,si,vol,t,nvtq,wk0,vist3d,vk0,
     .                  bcj,bck,bci,volk0,
     .                  nou,bou,nbuf,ibufdim,iadv)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Calculate right-hand residual contributions in the
c     K-direction due to the viscous terms missing in the thin-layer
c     formulation.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension sk(jdim*kdim,idim-1,5),vol(jdim,kdim,idim-1)
      dimension si(jdim,kdim,idim,5),sj(jdim,kdim,idim-1,5)

      dimension q(jdim,kdim,idim,5)
      dimension qj0(kdim,idim-1,5,4), qk0(jdim,idim-1,5,4),
     .          qi0(jdim,kdim,5,4)
      dimension res(jdim,kdim,idim-1,5),vist3d(jdim,kdim,idim),
     .          vk0(jdim,idim-1,1,4), volk0(jdim,idim-1,4)
      dimension t(nvtq,32),wk0(npl*(jdim-1),22)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
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
      common /twod/ i2d
      common /constit/ i_nonlin,c_nonlin,snonlin_lim,i_tauijs,i_qcr2000,
     .                 i_qcr2013,i_qcr2013v
      REAL :: coef_eddy
c
c******************************************************************************
c     written mid-2005   V. Venkatakrishnan & C. Rumsey
c******************************************************************************
c
      coef_eddy = 1.0
      if(ivisc(3)>=70) coef_eddy = 0.0
      if(i_tauijs .eq. 1) coef_eddy=0.
      kdim1 = kdim-1
      jdim1 = jdim-1
      idim1 = idim-1
c
c n  : number of cell centers for npl planes
c l0 : number of cell interfaces for npl planes
c jv : number of cell centers (and interfaces) on a k=constant plane
c
      n     = npl*jdim1*kdim1
      l0    = npl*jdim1*kdim
      jv    = npl*jdim1
      js    = jv+1
      nn    = n-jv
c
      xmre  = xmach/reue
      gpr   = gamma/pr
      gm1pr = gm1*pr
      prtr  = pr/prt
      gprgm1 = gpr/gm1
c
      if (isklton.gt.0 .and. i.eq.1 .and. iadv.ge.0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1443)
      end if
 1443 format(52h   computing cross-derivative viscous fluxes, K-dir.)
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
      l1 = jdim*kdim-1
      do 10 ipl=1,npl
      ii = i+ipl-1
      do 50 k=1,kdim1
      jkv = (k-1)*npl*jdim1 + (ipl-1)*jdim1 + 1
cdir$ ivdep
      do 1002 izz=1,jdim1
      t(izz+jkv-1,22) = q(izz,k,ii,1)
      t(izz+jkv-1,18) = q(izz,k,ii,2)
      t(izz+jkv-1,19) = q(izz,k,ii,3)
      t(izz+jkv-1,20) = q(izz,k,ii,4)
      t(izz+jkv-1,16) = q(izz,k,ii,5) 
 1002 continue
   50 continue
   10 continue
c
c******************************************************************************
c
c  store dependent variables on k=0 and k=kdim boundaries in wk0
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
      imin = i
      imax = i+npl-1
      jmin = 1
      jmax = jdim-1
c
      do 6120 ii=imin,imax
      jc         = (ii-i)*jdim1 + jmin-1
      do 6120 jj=jmin,jmax
      jc         = jc + 1
      wk0(jc,m+12) = qk0(jj,ii,1,mm)
      wk0(jc,m+8)  = qk0(jj,ii,2,mm)
      wk0(jc,m+9)  = qk0(jj,ii,3,mm)
      wk0(jc,m+10) = qk0(jj,ii,4,mm) 
      wk0(jc,m+6)  = qk0(jj,ii,5,mm) 
c     m+5 flag indicates whether ghost cell or interface values stored in wk0
      wk0(jc,m+5)  =  bck(jj,ii,mb)
 6120 continue
 1127 continue
c
c******************************************************************************
c
c  t(17)/wk0(m+7)   : (u*u+v*v+w*w)/2
c  t(1)/wk0(m+1)    : 1/density
c  t(16)/wk0(m+6)   : c*c/(gm1*pr)
c  note: m=0 gives values at k=0, m=10 gives values at k=kdim
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
      do 1005 izz=1,jv
      wk0(izz,m+7) = 0.5e0*(wk0(izz,m+8)*wk0(izz,m+8)
     .                     +wk0(izz,m+9)*wk0(izz,m+9)
     .                     +wk0(izz,m+10)*wk0(izz,m+10))
      wk0(izz,m+1)  = 1.e0/wk0(izz,m+12) 
      wk0(izz,m+6) = gpr*wk0(izz,m+6)*wk0(izz,m+1)/gm1
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
      t(izz+js-1,14) = (t(izz,7)+t(izz+js-1,7))*.5e0
 1007 continue
c
c     k=0 and k=kdim interfaces
cdir$ ivdep
      do 1008 izz=1,jv
      ab=1.+wk0(izz,5)
      bb=1.-wk0(izz,5)
      wk05        = gm1pr*0.5*(ab*wk0(izz,6)+bb*t(izz,16))
      wk06        = sqrt(wk05) 
      t(izz,14)   = c2bp*wk05*wk06/(c2b+wk05)
      ab=1.+wk0(izz,15)
      bb=1.-wk0(izz,15)
      wk05        = gm1pr*0.5*(ab*wk0(izz,16)+bb*t(izz+n-jv,16))
      wk06        = sqrt(wk05) 
      t(izz+n,14) = c2bp*wk05*wk06/(c2b+wk05)
 1008 continue
c
c******************************************************************************
c
c t(15) : average jacobian (inverse volume) at cell interface 
c
      n1  = jdim*kdim1 + 1
      l1  = jdim*kdim-1
      do 200 ipl=1,npl
      ii  = i+ipl-1
      do 200 k=1,kdim
      jkv = (k-1)*npl*jdim1 + (ipl-1)*jdim1 
      jk  = (k-1)*jdim 
      if (k.eq.1) then 
c     inverse volume at k=0 interface
        do 7013 j=1,jdim1 
        t(j+jkv,15) = 2./(volk0(j,ii,1)+vol(j,1,ii))
 7013   continue
      else if (k.eq.kdim) then 
c     inverse volume at k=kdim interface
        do 7113 j=1,jdim1 
        t(j+jkv,15) = 2./(volk0(j,ii,3)+vol(j,kdim1,ii))
 7113   continue
      else
c     inverse volume at interior interfaces
        do 1013 j=1,jdim1 
        t(j+jkv,15) = 2./(vol(j,k,ii)+vol(j,k-1,ii))
 1013   continue
      end if
c
c******************************************************************************
c
c  t(25) - t(27) : components of grad(zeta)
c  t1    : grad(zeta)/J 
c  t(25) : d(zeta)/dx
c  t(26) : d(zeta)/dy
c  t(27) : d(zeta)/dz
c
      do 1014 j=1,jdim1 
      t1          = sk(j+jk,ii,4)*t(j+jkv,15) 
      t(j+jkv,25) = sk(j+jk,ii,1)*t1
      t(j+jkv,26) = sk(j+jk,ii,2)*t1
      t(j+jkv,27) = sk(j+jk,ii,3)*t1
 1014 continue
  200 continue
c
c t(28)  - t(30) : Components of grad (xie)
c grad(xie) = 1/4 [(si(j,k,i,l)*si(j,k,i,4) + si(j,k-1,i,l)*si(j,k-1,i,4) + 
c si(j,k,i+1,l)*si(j,k,i+1,4) + si(j,k-1,i+1,l)*si(j,k-1,i+1,4)] /V_int,
c j=1,jdim1 , k=1,kdim, i=1,npl, l =1,3.  At the boundaries in k, 
c since we don't have si's, we assume that them to be equal to the interior one.
c
c Note: the code for k=1 is obtained from the general case by substituting
c k for k-1, and likewise the code for k=kdim by substituting k-1 for k.
c
      do 210 ipl=1,npl
      ii  = i+ipl-1
      do 210 k=1,kdim
      jkv  = (k-1)*npl*jdim1 + (ipl-1)*jdim1 
      if (k.eq.1) then 
       do 204 j = 1,jdim1
        t(jkv+j,28) = 0.25*t(j+jkv,15)*(si(j,k,ii,1)*si(j,k,ii,4)
     .              + si(j,k,ii,1)*si(j,k,ii,4)
     .              + si(j,k,ii+1,1)*si(j,k,ii+1,4)
     .              + si(j,k,ii+1,1)*si(j,k,ii+1,4))
        t(jkv+j,29) = 0.25*t(j+jkv,15)*(si(j,k,ii,2)*si(j,k,ii,4)
     .              + si(j,k,ii,2)*si(j,k,ii,4)
     .              + si(j,k,ii+1,2)*si(j,k,ii+1,4)
     .              + si(j,k,ii+1,2)*si(j,k,ii+1,4))
        t(jkv+j,30) = 0.25*t(j+jkv,15)*(si(j,k,ii,3)*si(j,k,ii,4)
     .              + si(j,k,ii,3)*si(j,k,ii,4)
     .              + si(j,k,ii+1,3)*si(j,k,ii+1,4)
     .              + si(j,k,ii+1,3)*si(j,k,ii+1,4))
  204  continue
      else if (k.eq.kdim) then
       do 206 j = 1,jdim1
        t(jkv+j,28) = 0.25*t(j+jkv,15)*(si(j,k-1,ii,1)*si(j,k-1,ii,4)
     .              + si(j,k-1,ii,1)*si(j,k-1,ii,4)
     .              + si(j,k-1,ii+1,1)*si(j,k-1,ii+1,4)
     .              + si(j,k-1,ii+1,1)*si(j,k-1,ii+1,4))
        t(jkv+j,29) = 0.25*t(j+jkv,15)*(si(j,k-1,ii,2)*si(j,k-1,ii,4)
     .              + si(j,k-1,ii,2)*si(j,k-1,ii,4)
     .              + si(j,k-1,ii+1,2)*si(j,k-1,ii+1,4)
     .              + si(j,k-1,ii+1,2)*si(j,k-1,ii+1,4))
        t(jkv+j,30) = 0.25*t(j+jkv,15)*(si(j,k-1,ii,3)*si(j,k-1,ii,4)
     .              + si(j,k-1,ii,3)*si(j,k-1,ii,4)
     .              + si(j,k-1,ii+1,3)*si(j,k-1,ii+1,4)
     .              + si(j,k-1,ii+1,3)*si(j,k-1,ii+1,4))
  206  continue
      else
c
c-- general case
c
       do 208 j = 1,jdim1
        t(jkv+j,28) = 0.25*t(j+jkv,15)*(si(j,k,ii,1)*si(j,k,ii,4)
     .              + si(j,k-1,ii,1)*si(j,k-1,ii,4)
     .              + si(j,k,ii+1,1)*si(j,k,ii+1,4)
     .              + si(j,k-1,ii+1,1)*si(j,k-1,ii+1,4))
        t(jkv+j,29) = 0.25*t(j+jkv,15)*(si(j,k,ii,2)*si(j,k,ii,4)
     .              + si(j,k-1,ii,2)*si(j,k-1,ii,4)
     .              + si(j,k,ii+1,2)*si(j,k,ii+1,4)
     .              + si(j,k-1,ii+1,2)*si(j,k-1,ii+1,4))
        t(jkv+j,30) = 0.25*t(j+jkv,15)*(si(j,k,ii,3)*si(j,k,ii,4)
     .              + si(j,k-1,ii,3)*si(j,k-1,ii,4)
     .              + si(j,k,ii+1,3)*si(j,k,ii+1,4)
     .              + si(j,k-1,ii+1,3)*si(j,k-1,ii+1,4))
  208  continue
      end if
  210 continue
c
c t(31),t(32),t(2)  : Components of grad (eta)
c grad(eta) = 1/4 [(sj(j,k,i,l)*sj(j,k,i,4) + sj(j,k-1,i,l)*sj(j,k-1,i,4) +
c sj(j+1,k,i,l)*sj(j+1,k,i,4) + sj(j+1,k-1,i,l)*sj(j+1,k-1,i,4)] /V_int,
c k=1,kdim , j=1,jdim1, i=1,npl,l = 1,3.  At the boundaries in j, 
c since we don't have si's, we assume that them to be equal to the interior one.
c
c Note: the code for j=1 is obtained from the general case by substituting
c j for j-1, and likewise the code for j=dim by substituting j-1 for j.
c
      do 220 ipl=1,npl
      ii  = i+ipl-1
      do 220 k=1,kdim
      jkv  = (k-1)*npl*jdim1 + (ipl-1)*jdim1 
      if (k.eq.1) then 
       do 214 j = 1,jdim1
        t(jkv+j,31) = 0.25*t(j+jkv,15)*(sj(j,k,ii,1)*sj(j,k,ii,4)
     .              + sj(j,k,ii,1)*sj(j,k,ii,4)
     .              + sj(j+1,k,ii,1)*sj(j+1,k,ii,4)
     .              + sj(j+1,k,ii,1)*sj(j+1,k,ii,4))
        t(jkv+j,32) = 0.25*t(j+jkv,15)*(sj(j,k,ii,2)*sj(j,k,ii,4)
     .              + sj(j,k,ii,2)*sj(j,k,ii,4)
     .              + sj(j+1,k,ii,2)*sj(j+1,k,ii,4)
     .              + sj(j+1,k,ii,2)*sj(j+1,k,ii,4))
        t(jkv+j,2)  = 0.25*t(j+jkv,15)*(sj(j,k,ii,3)*sj(j,k,ii,4)
     .              + sj(j,k,ii,3)*sj(j,k,ii,4)
     .              + sj(j+1,k,ii,3)*sj(j+1,k,ii,4)
     .              + sj(j+1,k,ii,3)*sj(j+1,k,ii,4))
  214  continue
      else if (k.eq.kdim) then
       do 216 j = 1,jdim1
        t(jkv+j,31) = 0.25*t(j+jkv,15)*(sj(j,k-1,ii,1)*sj(j,k-1,ii,4)
     .              + sj(j,k-1,ii,1)*sj(j,k-1,ii,4)
     .              + sj(j+1,k-1,ii,1)*sj(j+1,k-1,ii,4)
     .              + sj(j+1,k-1,ii,1)*sj(j+1,k-1,ii,4))
        t(jkv+j,32) = 0.25*t(j+jkv,15)*(sj(j,k-1,ii,2)*sj(j,k-1,ii,4)
     .              + sj(j,k-1,ii,2)*sj(j,k-1,ii,4)
     .              + sj(j+1,k-1,ii,2)*sj(j+1,k-1,ii,4)
     .              + sj(j+1,k-1,ii,2)*sj(j+1,k-1,ii,4))
        t(jkv+j,2)  = 0.25*t(j+jkv,15)*(sj(j,k-1,ii,3)*sj(j,k-1,ii,4)
     .              + sj(j,k-1,ii,3)*sj(j,k-1,ii,4)
     .              + sj(j+1,k-1,ii,3)*sj(j+1,k-1,ii,4)
     .              + sj(j+1,k-1,ii,3)*sj(j+1,k-1,ii,4))
  216  continue
      else
c
c-- general case
c
       do 218 j = 1,jdim1
        t(jkv+j,31) = 0.25*t(j+jkv,15)*(sj(j,k,ii,1)*sj(j,k,ii,4)
     .              + sj(j,k-1,ii,1)*sj(j,k-1,ii,4)
     .              + sj(j+1,k,ii,1)*sj(j+1,k,ii,4)
     .              + sj(j+1,k-1,ii,1)*sj(j+1,k-1,ii,4))
        t(jkv+j,32) = 0.25*t(j+jkv,15)*(sj(j,k,ii,2)*sj(j,k,ii,4)
     .              + sj(j,k-1,ii,2)*sj(j,k-1,ii,4)
     .              + sj(j+1,k,ii,2)*sj(j+1,k,ii,4)
     .              + sj(j+1,k-1,ii,2)*sj(j+1,k-1,ii,4))
        t(jkv+j,2)  = 0.25*t(j+jkv,15)*(sj(j,k,ii,3)*sj(j,k,ii,4)
     .              + sj(j,k-1,ii,3)*sj(j,k-1,ii,4)
     .              + sj(j+1,k,ii,3)*sj(j+1,k,ii,4)
     .              + sj(j+1,k-1,ii,3)*sj(j+1,k-1,ii,4))
  218  continue
      end if
  220 continue
c
c
c******************************************************************************
c
c   gradients at cell interfaces (needed for full ns terms). Note that 
c  we don't need derivatives in zeta direction, because they have already
c  been accounted for in hfluxv.
c  Recall bcj=0 => cell data, bcj=1 => face data
c
c  t(7)       : d( c*c/(gm1*pr) )/d(zeta)
c  t(8)       : d(u)/d(eta)
c  t(9)       : d(v)/d(eta)
c  t(10)      : d(w)/d(eta)
c
c     interior interfaces
c
      do 230 ipl = 1,npl
      ii = i + ipl -1
      do 230 k = 1,kdim
      jkv  = (k-1)*npl*jdim1 + (ipl-1)*jdim1 
c
      if (k.eq.1) then
c
c replace (k-1) data by boundary data
c for missing data at corners use one-sided differencing of known boundary values
c
        do 225 j = 1,jdim1
          if (j. eq. 1) then
c
c replace (j-1) data by boundary data 
c
           t(jkv+j,7) = 0.25*gprgm1*(
     .       q(j+1,k,ii,5)/ q(j+1,k,ii,1) -
     .       (1.+bcj(k,ii,1))*qj0(k,ii,5,1)/qj0(k,ii,1,1) +
     .       bcj(k,ii,1)*q(j,k,ii,5)/q(j,k,ii,1) +
     .       2.*(1.+bck(j+1,ii,1))*qk0(j+1,ii,5,1)/qk0(j+1,ii,1,1) -
     .       2.*bck(j+1,ii,1)*q(j+1,k,ii,5)/q(j+1,k,ii,1) -
     .       2.*(1.+bck(j,ii,1))*qk0(j,ii,5,1)/qk0(j,ii,1,1) +
     .       2.*bck(j,ii,1)*q(j,k,ii,5)/q(j,k,ii,1) )

           t(jkv+j,8) = 0.25*(
     .       q(j+1,k,ii,2) -
     .       (1.+bcj(k,ii,1))*qj0(k,ii,2,1) +
     .       bcj(k,ii,1)*q(j,k,ii,2) +
     .       2.*(1.+bck(j+1,ii,1))*qk0(j+1,ii,2,1) -
     .       2.*bck(j+1,ii,1)*q(j+1,k,ii,2) -
     .       2.*(1.+bck(j,ii,1))*qk0(j,ii,2,1) +
     .       2.*bck(j,ii,1)*q(j,k,ii,2) )

           t(jkv+j,9) = 0.25*(
     .       q(j+1,k,ii,3) -
     .       (1.+bcj(k,ii,1))*qj0(k,ii,3,1) +
     .       bcj(k,ii,1)*q(j,k,ii,3) +
     .       2.*(1.+bck(j+1,ii,1))*qk0(j+1,ii,3,1) -
     .       2.*bck(j+1,ii,1)*q(j+1,k,ii,3) -
     .       2.*(1.+bck(j,ii,1))*qk0(j,ii,3,1) +
     .       2.*bck(j,ii,1)*q(j,k,ii,3) )

           t(jkv+j,10) = 0.25*(
     .       q(j+1,k,ii,4) -
     .       (1.+bcj(k,ii,1))*qj0(k,ii,4,1) +
     .       bcj(k,ii,1)*q(j,k,ii,4) +
     .       2.*(1.+bck(j+1,ii,1))*qk0(j+1,ii,4,1) -
     .       2.*bck(j+1,ii,1)*q(j+1,k,ii,4) -
     .       2.*(1.+bck(j,ii,1))*qk0(j,ii,4,1) +
     .       2.*bck(j,ii,1)*q(j,k,ii,4) )

          else if (j .eq. jdim1) then
c
c replace (j+1) data by boundary data
c
           t(jkv+j,7) = 0.25*gprgm1*(
     .       (1.+bcj(k,ii,2))*qj0(k,ii,5,3)/qj0(k,ii,1,3) -
     .       bcj(k,ii,2)*q(j,k,ii,5)/q(j,k,ii,1) -
     .       q(j-1,k,ii,5)/q(j-1,k,ii,1)  +
     .       2.*(1.+bck(j,ii,1))*qk0(j,ii,5,1)/qk0(j,ii,1,1) -
     .       2.*bck(j,ii,1)*q(j,k,ii,5)/q(j,k,ii,1) -
     .       2.*(1.+bck(j-1,ii,1))*qk0(j-1,ii,5,1)/qk0(j-1,ii,1,1) +
     .       2.*bck(j-1,ii,1)*q(j-1,k,ii,5)/q(j-1,k,ii,1) )

           t(jkv+j,8) = 0.25*(
     .       (1.+bcj(k,ii,2))*qj0(k,ii,2,3) -
     .       bcj(k,ii,2)*q(j,k,ii,2) -
     .       q(j-1,k,ii,2)  +
     .       2.*(1.+bck(j,ii,1))*qk0(j,ii,2,1) -
     .       2.*bck(j,ii,1)*q(j,k,ii,2) -
     .       2.*(1.+bck(j-1,ii,1))*qk0(j-1,ii,2,1) +
     .       2.*bck(j-1,ii,1)*q(j-1,k,ii,2) )

           t(jkv+j,9) = 0.25*(
     .       (1.+bcj(k,ii,2))*qj0(k,ii,3,3) -
     .       bcj(k,ii,2)*q(j,k,ii,3) -
     .       q(j-1,k,ii,3)  +
     .       2.*(1.+bck(j,ii,1))*qk0(j,ii,3,1) -
     .       2.*bck(j,ii,1)*q(j,k,ii,3) -
     .       2.*(1.+bck(j-1,ii,1))*qk0(j-1,ii,3,1) +
     .       2.*bck(j-1,ii,1)*q(j-1,k,ii,3) )

           t(jkv+j,10) = 0.25*(
     .       (1.+bcj(k,ii,2))*qj0(k,ii,4,3) -
     .       bcj(k,ii,2)*q(j,k,ii,4) - 
     .       q(j-1,k,ii,4)  +
     .       2.*(1.+bck(j,ii,1))*qk0(j,ii,4,1) -
     .       2.*bck(j,ii,1)*q(j,k,ii,4) -
     .       2.*(1.+bck(j-1,ii,1))*qk0(j-1,ii,4,1) +
     .       2.*bck(j-1,ii,1)*q(j-1,k,ii,4) )

          else
c
c-- general case for k = 1
c
            t(jkv+j,7) = 0.25*gprgm1*(
     .       q(j+1,k,ii,5)/q(j+1,k,ii,1) -
     .       q(j-1,k,ii,5)/q(j-1,k,ii,1) +
     .       (1.+bck(j+1,ii,1))*qk0(j+1,ii,5,1)/qk0(j+1,ii,1,1) -
     .       bck(j+1,ii,1)*q(j+1,k,ii,5)/q(j+1,k,ii,1) -
     .       (1.+bck(j-1,ii,1))*qk0(j-1,ii,5,1)/qk0(j-1,ii,1,1) +
     .       bck(j-1,ii,1)*q(j-1,k,ii,5)/q(j-1,k,ii,1) )

            t(jkv+j,8) = 0.25*(
     .       q(j+1,k,ii,2) -
     .       q(j-1,k,ii,2) +
     .       (1.+bck(j+1,ii,1))*qk0(j+1,ii,2,1) -
     .       bck(j+1,ii,1)*q(j+1,k,ii,2) -
     .       (1.+bck(j-1,ii,1))*qk0(j-1,ii,2,1) +
     .       bck(j-1,ii,1)*q(j-1,k,ii,2) )

            t(jkv+j,9) = 0.25*(
     .       q(j+1,k,ii,3) -
     .       q(j-1,k,ii,3) +
     .       (1.+bck(j+1,ii,1))*qk0(j+1,ii,3,1) -
     .       bck(j+1,ii,1)*q(j+1,k,ii,3) -
     .       (1.+bck(j-1,ii,1))*qk0(j-1,ii,3,1) +
     .       bck(j-1,ii,1)*q(j-1,k,ii,3) )

            t(jkv+j,10)= 0.25*(
     .       q(j+1,k,ii,4) -
     .       q(j-1,k,ii,4) +
     .       (1.+bck(j+1,ii,1))*qk0(j+1,ii,4,1) -
     .       bck(j+1,ii,1)*q(j+1,k,ii,4) -
     .       (1.+bck(j-1,ii,1))*qk0(j-1,ii,4,1) +
     .       bck(j-1,ii,1)*q(j-1,k,ii,4) )
          end if
  225   continue

c
      else if (k.eq.kdim) then
c
c replace k data by boundary data
c for missing data at corners use one-sided differencing of known boundary values
c
        do 227 j = 1,jdim1
          if (j. eq. 1) then
c
c replace (j-1) data by boundary data 
c
           t(jkv+j,7) = 0.25*gprgm1*(
     .       2.*(1.+bck(j+1,ii,2))*qk0(j+1,ii,5,3)/qk0(j+1,ii,1,3) -
     .       2.*bck(j+1,ii,2)*q(j+1,k-1,ii,5)/q(j+1,k-1,ii,1) -
     .       2.*(1.+bck(j,ii,2))*qk0(j,ii,5,3)/qk0(j,ii,1,3) +
     .       2.*bck(j,ii,2)*q(j,k-1,ii,5)/q(j,k-1,ii,1) +
     .       q(j+1,k-1,ii,5)/q(j+1,k-1,ii,1) -
     .       (1.+bcj(k-1,ii,1))*qj0(k-1,ii,5,1)/qj0(k-1,ii,1,1) +
     .       bcj(k-1,ii,1)*q(j,k-1,ii,5)/q(j,k-1,ii,1) )

           t(jkv+j,8) = 0.25*(
     .       2.*(1.+bck(j+1,ii,2))*qk0(j+1,ii,2,3) -
     .       2.*bck(j+1,ii,2)*q(j+1,k-1,ii,2) -
     .       2.*(1.+bck(j,ii,2))*qk0(j,ii,2,3) +
     .       2.*bck(j,ii,2)*q(j,k-1,ii,2) +
     .       q(j+1,k-1,ii,2) -
     .       (1.+bcj(k-1,ii,1))*qj0(k-1,ii,2,1) +
     .       bcj(k-1,ii,1)*q(j,k-1,ii,2) )

           t(jkv+j,9) = 0.25*(
     .       2.*(1.+bck(j+1,ii,2))*qk0(j+1,ii,3,3) -
     .       2.*bck(j+1,ii,2)*q(j+1,k-1,ii,3) -
     .       2.*(1.+bck(j,ii,2))*qk0(j,ii,3,3) +
     .       2.*bck(j,ii,2)*q(j,k-1,ii,3) +
     .       q(j+1,k-1,ii,3) -
     .       (1.+bcj(k-1,ii,1))*qj0(k-1,ii,3,1) +
     .       bcj(k-1,ii,1)*q(j,k-1,ii,3) )

           t(jkv+j,10) = 0.25*(
     .       2.*(1.+bck(j+1,ii,2))*qk0(j+1,ii,4,3) -
     .       2.*bck(j+1,ii,2)*q(j+1,k-1,ii,4) -
     .       2.*(1.+bck(j,ii,2))*qk0(j,ii,4,3) +
     .       2.*bck(j,ii,2)*q(j,k-1,ii,4) +
     .       q(j+1,k-1,ii,4) -
     .       (1.+bcj(k-1,ii,1))*qj0(k-1,ii,4,1) +
     .       bcj(k-1,ii,1)*q(j,k-1,ii,4) )

          else if (j .eq. jdim1) then
c
c replace (j+1) data by boundary data 
c
           t(jkv+j,7) = 0.25*gprgm1*(
     .       2.*(1.+bck(j,ii,2))*qk0(j,ii,5,3)/qk0(j,ii,1,3) -
     .       2.*bck(j,ii,2)*q(j,k-1,ii,5)/q(j,k-1,ii,1) -
     .       2.*(1.+bck(j-1,ii,2))*qk0(j-1,ii,5,3)/qk0(j-1,ii,1,3) +
     .       2.*bck(j-1,ii,2)*q(j-1,k-1,ii,5)/q(j-1,k-1,ii,1) +
     .       (1.+bcj(k-1,ii,2))*qj0(k-1,ii,5,3)/qj0(k-1,ii,1,3) -
     .       bcj(k-1,ii,2)*q(j,k-1,ii,5)/q(j,k-1,ii,1) -
     .       q(j-1,k-1,ii,5)/q(j-1,k-1,ii,1) )

           t(jkv+j,8) = 0.25*(
     .       2.*(1.+bck(j,ii,2))*qk0(j,ii,2,3) -
     .       2.*bck(j,ii,2)*q(j,k-1,ii,2) -
     .       2.*(1.+bck(j-1,ii,2))*qk0(j-1,ii,2,3) +
     .       2.*bck(j-1,ii,2)*q(j-1,k-1,ii,2) +
     .       (1.+bcj(k-1,ii,2))*qj0(k-1,ii,2,3) -
     .       bcj(k-1,ii,2)*q(j,k-1,ii,2) -
     .       q(j-1,k-1,ii,2) )

           t(jkv+j,9) = 0.25*(
     .       2.*(1.+bck(j,ii,2))*qk0(j,ii,3,3) -
     .       2.*bck(j,ii,2)*q(j,k-1,ii,3) -
     .       2.*(1.+bck(j-1,ii,2))*qk0(j-1,ii,3,3) +
     .       2.*bck(j-1,ii,2)*q(j-1,k-1,ii,3) +
     .       (1.+bcj(k-1,ii,2))*qj0(k-1,ii,3,3) -
     .       bcj(k-1,ii,2)*q(j,k-1,ii,3) -
     .       q(j-1,k-1,ii,3) )

           t(jkv+j,10) = 0.25*(
     .       2.*(1.+bck(j,ii,2))*qk0(j,ii,4,3) -
     .       2.*bck(j,ii,2)*q(j,k-1,ii,4) -
     .       2.*(1.+bck(j-1,ii,2))*qk0(j-1,ii,4,3) +
     .       2.*bck(j-1,ii,2)*q(j-1,k-1,ii,4) +
     .       (1.+bcj(k-1,ii,2))*qj0(k-1,ii,4,3) -
     .       bcj(k-1,ii,2)*q(j,k-1,ii,4) - 
     .       q(j-1,k-1,ii,4) )

          else
c
c-- general case for k = kdim
c
           t(jkv+j,7) = 0.25*gprgm1*(
     .       (1.+bck(j+1,ii,2))*qk0(j+1,ii,5,3)/qk0(j+1,ii,1,3) -
     .       bck(j+1,ii,2)*q(j+1,k-1,ii,5)/q(j+1,k-1,ii,1) -
     .       (1.+bck(j-1,ii,2))*qk0(j-1,ii,5,3)/qk0(j-1,ii,1,3) +
     .       bck(j-1,ii,2)*q(j-1,k-1,ii,5)/q(j-1,k-1,ii,1) +
     .       q(j+1,k-1,ii,5)/q(j+1,k-1,ii,1) -
     .       q(j-1,k-1,ii,5)/q(j-1,k-1,ii,1) )

           t(jkv+j,8) = 0.25*(
     .       (1.+bck(j+1,ii,2))*qk0(j+1,ii,2,3) -
     .       bck(j+1,ii,2)*q(j+1,k-1,ii,2) -
     .       (1.+bck(j-1,ii,2))*qk0(j-1,ii,2,3) +
     .       bck(j-1,ii,2)*q(j-1,k-1,ii,2) +
     .       q(j+1,k-1,ii,2) -
     .       q(j-1,k-1,ii,2) )

           t(jkv+j,9) = 0.25*(
     .       (1.+bck(j+1,ii,2))*qk0(j+1,ii,3,3) -
     .       bck(j+1,ii,2)*q(j+1,k-1,ii,3) -
     .       (1.+bck(j-1,ii,2))*qk0(j-1,ii,3,3) +
     .       bck(j-1,ii,2)*q(j-1,k-1,ii,3) +
     .       q(j+1,k-1,ii,3) -
     .       q(j-1,k-1,ii,3) )

           t(jkv+j,10)= 0.25*(
     .       (1.+bck(j+1,ii,2))*qk0(j+1,ii,4,3) -
     .       bck(j+1,ii,2)*q(j+1,k-1,ii,4) -
     .       (1.+bck(j-1,ii,2))*qk0(j-1,ii,4,3) +
     .       bck(j-1,ii,2)*q(j-1,k-1,ii,4) +
     .       q(j+1,k-1,ii,4) -
     .       q(j-1,k-1,ii,4) )
          end if
  227   continue
      else
c
c-- general case for k between 1 and kdim
c
        do 229 j = 1,jdim1
          if (j. eq. 1) then
c
c replace (j-1) data by boundary data 
c
           t(jkv+j,7) = 0.25*gprgm1*(
     .       q(j+1,k,ii,5)/q(j+1,k,ii,1) -
     .       (1.+bcj(k,ii,1))*qj0(k,ii,5,1)/qj0(k,ii,1,1) +
     .       bcj(k,ii,1)*q(j,k,ii,5)/q(j,k,ii,1) +
     .       q(j+1,k-1,ii,5)/q(j+1,k-1,ii,1) -
     .       (1.+bcj(k-1,ii,1))*qj0(k-1,ii,5,1)/qj0(k-1,ii,1,1) +
     .       bcj(k-1,ii,1)*q(j,k-1,ii,5)/q(j,k-1,ii,1) )

           t(jkv+j,8) = 0.25*(
     .       q(j+1,k,ii,2) -
     .       (1.+bcj(k,ii,1))*qj0(k,ii,2,1) +
     .       bcj(k,ii,1)*q(j,k,ii,2) +
     .       q(j+1,k-1,ii,2) -
     .       (1.+bcj(k-1,ii,1))*qj0(k-1,ii,2,1) +
     .       bcj(k-1,ii,1)*q(j,k-1,ii,2) )

           t(jkv+j,9) = 0.25*(
     .       q(j+1,k,ii,3) -
     .       (1.+bcj(k,ii,1))*qj0(k,ii,3,1) +
     .       bcj(k,ii,1)*q(j,k,ii,3) +
     .       q(j+1,k-1,ii,3) -
     .       (1.+bcj(k-1,ii,1))*qj0(k-1,ii,3,1) +
     .       bcj(k-1,ii,1)*q(j,k-1,ii,3) )

           t(jkv+j,10) = 0.25*(
     .       q(j+1,k,ii,4) -
     .       (1.+bcj(k,ii,1))*qj0(k,ii,4,1) +
     .       bcj(k,ii,1)*q(j,k,ii,4) +
     .       q(j+1,k-1,ii,4) -
     .       (1.+bcj(k-1,ii,1))*qj0(k-1,ii,4,1) +
     .       bcj(k-1,ii,1)*q(j,k-1,ii,4) )

          else if (j .eq. jdim1) then
c
c replace (j+1) data by boundary data 
c
           t(jkv+j,7) = 0.25*gprgm1*(
     .       (1.+bcj(k,ii,2))*qj0(k,ii,5,3)/qj0(k,ii,1,3) -
     .       bcj(k,ii,2)*q(j,k,ii,5)/q(j,k,ii,1) -
     .       q(j-1,k,ii,5)/q(j-1,k,ii,1)  +
     .       (1.+bcj(k-1,ii,2))*qj0(k-1,ii,5,3)/qj0(k-1,ii,1,3) -
     .       bcj(k-1,ii,2)*q(j,k-1,ii,5)/q(j,k-1,ii,1) -
     .       q(j-1,k-1,ii,5)/q(j-1,k-1,ii,1) )

           t(jkv+j,8) = 0.25*(
     .       (1.+bcj(k,ii,2))*qj0(k,ii,2,3) -
     .       bcj(k,ii,2)*q(j,k,ii,2) -
     .       q(j-1,k,ii,2)  +
     .       (1.+bcj(k-1,ii,2))*qj0(k-1,ii,2,3) -
     .       bcj(k-1,ii,2)*q(j,k-1,ii,2) -
     .       q(j-1,k-1,ii,2) )

           t(jkv+j,9) = 0.25*(
     .       (1.+bcj(k,ii,2))*qj0(k,ii,3,3) -
     .       bcj(k,ii,2)*q(j,k,ii,3) -
     .       q(j-1,k,ii,3)  +
     .       (1.+bcj(k-1,ii,2))*qj0(k-1,ii,3,3) -
     .       bcj(k-1,ii,2)*q(j,k-1,ii,3) -
     .       q(j-1,k-1,ii,3) )

           t(jkv+j,10) = 0.25*(
     .       (1.+bcj(k,ii,2))*qj0(k,ii,4,3) -
     .       bcj(k,ii,2)*q(j,k,ii,4) -
     .       q(j-1,k,ii,4)  +
     .       (1.+bcj(k-1,ii,2))*qj0(k-1,ii,4,3) -
     .       bcj(k-1,ii,2)*q(j,k-1,ii,4) -
     .       q(j-1,k-1,ii,4) )
          else
c
c-- general case
c
            t(jkv+j,7) = 0.25*gprgm1*(
     .        q(j+1,k,ii,5)/q(j+1,k,ii,1) - 
     .        q(j-1,k,ii,5)/q(j-1,k,ii,1) +
     .        q(j+1,k-1,ii,5)/q(j+1,k-1,ii,1) -
     .        q(j-1,k-1,ii,5)/q(j-1,k-1,ii,1) )

            t(jkv+j,8) = 0.25*(
     .        q(j+1,k,ii,2)   - q(j-1,k,ii,2) +
     .        q(j+1,k-1,ii,2) - q(j-1,k-1,ii,2) )

            t(jkv+j,9) = 0.25*(
     .        q(j+1,k,ii,3)   - q(j-1,k,ii,3) +
     .        q(j+1,k-1,ii,3) - q(j-1,k-1,ii,3) )

            t(jkv+j,10) = 0.25*(
     .        q(j+1,k,ii,4)   - q(j-1,k,ii,4) +
     .        q(j+1,k-1,ii,4) - q(j-1,k-1,ii,4) )

          end if
  229   continue
      end if
  230 continue
c******************************************************************************
c
c   gradients at cell interfaces (needed for full ns terms). Note that 
c  we don't need derivatives in zeta direction, because they have already
c  been accounted for in hfluxv.
c  Recall bci=0 => cell data, bci=1 => face data
c
c  t(3)       : d( c*c/(gm1*pr) )/d(xie)
c  t(4)       : d(u)/d(xie)
c  t(5)       : d(v)/d(xie)
c  t(6)       : d(w)/d(xie)
c
c     interior interfaces
c
      if (i2d .eq. 0 .and. idim .gt. 2) then
      do 380 ipl = 1,npl
      ii = i + ipl -1
      if (ii .eq. 1) then
c
c replace (ii-1) data by boundary data
c
        do 330 k = 1,kdim
          jkv     = (k-1)*npl*jdim1 + (ipl-1)*jdim1 
          if (k .eq. 1) then
c
c replace (k-1) data by boundary data
c for missing data at corners use one-sided differencing of known boundary values
c
            do 323 j = 1,jdim1

            t(jkv+j,3) = 0.25*gprgm1*(
     .        q(j,k,ii+1,5)/q(j,k,ii+1,1)  -
     .        (1.+bci(j,k,1))*qi0(j,k,5,1)/qi0(j,k,1,1) +
     .        bci(j,k,1)*q(j,k,ii,5)/q(j,k,ii,1) +
     .        2.*(1.+bck(j,ii+1,1))*qk0(j,ii+1,5,1)/qk0(j,ii+1,1,1) -
     .        2.*bck(j,ii+1,1)*q(j,k,ii+1,5)/q(j,k,ii+1,1) -
     .        2.*(1.+bck(j,ii,1))*qk0(j,ii,5,1)/qk0(j,ii,1,1) +
     .        2.*bck(j,ii,1)*q(j,k,ii,5)/q(j,k,ii,1) )

            t(jkv+j,4) = 0.25*(
     .        q(j,k,ii+1,2)  -
     .        (1.+bci(j,k,1))*qi0(j,k,2,1) +
     .        bci(j,k,1)*q(j,k,ii,2) +
     .        2.*(1.+bck(j,ii+1,1))*qk0(j,ii+1,2,1) -
     .        2.*bck(j,ii+1,1)*q(j,k,ii+1,2) -
     .        2.*(1.+bck(j,ii,1))*qk0(j,ii,2,1) +
     .        2.*bck(j,ii,1)*q(j,k,ii,2) )

            t(jkv+j,5) = 0.25*(
     .        q(j,k,ii+1,3)  -
     .        (1.+bci(j,k,1))*qi0(j,k,3,1) +
     .        bci(j,k,1)*q(j,k,ii,3) +
     .        2.*(1.+bck(j,ii+1,1))*qk0(j,ii+1,3,1) -
     .        2.*bck(j,ii+1,1)*q(j,k,ii+1,3) -
     .        2.*(1.+bck(j,ii,1))*qk0(j,ii,3,1) +
     .        2.*bck(j,ii,1)*q(j,k,ii,3) )

            t(jkv+j,6) = 0.25*(
     .        q(j,k,ii+1,4)  -
     .        (1.+bci(j,k,1))*qi0(j,k,4,1) +
     .        bci(j,k,1)*q(j,k,ii,4) +
     .        2.*(1.+bck(j,ii+1,1))*qk0(j,ii+1,4,1) -
     .        2.*bck(j,ii+1,1)*q(j,k,ii+1,4) -
     .        2.*(1.+bck(j,ii,1))*qk0(j,ii,4,1) +
     .        2.*bck(j,ii,1)*q(j,k,ii,4) )

  323       continue
          else if (k .eq. kdim) then
c
c replace k data by boundary data
c for missing data at corners use one-sided differencing of known boundary values
c
            do 325 j = 1,jdim1

            t(jkv+j,3) = 0.25*gprgm1*(
     .        2.*(1.+bck(j,ii+1,2))*qk0(j,ii+1,5,3)/qk0(j,ii+1,1,3) -
     .        2.*bck(j,ii+1,2)*q(j,k-1,ii+1,5)/q(j,k-1,ii+1,1) -
     .        2.*(1.+bck(j,ii,2))*qk0(j,ii,5,3)/qk0(j,ii,1,3) +
     .        2.*bck(j,ii,2)*q(j,k-1,ii,5)/q(j,k-1,ii,1) +
     .        q(j,k-1,ii+1,5)/q(j,k-1,ii+1,1) -
     .        (1.+bci(j,k-1,1))*qi0(j,k-1,5,1)/qi0(j,k-1,1,1) +
     .        bci(j,k-1,1)*q(j,k-1,ii,5)/q(j,k-1,ii,1) )

            t(jkv+j,4) = 0.25*(
     .        2.*(1.+bck(j,ii+1,2))*qk0(j,ii+1,2,3) -
     .        2.*bck(j,ii+1,2)*q(j,k-1,ii+1,2) -
     .        2.*(1.+bck(j,ii,2))*qk0(j,ii,2,3) +
     .        2.*bck(j,ii,2)*q(j,k-1,ii,2) +
     .        q(j,k-1,ii+1,2) -
     .        (1.+bci(j,k-1,1))*qi0(j,k-1,2,1) +
     .        bci(j,k-1,1)*q(j,k-1,ii,2) )

            t(jkv+j,5) = 0.25*(
     .        2.*(1.+bck(j,ii+1,2))*qk0(j,ii+1,3,3) -
     .        2.*bck(j,ii+1,2)*q(j,k-1,ii+1,3) -
     .        2.*(1.+bck(j,ii,2))*qk0(j,ii,3,3) +
     .        2.*bck(j,ii,2)*q(j,k-1,ii,3) +
     .        q(j,k-1,ii+1,3) -
     .        (1.+bci(j,k-1,1))*qi0(j,k-1,3,1) +
     .        bci(j,k-1,1)*q(j,k-1,ii,3) )

            t(jkv+j,6) = 0.25*(
     .        2.*(1.+bck(j,ii+1,2))*qk0(j,ii+1,4,3) -
     .        2.*bck(j,ii+1,2)*q(j,k-1,ii+1,4) -
     .        2.*(1.+bck(j,ii,2))*qk0(j,ii,4,3) +
     .        2.*bck(j,ii,2)*q(j,k-1,ii,4) +
     .        q(j,k-1,ii+1,4) -
     .        (1.+bci(j,k-1,1))*qi0(j,k-1,4,1) +
     .        bci(j,k-1,1)*q(j,k-1,ii,4) )

  325       continue
          else
c
c-- general case for ii = 1
c
            do 327 j = 1,jdim1

            t(jkv+j,3) = 0.25*gprgm1*(
     .        q(j,k,ii+1,5)/q(j,k,ii+1,1) -
     .        (1.+bci(j,k,1))*qi0(j,k,5,1)/qi0(j,k,1,1) +
     .        bci(j,k,1)*q(j,k,ii,5)/q(j,k,ii,1) +
     .        q(j,k-1,ii+1,5)/q(j,k-1,ii+1,1) -
     .        (1.+bci(j,k-1,1))*qi0(j,k-1,5,1)/qi0(j,k-1,1,1) +
     .        bci(j,k-1,1)*q(j,k-1,ii,5)/q(j,k-1,ii,1) )

            t(jkv+j,4) = 0.25*(
     .        q(j,k,ii+1,2) -
     .        (1.+bci(j,k,1))*qi0(j,k,2,1) +
     .        bci(j,k,1)*q(j,k,ii,2) +
     .        q(j,k-1,ii+1,2) -
     .        (1.+bci(j,k-1,1))*qi0(j,k-1,2,1) +
     .        bci(j,k-1,1)*q(j,k-1,ii,2) )

            t(jkv+j,5) = 0.25*(
     .        q(j,k,ii+1,3) -
     .        (1.+bci(j,k,1))*qi0(j,k,3,1) +
     .        bci(j,k,1)*q(j,k,ii,3) +
     .        q(j,k-1,ii+1,3) -
     .        (1.+bci(j,k-1,1))*qi0(j,k-1,3,1) +
     .        bci(j,k-1,1)*q(j,k-1,ii,3) )

            t(jkv+j,6) = 0.25*(
     .        q(j,k,ii+1,4) -
     .        (1.+bci(j,k,1))*qi0(j,k,4,1) +
     .        bci(j,k,1)*q(j,k,ii,4) +
     .        q(j,k-1,ii+1,4) -
     .        (1.+bci(j,k-1,1))*qi0(j,k-1,4,1) +
     .        bci(j,k-1,1)*q(j,k-1,ii,4) )

  327       continue
          end if
  330   continue

      else if (ii .eq. idim1) then
c
c replace (ii+1) data by boundary data
c
        do 340 k = 1,kdim
          jkv     = (k-1)*npl*jdim1 + (ipl-1)*jdim1 
          if (k .eq. 1) then
c
c replace (k-1) data by boundary data
c for missing data at corners use one-sided differencing of known boundary values
c
            do 333 j = 1,jdim1

            t(jkv+j,3) = 0.25*gprgm1*(
     .        (1.+bci(j,k,2))*qi0(j,k,5,3)/qi0(j,k,1,3) -
     .        bci(j,k,2)*q(j,k,ii,5)/q(j,k,ii,1) -
     .        q(j,k,ii-1,5)/q(j,k,ii-1,1) +
     .        2.*(1.+bck(j,ii,1))*qk0(j,ii,5,1)/qk0(j,ii,1,1) -
     .        2.*bck(j,ii,1)*q(j,k,ii,5)/q(j,k,ii,1) -
     .        2.*(1.+bck(j,ii-1,1))*qk0(j,ii-1,5,1)/qk0(j,ii-1,1,1) +
     .        2.*bck(j,ii-1,1)*q(j,k,ii-1,5)/q(j,k,ii-1,1) )

            t(jkv+j,4) = 0.25*(
     .        (1.+bci(j,k,2))*qi0(j,k,2,3) -
     .        bci(j,k,2)*q(j,k,ii,2) -
     .        q(j,k,ii-1,2) +
     .        2.*(1.+bck(j,ii,1))*qk0(j,ii,2,1) -
     .        2.*bck(j,ii,1)*q(j,k,ii,2) -
     .        2.*(1.+bck(j,ii-1,1))*qk0(j,ii-1,2,1) +
     .        2.*bck(j,ii-1,1)*q(j,k,ii-1,2) )

            t(jkv+j,5) = 0.25*(
     .        (1.+bci(j,k,2))*qi0(j,k,3,3) -
     .        bci(j,k,2)*q(j,k,ii,3) -
     .        q(j,k,ii-1,3) +
     .        2.*(1.+bck(j,ii,1))*qk0(j,ii,3,1) -
     .        2.*bck(j,ii,1)*q(j,k,ii,3) -
     .        2.*(1.+bck(j,ii-1,1))*qk0(j,ii-1,3,1) +
     .        2.*bck(j,ii-1,1)*q(j,k,ii-1,3) )

            t(jkv+j,6) = 0.25*(
     .        (1.+bci(j,k,2))*qi0(j,k,4,3) -
     .        bci(j,k,2)*q(j,k,ii,4) -
     .        q(j,k,ii-1,4) +
     .        2.*(1.+bck(j,ii,1))*qk0(j,ii,4,1) -
     .        2.*bck(j,ii,1)*q(j,k,ii,4) -
     .        2.*(1.+bck(j,ii-1,1))*qk0(j,ii-1,4,1) +
     .        2.*bck(j,ii-1,1)*q(j,k,ii-1,4) )

  333       continue
          else if (k .eq. kdim) then
c
c replace k data by boundary data
c for missing data at corners use one-sided differencing of known boundary values
c
            do 335 j = 1,jdim1

            t(jkv+j,3) = 0.25*gprgm1*(
     .        2.*(1.+bck(j,ii,2))*qk0(j,ii,5,3)/qk0(j,ii,1,3) -
     .        2.*bck(j,ii,2)*q(j,k-1,ii,5)/q(j,k-1,ii,1) -
     .        2.*(1.+bck(j,ii-1,2))*qk0(j,ii-1,5,3)/qk0(j,ii-1,1,3) +
     .        2.*bck(j,ii-1,2)*q(j,k-1,ii-1,5)/q(j,k-1,ii-1,1) +
     .        (1.+bci(j,k-1,2))*qi0(j,k-1,5,3)/qi0(j,k-1,1,3) -
     .        bci(j,k-1,2)*q(j,k-1,ii,5)/q(j,k-1,ii,1) -
     .        q(j,k-1,ii-1,5)/q(j,k-1,ii-1,1) )

            t(jkv+j,4) = 0.25*(
     .        2.*(1.+bck(j,ii,2))*qk0(j,ii,2,3) -
     .        2.*bck(j,ii,2)*q(j,k-1,ii,2) -
     .        2.*(1.+bck(j,ii-1,2))*qk0(j,ii-1,2,3) +
     .        2.*bck(j,ii-1,2)*q(j,k-1,ii-1,2) +
     .        (1.+bci(j,k-1,2))*qi0(j,k-1,2,3) -
     .        bci(j,k-1,2)*q(j,k-1,ii,2) -
     .        q(j,k-1,ii-1,2) )

            t(jkv+j,5) = 0.25*(
     .        2.*(1.+bck(j,ii,2))*qk0(j,ii,3,3) -
     .        2.*bck(j,ii,2)*q(j,k-1,ii,3) -
     .        2.*(1.+bck(j,ii-1,2))*qk0(j,ii-1,3,3) +
     .        2.*bck(j,ii-1,2)*q(j,k-1,ii-1,3) +
     .        (1.+bci(j,k-1,2))*qi0(j,k-1,3,3) -
     .        bci(j,k-1,2)*q(j,k-1,ii,3) -
     .        q(j,k-1,ii-1,3) )

            t(jkv+j,6) = 0.25*(
     .        2.*(1.+bck(j,ii,2))*qk0(j,ii,4,3) -
     .        2.*bck(j,ii,2)*q(j,k-1,ii,4) -
     .        2.*(1.+bck(j,ii-1,2))*qk0(j,ii-1,4,3) +
     .        2.*bck(j,ii-1,2)*q(j,k-1,ii-1,4) +
     .        (1.+bci(j,k-1,2))*qi0(j,k-1,4,3) -
     .        bci(j,k-1,2)*q(j,k-1,ii,4) -
     .        q(j,k-1,ii-1,4) )

  335       continue
          else
c
c-- general case for ii = idim1
c
            do 337 j = 1,jdim1

            t(jkv+j,3) = 0.25*gprgm1*(
     .        (1.+bci(j,k,2))*qi0(j,k,5,3)/qi0(j,k,1,3) -
     .        bci(j,k,2)*q(j,k,ii,5)/q(j,k,ii,1) -
     .        q(j,k,ii-1,5)/q(j,k,ii-1,1) +
     .        (1.+bci(j,k-1,2))*qi0(j,k-1,5,3)/qi0(j,k-1,1,3) -
     .        bci(j,k-1,2)*q(j,k-1,ii,5)/q(j,k-1,ii,1) -
     .        q(j,k-1,ii-1,5)/q(j,k-1,ii-1,1) )

            t(jkv+j,4) = 0.25*(
     .        (1.+bci(j,k,2))*qi0(j,k,2,3) -
     .        bci(j,k,2)*q(j,k,ii,2) -
     .        q(j,k,ii-1,2) +
     .        (1.+bci(j,k-1,2))*qi0(j,k-1,2,3) -
     .        bci(j,k-1,2)*q(j,k-1,ii,2) -
     .        q(j,k-1,ii-1,2) )

            t(jkv+j,5) = 0.25*(
     .        (1.+bci(j,k,2))*qi0(j,k,3,3) -
     .        bci(j,k,2)*q(j,k,ii,3) -
     .        q(j,k,ii-1,3) +
     .        (1.+bci(j,k-1,2))*qi0(j,k-1,3,3) -
     .        bci(j,k-1,2)*q(j,k-1,ii,3) -
     .        q(j,k-1,ii-1,3) )

            t(jkv+j,6) = 0.25*(
     .        (1.+bci(j,k,2))*qi0(j,k,4,3) -
     .        bci(j,k,2)*q(j,k,ii,4) -
     .        q(j,k,ii-1,4) +
     .        (1.+bci(j,k-1,2))*qi0(j,k-1,4,3) -
     .        bci(j,k-1,2)*q(j,k-1,ii,4) -
     .        q(j,k-1,ii-1,4) )

  337       continue
          end if
  340   continue
      else
c
c-- general case (ii .ne 1. and ii .ne. idim1)
c
        do 350 k = 1,kdim
          jkv     = (k-1)*npl*jdim1 + (ipl-1)*jdim1 
          if (k .eq. 1) then
c
c replace (k-1) data by boundary data
c
            do 343 j = 1,jdim1

            t(jkv+j,3) = 0.25*gprgm1*(
     .        q(j,k,ii+1,5)/q(j,k,ii+1,1) -
     .        q(j,k,ii-1,5)/q(j,k,ii-1,1) +
     .        (1.+bck(j,ii+1,1))*qk0(j,ii+1,5,1)/qk0(j,ii+1,1,1) -
     .        bck(j,ii+1,1)*q(j,k,ii+1,5)/q(j,k,ii+1,1) -
     .        (1.+bck(j,ii-1,1))*qk0(j,ii-1,5,1)/qk0(j,ii-1,1,1) +
     .        bck(j,ii-1,1)*q(j,k,ii-1,5)/q(j,k,ii-1,1) )

            t(jkv+j,4) = 0.25*(
     .        q(j,k,ii+1,2) -q(j,k,ii-1,2) +
     .        (1.+bck(j,ii+1,1))*qk0(j,ii+1,2,1) -
     .        bck(j,ii+1,1)*q(j,k,ii+1,2) -
     .        (1.+bck(j,ii-1,1))*qk0(j,ii-1,2,1) +
     .        bck(j,ii-1,1)*q(j,k,ii-1,2) )

            t(jkv+j,5) = 0.25*(
     .        q(j,k,ii+1,3) -q(j,k,ii-1,3) +
     .        (1.+bck(j,ii+1,1))*qk0(j,ii+1,3,1) -
     .        bck(j,ii+1,1)*q(j,k,ii+1,3) -
     .        (1.+bck(j,ii-1,1))*qk0(j,ii-1,3,1) +
     .        bck(j,ii-1,1)*q(j,k,ii-1,3) )

            t(jkv+j,6) = 0.25*(
     .        q(j,k,ii+1,4) -q(j,k,ii-1,4) +
     .        (1.+bck(j,ii+1,1))*qk0(j,ii+1,4,1) -
     .        bck(j,ii+1,1)*q(j,k,ii+1,4) -
     .        (1.+bck(j,ii-1,1))*qk0(j,ii-1,4,1) +
     .        bck(j,ii-1,1)*q(j,k,ii-1,4) )

  343       continue
          else if (k .eq. kdim) then
c
c replace k data by boundary data
c
            do 345 j = 1,jdim1

            t(jkv+j,3) = 0.25*gprgm1*(
     .        (1.+bck(j,ii+1,2))*qk0(j,ii+1,5,3)/qk0(j,ii+1,1,3) -
     .        bck(j,ii+1,2)*q(j,k-1,ii+1,5)/q(j,k-1,ii+1,1) -
     .        (1.+bck(j,ii-1,2))*qk0(j,ii-1,5,3)/qk0(j,ii-1,1,3) +
     .        bck(j,ii-1,2)*q(j,k-1,ii-1,5)/q(j,k-1,ii-1,1) +
     .        q(j,k-1,ii+1,5)/q(j,k-1,ii+1,1) -
     .        q(j,k-1,ii-1,5)/q(j,k-1,ii-1,1) )

            t(jkv+j,4) = 0.25*(
     .        (1.+bck(j,ii+1,2))*qk0(j,ii+1,2,3) -
     .        bck(j,ii+1,2)*q(j,k-1,ii+1,2) -
     .        (1.+bck(j,ii-1,2))*qk0(j,ii-1,2,3) +
     .        bck(j,ii-1,2)*q(j,k-1,ii-1,2) +
     .        q(j,k-1,ii+1,2) -q(j,k-1,ii-1,2) )

            t(jkv+j,5) = 0.25*(
     .        (1.+bck(j,ii+1,2))*qk0(j,ii+1,3,3) -
     .        bck(j,ii+1,2)*q(j,k-1,ii+1,3) -
     .        (1.+bck(j,ii-1,2))*qk0(j,ii-1,3,3) +
     .        bck(j,ii-1,2)*q(j,k-1,ii-1,3) +
     .        q(j,k-1,ii+1,3) -q(j,k-1,ii-1,3) )

            t(jkv+j,6) = 0.25*(
     .        (1.+bck(j,ii+1,2))*qk0(j,ii+1,4,3) -
     .        bck(j,ii+1,2)*q(j,k-1,ii+1,4) -
     .        (1.+bck(j,ii-1,2))*qk0(j,ii-1,4,3) +
     .        bck(j,ii-1,2)*q(j,k-1,ii-1,4) +
     .        q(j,k-1,ii+1,4) -q(j,k-1,ii-1,4) )

  345       continue
          else
            do 347 j = 1,jdim1
c
c-- general case
c
            t(jkv+j,3) = 0.25*gprgm1*(
     .        q(j,k,ii+1,5)/q(j,k,ii+1,1) -
     .        q(j,k,ii-1,5)/q(j,k,ii-1,1) +
     .        q(j,k-1,ii+1,5)/q(j,k-1,ii+1,1) -
     .        q(j,k-1,ii-1,5)/q(j,k-1,ii-1,1) )

            t(jkv+j,4) = 0.25*(
     .        q(j,k,ii+1,2)  -q(j,k,ii-1,2) +
     .        q(j,k-1,ii+1,2)-q(j,k-1,ii-1,2) )

            t(jkv+j,5) = 0.25*(
     .        q(j,k,ii+1,3)  -q(j,k,ii-1,3) +
     .        q(j,k-1,ii+1,3)-q(j,k-1,ii-1,3) )

            t(jkv+j,6) = 0.25*(
     .        q(j,k,ii+1,4)  -q(j,k,ii-1,4) +
     .        q(j,k-1,ii+1,4)-q(j,k-1,ii-1,4) )

  347       continue
          end if
  350   continue
      end if
  380 continue
      end if
c
c******************************************************************************
c
c  t(24) : turbulent viscosity at interfaces (=0 for laminar flow)
c
      iviscc = ivisc(3)
      if (iviscc.gt.1) then
      do 1401 ipl=1,npl
      ii = i+ipl-1
      do 1401 k=1,kdim
      jkv=(k-1)*npl*jdim1 + (ipl-1)*jdim1
c     k=0 interface
      if (k.eq.1) then
        do 1017 j=1,jdim1
          t(j+jkv,24)=bck(j,ii,1)*vk0(j,ii,1,1)+
     +     (1.-bck(j,ii,1))*0.5*(vk0(j,ii,1,1)+vist3d(j,1,ii))
 1017   continue
c     k=kdim interface
      else if (k.eq.kdim) then
        do 1117 j=1,jdim1
          t(j+jkv,24)=bck(j,ii,2)*vk0(j,ii,1,3)+
     +      (1.-bck(j,ii,2))*0.5*(vk0(j,ii,1,3)+vist3d(j,kdim1,ii))
 1117   continue
      else
c     interior interfaces
        do 1217 j=1,jdim1
          t(j+jkv,24)=0.5*(vist3d(j,k,ii)+vist3d(j,k-1,ii))
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
      t25 = t(izz,24)/t(izz,14) 
      t24 = 1.e0+t25*coef_eddy
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
 1022 continue
c
      if (iadv .lt. 0) return
c
c******************************************************************************
c
c  t(11) : u at cell interfaces
c  t(12) : v at cell interfaces
c  t(13) : w at cell interfaces
c
c     interior interfaces
cdir$ ivdep
      do 1024 izz=1,nn
      t(izz+js-1,11) =0.5e0*(t(izz+js-1,18)+t(izz,18))
      t(izz+js-1,12) =0.5e0*(t(izz+js-1,19)+t(izz,19))
      t(izz+js-1,13) =0.5e0*(t(izz+js-1,20)+t(izz,20))
 1024 continue
c     interfaces at j=0 and j=jdim
cdir$ ivdep
      do 1025 izz=1,jv
      ab=1.+wk0(izz,15)
      bb=1.-wk0(izz,15)
      t(izz+n,11)=0.5*(ab*wk0(izz,18)+bb*t(izz+n-jv,18))
      t(izz+n,12)=0.5*(ab*wk0(izz,19)+bb*t(izz+n-jv,19))
      t(izz+n,13)=0.5*(ab*wk0(izz,20)+bb*t(izz+n-jv,20))
      ab=1.+wk0(izz,5)
      bb=1.-wk0(izz,5)
      t(izz,11) = 0.5*(ab*wk0(izz,8) +bb*t(izz,18))
      t(izz,12) = 0.5*(ab*wk0(izz,9) +bb*t(izz,19))
      t(izz,13) = 0.5*(ab*wk0(izz,10) +bb*t(izz,20))
 1025 continue
c
c******************************************************************************
c
c  calculate fluxes
c
c  viscous terms at interfaces
cdir$ ivdep
      do 1026 izz=1,l0
c
c-- form ux, uy, uz, vx, vy, vz, wx, wy, wz, tx,ty and tz 
c   excluding contributions in k/zeta direction (already accounted for in 
c   thin-layer approximation. Note that t = c*c/(gm1*pr)
c
       ux = t(izz,4)*t(izz,28) + t(izz,8)*t(izz,31)
       uy = t(izz,4)*t(izz,29) + t(izz,8)*t(izz,32)
       uz = t(izz,4)*t(izz,30) + t(izz,8)*t(izz,2)

       vx = t(izz,5)*t(izz,28) + t(izz,9)*t(izz,31)
       vy = t(izz,5)*t(izz,29) + t(izz,9)*t(izz,32)
       vz = t(izz,5)*t(izz,30) + t(izz,9)*t(izz,2)

       wx = t(izz,6)*t(izz,28) + t(izz,10)*t(izz,31)
       wy = t(izz,6)*t(izz,29) + t(izz,10)*t(izz,32)
       wz = t(izz,6)*t(izz,30) + t(izz,10)*t(izz,2)

       tx = t(izz,3)*t(izz,28) + t(izz,7)*t(izz,31)
       ty = t(izz,3)*t(izz,29) + t(izz,7)*t(izz,32)
       tz = t(izz,3)*t(izz,30) + t(izz,7)*t(izz,2)
c
c-- Note lambda = -2/3 mu, lambda + 2 mu = 4/3 mu 
c--  u-momentum flux
c 
      t(izz,16)      = t(izz,14)*(  t(izz,25)*(4./3.*ux-2./3.*(vy + wz))
     .                            + t(izz,26)*(uy + vx) 
     .                            + t(izz,27)*(uz + wx)) 

c
c--  v-momentum flux
c 
      t(izz,17)      = t(izz,14)*(  t(izz,25)*(vx + uy)
     .                            + t(izz,26)*(4./3.*vy-2./3.*(ux + wz))
     .                            + t(izz,27)*(vz + wy)) 

c
c--  w-momentum flux
c 
      t(izz,18)      = t(izz,14)*(  t(izz,25)*(wx+uz)
     .                            + t(izz,26)*(wy+vz)
     .                            + t(izz,27)*(4./3.*wz-2./3.*(vy+ux)))
c
c--  energy flux
c 
      t(izz,19) = t(izz,14)*( 
     .              t(izz,25)*((4./3.*ux - 2./3*(vy+wz))*t(izz,11) +
     .                         (uy+vx)*t(izz,12) + (uz+wx)*t(izz,13)) +
     .              t(izz,26)*((vx+uy)*t(izz,11) + (vz+wy)*t(izz,13) +
     .                         (4./3.*vy - 2./3*(ux+wz))*t(izz,12)) +
     .              t(izz,27)*((wx+uz)*t(izz,11) + (wy+vz)*t(izz,12) +
     .                         (4./3.*wz - 2./3*(ux+vy))*t(izz,13))  +
     .          t(izz,23)*(t(izz,25)*tx + t(izz,26)*ty + t(izz,27)*tz))
c   
 1026 continue
c  (-)viscous flux =  Gv(k-1/2)-Gv(k+1/2) (momentum eqs.)
cdir$ ivdep
      do 1027 izz=1,n
      t(izz,2)       = -t(izz+js-1,16)+t(izz,16)
      t(izz,3)       = -t(izz+js-1,17)+t(izz,17)
      t(izz,4)       = -t(izz+js-1,18)+t(izz,18)
      t(izz,5)       = -t(izz+js-1,19)+t(izz,19)
 1027 continue
c
c******************************************************************************
c
c  calculate residuals
c
      do 400 ipl=1,npl
      ii  = i+ipl-1
      do 400 k=1,kdim1
      jkv = (k-1)*jdim1*npl + (ipl-1)*jdim1 + 1
c   for 2-D, t(3) should be identically zero
      if (i2d .eq. 1) then
        do izz=1,jdim1
          t(izz+jkv-1,3)=0.
        enddo
      end if
      do 400 l=2,5
cdir$ ivdep
      do 1030 izz=1,jdim1
      res(izz,k,ii,l) = res(izz,k,ii,l)+t(izz+jkv-1,l)
 1030 continue
  400 continue
      return
      end
