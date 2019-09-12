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
      subroutine gfluxv1(i,npl,jdim,kdim,idim,res,q,qj0,qk0,qi0,
     .                  sj,sk,si,vol,t,nvtq,wj0,vist3d,vj0,
     .                  bcj,bck,bci,volj0,
     .                  nou,bou,nbuf,ibufdim,iadv)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Calculate right-hand residual contributions in the 
c     J-direction due to the viscous terms missing in the thin-layer 
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
      dimension sj(jdim,kdim,idim-1,5),vol(jdim,kdim,idim-1)
      dimension si(jdim,kdim,idim,5),sk(jdim,kdim,idim-1,5)

      dimension q(jdim,kdim,idim,5)
      dimension qj0(kdim,idim-1,5,4), qk0(jdim,idim-1,5,4),
     .          qi0(jdim,kdim,5,4)
      dimension res(jdim,kdim,idim-1,5),vist3d(jdim,kdim,idim),
     .          vj0(kdim,idim-1,1,4),volj0(kdim,idim-1,4)
      dimension t(nvtq,32),wj0(npl*(kdim-1),22)
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
      if(ivisc(2)>=70) coef_eddy = 0.0
      if(i_tauijs .eq. 1) coef_eddy=0.
      kdim1 = kdim-1
      jdim1 = jdim-1
      idim1 = idim-1
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
      gprgm1 = gpr/gm1
c
      if (isklton.gt.0 .and. i.eq.1 .and. iadv.ge.0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1442)
      end if
 1442 format(52h   computing cross-derivative viscous fluxes, J-dir.)
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
       t1          =  sj(j,k,ii,4)*t(k+jkv,15) 
       t(k+jkv,25) =  sj(j,k,ii,1)*t1
       t(k+jkv,26) =  sj(j,k,ii,2)*t1
       t(k+jkv,27) =  sj(j,k,ii,3)*t1
 1014 continue
  200 continue
c
c t(28)  - t(30) : Components of grad (zeta)
c grad(zeta) = 1/4 [(sk(j,k,i,l)*sk(j,k,i,4) + sk(j-1,k,i,l)*sk(j-1,k,i,4) + 
c sk(j,k+1,i,l)*sk(j,k+1,i,4) + sk(j-1,k+1,i,l)*sk(j-1,k+1,i,4)] /V_int,
c j=1,jdim , k=1,kdim1, i=1,npl, l =1,3.  At the boundaries in j, 
c since we don't have sk's, we assume that them to be equal to the interior one.
c
c Note: the code for j=1 is obtained from the general case by substituting
c j for j-1, and likewise the code for j=jdim by substituting j-1 for j.
c
      do 210 ipl=1,npl
      ii  = i+ipl-1
      do 210 j=1,jdim
      jkv  = (j-1)*npl*kdim1 + (ipl-1)*kdim1 
      if (j.eq.1) then 
       do 204 k = 1,kdim1
        t(jkv+k,28) = 0.25*t(k+jkv,15)*(sk(j,k,ii,1)*sk(j,k,ii,4)
     .              + sk(j,k,ii,1)*sk(j,k,ii,4)
     .              + sk(j,k+1,ii,1)*sk(j,k+1,ii,4)
     .              + sk(j,k+1,ii,1)*sk(j,k+1,ii,4))
        t(jkv+k,29) = 0.25*t(k+jkv,15)*(sk(j,k,ii,2)*sk(j,k,ii,4)
     .              + sk(j,k,ii,2)*sk(j,k,ii,4)
     .              + sk(j,k+1,ii,2)*sk(j,k+1,ii,4)
     .              + sk(j,k+1,ii,2)*sk(j,k+1,ii,4))
        t(jkv+k,30) = 0.25*t(k+jkv,15)*(sk(j,k,ii,3)*sk(j,k,ii,4)
     .              + sk(j,k,ii,3)*sk(j,k,ii,4)
     .              + sk(j,k+1,ii,3)*sk(j,k+1,ii,4)
     .              + sk(j,k+1,ii,3)*sk(j,k+1,ii,4))
  204  continue
      else if (j.eq.jdim) then
       do 206 k = 1,kdim1
        t(jkv+k,28) = 0.25*t(k+jkv,15)*(sk(j-1,k,ii,1)*sk(j-1,k,ii,4)
     .              + sk(j-1,k,ii,1)*sk(j-1,k,ii,4)
     .              + sk(j-1,k+1,ii,1)*sk(j-1,k+1,ii,4)
     .              + sk(j-1,k+1,ii,1)*sk(j-1,k+1,ii,4))
        t(jkv+k,29) = 0.25*t(k+jkv,15)*(sk(j-1,k,ii,2)*sk(j-1,k,ii,4)
     .              + sk(j-1,k,ii,2)*sk(j-1,k,ii,4)
     .              + sk(j-1,k+1,ii,2)*sk(j-1,k+1,ii,4)
     .              + sk(j-1,k+1,ii,2)*sk(j-1,k+1,ii,4))
        t(jkv+k,30) = 0.25*t(k+jkv,15)*(sk(j-1,k,ii,3)*sk(j-1,k,ii,4)
     .              + sk(j-1,k,ii,3)*sk(j-1,k,ii,4)
     .              + sk(j-1,k+1,ii,3)*sk(j-1,k+1,ii,4)
     .              + sk(j-1,k+1,ii,3)*sk(j-1,k+1,ii,4))
  206  continue
      else
c
c-- general case
c
       do 208 k = 1,kdim1
        t(jkv+k,28) = 0.25*t(k+jkv,15)*(sk(j,k,ii,1)*sk(j,k,ii,4)
     .              + sk(j-1,k,ii,1)*sk(j-1,k,ii,4)
     .              + sk(j,k+1,ii,1)*sk(j,k+1,ii,4)
     .              + sk(j-1,k+1,ii,1)*sk(j-1,k+1,ii,4))
        t(jkv+k,29) = 0.25*t(k+jkv,15)*(sk(j,k,ii,2)*sk(j,k,ii,4)
     .              + sk(j-1,k,ii,2)*sk(j-1,k,ii,4)
     .              + sk(j,k+1,ii,2)*sk(j,k+1,ii,4)
     .              + sk(j-1,k+1,ii,2)*sk(j-1,k+1,ii,4))
        t(jkv+k,30) = 0.25*t(k+jkv,15)*(sk(j,k,ii,3)*sk(j,k,ii,4)
     .              + sk(j-1,k,ii,3)*sk(j-1,k,ii,4)
     .              + sk(j,k+1,ii,3)*sk(j,k+1,ii,4)
     .              + sk(j-1,k+1,ii,3)*sk(j-1,k+1,ii,4))
  208  continue
      end if
  210 continue
c
c t(31),t(32),t(2)  : Components of grad (xie)
c grad(xie) = 1/4 [(si(j,k,i,l)*si(j,k,i,4) + si(j-1,k,i,l)*si(j-1,k,i,4) +
c si(j,k,i+1,l)*si(j,k,i+1,4) + si(j-1,k,i+1,l)*si(j-1,k,i+1,4)] /V_int,
c j=1,jdim , k=1,kdim1, i=1,npl,l = 1,3.  At the boundaries in j, 
c since we don't have si's, we assume that them to be equal to the interior one.
c
c Note: the code for j=1 is obtained from the general case by substituting
c j for j-1, and likewise the code for j=jdim by substituting j-1 for j.
c
      do 220 ipl=1,npl
      ii  = i+ipl-1
      do 220 j=1,jdim
      jkv  = (j-1)*npl*kdim1 + (ipl-1)*kdim1 
      if (j.eq.1) then 
       do 214 k = 1,kdim1
        t(jkv+k,31) = 0.25*t(k+jkv,15)*(si(j,k,ii,1)*si(j,k,ii,4)
     .              + si(j,k,ii,1)*si(j,k,ii,4)
     .              + si(j,k,ii+1,1)*si(j,k,ii+1,4)
     .              + si(j,k,ii+1,1)*si(j,k,ii+1,4))
        t(jkv+k,32) = 0.25*t(k+jkv,15)*(si(j,k,ii,2)*si(j,k,ii,4)
     .              + si(j,k,ii,2)*si(j,k,ii,4)
     .              + si(j,k,ii+1,2)*si(j,k,ii+1,4)
     .              + si(j,k,ii+1,2)*si(j,k,ii+1,4))
        t(jkv+k,2)  = 0.25*t(k+jkv,15)*(si(j,k,ii,3)*si(j,k,ii,4)
     .              + si(j,k,ii,3)*si(j,k,ii,4)
     .              + si(j,k,ii+1,3)*si(j,k,ii+1,4)
     .              + si(j,k,ii+1,3)*si(j,k,ii+1,4))
  214  continue
      else if (j.eq.jdim) then
       do 216 k = 1,kdim1
        t(jkv+k,31) = 0.25*t(k+jkv,15)*(si(j-1,k,ii,1)*si(j-1,k,ii,4)
     .              + si(j-1,k,ii,1)*si(j-1,k,ii,4)
     .              + si(j-1,k,ii+1,1)*si(j-1,k,ii+1,4)
     .              + si(j-1,k,ii+1,1)*si(j-1,k,ii+1,4))
        t(jkv+k,32) = 0.25*t(k+jkv,15)*(si(j-1,k,ii,2)*si(j-1,k,ii,4)
     .              + si(j-1,k,ii,2)*si(j-1,k,ii,4)
     .              + si(j-1,k,ii+1,2)*si(j-1,k,ii+1,4)
     .              + si(j-1,k,ii+1,2)*si(j-1,k,ii+1,4))
        t(jkv+k,2)  = 0.25*t(k+jkv,15)*(si(j-1,k,ii,3)*si(j-1,k,ii,4)
     .              + si(j-1,k,ii,3)*si(j-1,k,ii,4)
     .              + si(j-1,k,ii+1,3)*si(j-1,k,ii+1,4)
     .              + si(j-1,k,ii+1,3)*si(j-1,k,ii+1,4))
  216  continue
      else
c
c-- general case
c
       do 218 k = 1,kdim1
        t(jkv+k,31) = 0.25*t(k+jkv,15)*(si(j,k,ii,1)*si(j,k,ii,4)
     .              + si(j-1,k,ii,1)*si(j-1,k,ii,4)
     .              + si(j,k,ii+1,1)*si(j,k,ii+1,4)
     .              + si(j-1,k,ii+1,1)*si(j-1,k,ii+1,4))
        t(jkv+k,32) = 0.25*t(k+jkv,15)*(si(j,k,ii,2)*si(j,k,ii,4)
     .              + si(j-1,k,ii,2)*si(j-1,k,ii,4)
     .              + si(j,k,ii+1,2)*si(j,k,ii+1,4)
     .              + si(j-1,k,ii+1,2)*si(j-1,k,ii+1,4))
        t(jkv+k,2)  = 0.25*t(k+jkv,15)*(si(j,k,ii,3)*si(j,k,ii,4)
     .              + si(j-1,k,ii,3)*si(j-1,k,ii,4)
     .              + si(j,k,ii+1,3)*si(j,k,ii+1,4)
     .              + si(j-1,k,ii+1,3)*si(j-1,k,ii+1,4))
  218  continue
      end if
  220 continue
c
c
c******************************************************************************
c
c   gradients at cell interfaces (needed for full ns terms). Note that 
c  we don't need derivatives in eta direction, because they have already
c  been accounted for in gfluxv.
c  Recall bck=0 => cell data, bck=1 => face data
c
c  t(3)       : d( c*c/(gm1*pr) )/d(zeta)
c  t(4)       : d(u)/d(zeta)
c  t(5)       : d(v)/d(zeta)
c  t(6)       : d(w)/d(zeta)
c
c     interior interfaces
c
c
      do 230 ipl = 1,npl
      ii = i + ipl -1
      do 230 j = 1,jdim
      jkv  = (j-1)*npl*kdim1 + (ipl-1)*kdim1
c
      if (j.eq.1) then
c
c replace (j-1) data by boundary data 
c for missing data at corners use one-sided differencing of known boundary values
c
        do 225 k = 1,kdim1
          if (k. eq. 1) then
c
c replace (k-1) data by boundary data 
c
           t(jkv+k,3) = 0.25*gprgm1*(
     .       q(j,k+1,ii,5)/q(j,k+1,ii,1) -
     .       (1.+bck(j,ii,1))*qk0(j,ii,5,1)/qk0(j,ii,1,1) +
     .       bck(j,ii,1)*q(j,k,ii,5)/q(j,k,ii,1) +
     .       2.*(1.+bcj(k+1,ii,1))*qj0(k+1,ii,5,1)/qj0(k+1,ii,1,1) -
     .       2.*bcj(k+1,ii,1)*q(j,k+1,ii,5)/q(j,k+1,ii,1) -
     .       2.*(1.+bcj(k,ii,1))*qj0(k,ii,5,1)/qj0(k,ii,1,1) +
     .       2.*bcj(k,ii,1)*q(j,k,ii,5)/q(j,k,ii,1))

           t(jkv+k,4) = 0.25*(
     .       q(j,k+1,ii,2) -
     .       (1.+bck(j,ii,1))*qk0(j,ii,2,1) +
     .       bck(j,ii,1)*q(j,k,ii,2) +
     .       2.*(1.+bcj(k+1,ii,1))*qj0(k+1,ii,2,1) -
     .       2.*bcj(k+1,ii,1)*q(j,k+1,ii,2) -
     .       2.*(1.+bcj(k,ii,1))*qj0(k,ii,2,1) +
     .       2.*bcj(k,ii,1)*q(j,k,ii,2))

           t(jkv+k,5) = 0.25*(
     .       q(j,k+1,ii,3) -
     .       (1.+bck(j,ii,1))*qk0(j,ii,3,1) +
     .       bck(j,ii,1)*q(j,k,ii,3) +
     .       2.*(1.+bcj(k+1,ii,1))*qj0(k+1,ii,3,1) -
     .       2.*bcj(k+1,ii,1)*q(j,k+1,ii,3) -
     .       2.*(1.+bcj(k,ii,1))*qj0(k,ii,3,1) +
     .       2.*bcj(k,ii,1)*q(j,k,ii,3))

           t(jkv+k,6) = 0.25*(
     .       q(j,k+1,ii,4) -
     .       (1.+bck(j,ii,1))*qk0(j,ii,4,1) +
     .       bck(j,ii,1)*q(j,k,ii,4) +
     .       2.*(1.+bcj(k+1,ii,1))*qj0(k+1,ii,4,1) -
     .       2.*bcj(k+1,ii,1)*q(j,k+1,ii,4) -
     .       2.*(1.+bcj(k,ii,1))*qj0(k,ii,4,1) +
     .       2.*bcj(k,ii,1)*q(j,k,ii,4))

          else if (k .eq. kdim1) then
c
c replace (k+1) data by boundary data 
c
           t(jkv+k,3) = 0.25*gprgm1*(
     .       (1.+bck(j,ii,2))*qk0(j,ii,5,3)/qk0(j,ii,1,3) -
     .       bck(j,ii,2)*q(j,k,ii,5)/q(j,k,ii,1) -
     .       q(j,k-1,ii,5)/q(j,k-1,ii,1) +
     .       2.*(1.+bcj(k,ii,1))*qj0(k,ii,5,1)/qj0(k,ii,1,1) -
     .       2.*bcj(k,ii,1)*q(j,k,ii,5)/q(j,k,ii,1) -
     .       2.*(1.+bcj(k-1,ii,1))*qj0(k-1,ii,5,1)/qj0(k-1,ii,1,1) +
     .       2.*bcj(k-1,ii,1)*q(j,k-1,ii,5)/q(j,k-1,ii,1) )

           t(jkv+k,4) = 0.25*(
     .       (1.+bck(j,ii,2))*qk0(j,ii,2,3) -
     .       bck(j,ii,2)*q(j,k,ii,2) -
     .       q(j,k-1,ii,2) +
     .       2.*(1.+bcj(k,ii,1))*qj0(k,ii,2,1) -
     .       2.*bcj(k,ii,1)*q(j,k,ii,2) -
     .       2.*(1.+bcj(k-1,ii,1))*qj0(k-1,ii,2,1) +
     .       2.*bcj(k-1,ii,1)*q(j,k-1,ii,2) )
 
           t(jkv+k,5) = 0.25*(
     .       (1.+bck(j,ii,2))*qk0(j,ii,3,3) -
     .       bck(j,ii,2)*q(j,k,ii,3) -
     .       q(j,k-1,ii,3) +
     .       2.*(1.+bcj(k,ii,1))*qj0(k,ii,3,1) -
     .       2.*bcj(k,ii,1)*q(j,k,ii,3) -
     .       2.*(1.+bcj(k-1,ii,1))*qj0(k-1,ii,3,1) +
     .       2.*bcj(k-1,ii,1)*q(j,k-1,ii,3) )

           t(jkv+k,6) = 0.25*(
     .       (1.+bck(j,ii,2))*qk0(j,ii,4,3) -
     .       bck(j,ii,2)*q(j,k,ii,4) -
     .       q(j,k-1,ii,4) +
     .       2.*(1.+bcj(k,ii,1))*qj0(k,ii,4,1) -
     .       2.*bcj(k,ii,1)*q(j,k,ii,4) -
     .       2.*(1.+bcj(k-1,ii,1))*qj0(k-1,ii,4,1) +
     .       2.*bcj(k-1,ii,1)*q(j,k-1,ii,4) )
 
          else
c
c-- general case for j = 1
c
           t(jkv+k,3) = 0.25*gprgm1*(
     .       q(j,k+1,ii,5)/q(j,k+1,ii,1) -
     .       q(j,k-1,ii,5)/q(j,k-1,ii,1) +
     .       (1.+bcj(k+1,ii,1))*qj0(k+1,ii,5,1)/qj0(k+1,ii,1,1) -
     .       bcj(k+1,ii,1)*q(j,k+1,ii,5)/q(j,k+1,ii,1) -
     .       (1.+bcj(k-1,ii,1))*qj0(k-1,ii,5,1)/qj0(k-1,ii,1,1) +
     .       bcj(k-1,ii,1)*q(j,k-1,ii,5)/q(j,k-1,ii,1) )

           t(jkv+k,4) = 0.25*(
     .       q(j,k+1,ii,2) -
     .       q(j,k-1,ii,2) +
     .       (1.+bcj(k+1,ii,1))*qj0(k+1,ii,2,1) -
     .       bcj(k+1,ii,1)*q(j,k+1,ii,2) -
     .       (1.+bcj(k-1,ii,1))*qj0(k-1,ii,2,1) +
     .       bcj(k-1,ii,1)*q(j,k-1,ii,2) )

           t(jkv+k,5) = 0.25*(
     .       q(j,k+1,ii,3) -
     .       q(j,k-1,ii,3) +
     .       (1.+bcj(k+1,ii,1))*qj0(k+1,ii,3,1) -
     .       bcj(k+1,ii,1)*q(j,k+1,ii,3) -
     .       (1.+bcj(k-1,ii,1))*qj0(k-1,ii,3,1) +
     .       bcj(k-1,ii,1)*q(j,k-1,ii,3) )

           t(jkv+k,6) = 0.25*(
     .       q(j,k+1,ii,4) -
     .       q(j,k-1,ii,4) +
     .       (1.+bcj(k+1,ii,1))*qj0(k+1,ii,4,1) -
     .       bcj(k+1,ii,1)*q(j,k+1,ii,4) -
     .       (1.+bcj(k-1,ii,1))*qj0(k-1,ii,4,1) +
     .       bcj(k-1,ii,1)*q(j,k-1,ii,4) )
          end if
  225   continue
c
      else if (j.eq.jdim) then
c
c replace j data by boundary data
c for missing data at corners use one-sided differencing of known boundary values
c
        do 227 k = 1,kdim1
          if (k. eq. 1) then
c
c replace (k-1) data by boundary data 
c
           t(jkv+k,3) = 0.25*gprgm1*(
     .       2.*(1.+bcj(k+1,ii,2))*qj0(k+1,ii,5,3)/qj0(k+1,ii,1,3) -
     .       2.*bcj(k+1,ii,2)*q(j-1,k+1,ii,5)/q(j-1,k+1,ii,1) -
     .       2.*(1.+bcj(k,ii,2))*qj0(k,ii,5,3)/qj0(k,ii,1,3) +
     .       2.*bcj(k,ii,2)*q(j-1,k,ii,5)/q(j-1,k,ii,1) +
     .       q(j-1,k+1,ii,5)/q(j-1,k+1,ii,1) -
     .       (1.+bck(j-1,ii,1))*qk0(j-1,ii,5,1)/qk0(j-1,ii,1,1) +
     .       bck(j-1,ii,1)*q(j-1,k,ii,5)/q(j-1,k,ii,1) )

           t(jkv+k,4) = 0.25*(
     .       2.*(1.+bcj(k+1,ii,2))*qj0(k+1,ii,2,3) -
     .       2.*bcj(k+1,ii,2)*q(j-1,k+1,ii,2) -
     .       2.*(1.+bcj(k,ii,2))*qj0(k,ii,2,3) +
     .       2.*bcj(k,ii,2)*q(j-1,k,ii,2) +
     .       q(j-1,k+1,ii,2) -
     .       (1.+bck(j-1,ii,1))*qk0(j-1,ii,2,1) +
     .       bck(j-1,ii,1)*q(j-1,k,ii,2) )

           t(jkv+k,5) = 0.25*(
     .       2.*(1.+bcj(k+1,ii,2))*qj0(k+1,ii,3,3) -
     .       2.*bcj(k+1,ii,2)*q(j-1,k+1,ii,3) -
     .       2.*(1.+bcj(k,ii,2))*qj0(k,ii,3,3) +
     .       2.*bcj(k,ii,2)*q(j-1,k,ii,3) +
     .       q(j-1,k+1,ii,3) -
     .       (1.+bck(j-1,ii,1))*qk0(j-1,ii,3,1) +
     .       bck(j-1,ii,1)*q(j-1,k,ii,3) )
 
           t(jkv+k,6) = 0.25*(
     .       2.*(1.+bcj(k+1,ii,2))*qj0(k+1,ii,4,3) -
     .       2.*bcj(k+1,ii,2)*q(j-1,k+1,ii,4) -
     .       2.*(1.+bcj(k,ii,2))*qj0(k,ii,4,3) +
     .       2.*bcj(k,ii,2)*q(j-1,k,ii,4) +
     .       q(j-1,k+1,ii,4) -
     .       (1.+bck(j-1,ii,1))*qk0(j-1,ii,4,1) +
     .       bck(j-1,ii,1)*q(j-1,k,ii,4) )

          else if (k .eq. kdim1) then
c
c replace (k+1) data by boundary data 
c
           t(jkv+k,3) = 0.25*gprgm1*(
     .       2.*(1.+bcj(k,ii,2))*qj0(k,ii,5,3)/qj0(k,ii,1,3) -
     .       2.*bcj(k,ii,2)*q(j-1,k,ii,5)/q(j-1,k,ii,1) -
     .       2.*(1.+bcj(k-1,ii,2))*qj0(k-1,ii,5,3)/qj0(k-1,ii,1,3) +
     .       2.*bcj(k-1,ii,2)*q(j-1,k-1,ii,5)/q(j-1,k-1,ii,1) +
     .       (1.+bck(j-1,ii,2))*qk0(j-1,ii,5,3)/qk0(j-1,ii,1,3) -
     .       bck(j-1,ii,2)*q(j-1,k,ii,5)/q(j-1,k,ii,1) -
     .       q(j-1,k-1,ii,5)/q(j-1,k-1,ii,1) )

           t(jkv+k,4) = 0.25*(
     .       2.*(1.+bcj(k,ii,2))*qj0(k,ii,2,3) -
     .       2.*bcj(k,ii,2)*q(j-1,k,ii,2) -
     .       2.*(1.+bcj(k-1,ii,2))*qj0(k-1,ii,2,3) +
     .       2.*bcj(k-1,ii,2)*q(j-1,k-1,ii,2) +
     .       (1.+bck(j-1,ii,2))*qk0(j-1,ii,2,3) -
     .       bck(j-1,ii,2)*q(j-1,k,ii,2) -
     .       q(j-1,k-1,ii,2) )
 
           t(jkv+k,5) = 0.25*(
     .       2.*(1.+bcj(k,ii,2))*qj0(k,ii,3,3) -
     .       2.*bcj(k,ii,2)*q(j-1,k,ii,3) -
     .       2.*(1.+bcj(k-1,ii,2))*qj0(k-1,ii,3,3) +
     .       2.*bcj(k-1,ii,2)*q(j-1,k-1,ii,3) +
     .       (1.+bck(j-1,ii,2))*qk0(j-1,ii,3,3) -
     .       bck(j-1,ii,2)*q(j-1,k,ii,3) -
     .       q(j-1,k-1,ii,3) )
 
           t(jkv+k,6) = 0.25*(
     .       2.*(1.+bcj(k,ii,2))*qj0(k,ii,4,3) -
     .       2.*bcj(k,ii,2)*q(j-1,k,ii,4) -
     .       2.*(1.+bcj(k-1,ii,2))*qj0(k-1,ii,4,3) +
     .       2.*bcj(k-1,ii,2)*q(j-1,k-1,ii,4) +
     .       (1.+bck(j-1,ii,2))*qk0(j-1,ii,4,3) -
     .       bck(j-1,ii,2)*q(j-1,k,ii,4) -
     .       q(j-1,k-1,ii,4) )
 
          else
c
c-- general case for j = jdim
c
           t(jkv+k,3) = 0.25*gprgm1*(
     .       (1.+bcj(k+1,ii,2))*qj0(k+1,ii,5,3)/qj0(k+1,ii,1,3) -
     .       bcj(k+1,ii,2)*q(j-1,k+1,ii,5)/q(j-1,k+1,ii,1) -
     .       (1.+bcj(k-1,ii,2))*qj0(k-1,ii,5,3)/qj0(k-1,ii,1,3) +
     .       bcj(k-1,ii,2)*q(j-1,k-1,ii,5)/q(j-1,k-1,ii,1) +
     .       q(j-1,k+1,ii,5)/q(j-1,k+1,ii,1)-
     .       q(j-1,k-1,ii,5)/q(j-1,k-1,ii,1) )
 
           t(jkv+k,4) = 0.25*(
     .       (1.+bcj(k+1,ii,2))*qj0(k+1,ii,2,3) -
     .       bcj(k+1,ii,2)*q(j-1,k+1,ii,2) -
     .       (1.+bcj(k-1,ii,2))*qj0(k-1,ii,2,3) +
     .       bcj(k-1,ii,2)*q(j-1,k-1,ii,2) +
     .       q(j-1,k+1,ii,2)-
     .       q(j-1,k-1,ii,2) )

           t(jkv+k,5) = 0.25*(
     .       (1.+bcj(k+1,ii,2))*qj0(k+1,ii,3,3) -
     .       bcj(k+1,ii,2)*q(j-1,k+1,ii,3) -
     .       (1.+bcj(k-1,ii,2))*qj0(k-1,ii,3,3) +
     .       bcj(k-1,ii,2)*q(j-1,k-1,ii,3) +
     .       q(j-1,k+1,ii,3)-
     .       q(j-1,k-1,ii,3) )
 
           t(jkv+k,6) = 0.25*(
     .       (1.+bcj(k+1,ii,2))*qj0(k+1,ii,4,3) -
     .       bcj(k+1,ii,2)*q(j-1,k+1,ii,4) -
     .       (1.+bcj(k-1,ii,2))*qj0(k-1,ii,4,3) +
     .       bcj(k-1,ii,2)*q(j-1,k-1,ii,4) +
     .       q(j-1,k+1,ii,4)-
     .       q(j-1,k-1,ii,4) )
          end if
  227   continue
      else
c
c-- general case for j between 1 and jdim
c
        do 229 k = 1,kdim1
          if (k. eq. 1) then
c
c replace (k-1) data by boundary data 
c
           t(jkv+k,3) = 0.25*gprgm1*(
     .       q(j,k+1,ii,5)/q(j,k+1,ii,1) -
     .       (1.+bck(j,ii,1))*qk0(j,ii,5,1)/qk0(j,ii,1,1) +
     .       bck(j,ii,1)*q(j,k,ii,5)/q(j,k,ii,1) +
     .       q(j-1,k+1,ii,5)/ q(j-1,k+1,ii,1) -
     .       (1.+bck(j-1,ii,1))*qk0(j-1,ii,5,1)/qk0(j-1,ii,1,1) +
     .       bck(j-1,ii,1)*q(j-1,k,ii,5)/q(j-1,k,ii,1) )
 
           t(jkv+k,4) = 0.25*(
     .       q(j,k+1,ii,2) -
     .       (1.+bck(j,ii,1))*qk0(j,ii,2,1) +
     .       bck(j,ii,1)*q(j,k,ii,2) +
     .       q(j-1,k+1,ii,2) -
     .       (1.+bck(j-1,ii,1))*qk0(j-1,ii,2,1) +
     .       bck(j-1,ii,1)*q(j-1,k,ii,2) )

           t(jkv+k,5) = 0.25*(
     .       q(j,k+1,ii,3) -
     .       (1.+bck(j,ii,1))*qk0(j,ii,3,1) +
     .       bck(j,ii,1)*q(j,k,ii,3) +
     .       q(j-1,k+1,ii,3) -
     .       (1.+bck(j-1,ii,1))*qk0(j-1,ii,3,1) +
     .       bck(j-1,ii,1)*q(j-1,k,ii,3) )

           t(jkv+k,6) = 0.25*(
     .       q(j,k+1,ii,4) -
     .       (1.+bck(j,ii,1))*qk0(j,ii,4,1) +
     .       bck(j,ii,1)*q(j,k,ii,4) +
     .       q(j-1,k+1,ii,4) -
     .       (1.+bck(j-1,ii,1))*qk0(j-1,ii,4,1) +
     .       bck(j-1,ii,1)*q(j-1,k,ii,4) )

          else if (k .eq. kdim1) then
c
c replace (k+1) data by boundary data 
c
           t(jkv+k,3) = 0.25*gprgm1*(
     .       (1.+bck(j,ii,2))*qk0(j,ii,5,3)/qk0(j,ii,1,3) -
     .       bck(j,ii,2)*q(j,k,ii,5)/q(j,k,ii,1) - 
     .       q(j,k-1,ii,5)/ q(j,k-1,ii,1)  +
     .       (1.+bck(j-1,ii,2))*qk0(j-1,ii,5,3)/qk0(j-1,ii,1,3) -
     .       bck(j-1,ii,2)*q(j-1,k,ii,5)/q(j-1,k,ii,1) - 
     .       q(j-1,k-1,ii,5)/q(j-1,k-1,ii,1) )

           t(jkv+k,4) = 0.25*(
     .       (1.+bck(j,ii,2))*qk0(j,ii,2,3) -
     .       bck(j,ii,2)*q(j,k,ii,2) - 
     .       q(j,k-1,ii,2) +
     .       (1.+bck(j-1,ii,2))*qk0(j-1,ii,2,3) -
     .       bck(j-1,ii,2)*q(j-1,k,ii,2) - 
     .       q(j-1,k-1,ii,2) )

           t(jkv+k,5) = 0.25*(
     .       (1.+bck(j,ii,2))*qk0(j,ii,3,3) -
     .       bck(j,ii,2)*q(j,k,ii,3) - 
     .       q(j,k-1,ii,3)  +
     .       (1.+bck(j-1,ii,2))*qk0(j-1,ii,3,3) -
     .       bck(j-1,ii,2)*q(j-1,k,ii,3) - 
     .       q(j-1,k-1,ii,3) )

           t(jkv+k,6) = 0.25*(
     .       (1.+bck(j,ii,2))*qk0(j,ii,4,3) -
     .       bck(j,ii,2)*q(j,k,ii,4) - 
     .       q(j,k-1,ii,4)  +
     .       (1.+bck(j-1,ii,2))*qk0(j-1,ii,4,3) -
     .       bck(j-1,ii,2)*q(j-1,k,ii,4) - 
     .       q(j-1,k-1,ii,4) )
          else
c
c-- general case
c
            t(jkv+k,3) = 0.25*gprgm1*(
     .        q(j,k+1,ii,5)/q(j,k+1,ii,1) -
     .        q(j,k-1,ii,5)/q(j,k-1,ii,1) +
     .        q(j-1,k+1,ii,5)/q(j-1,k+1,ii,1) -
     .        q(j-1,k-1,ii,5)/q(j-1,k-1,ii,1) )

            t(jkv+k,4) = 0.25*(
     .        q(j,k+1,ii,2)   - q(j,k-1,ii,2) +
     .        q(j-1,k+1,ii,2) - q(j-1,k-1,ii,2) )

            t(jkv+k,5) = 0.25*(
     .        q(j,k+1,ii,3)   - q(j,k-1,ii,3) +
     .        q(j-1,k+1,ii,3) - q(j-1,k-1,ii,3) )

            t(jkv+k,6) = 0.25*(
     .        q(j,k+1,ii,4)   - q(j,k-1,ii,4) +
     .        q(j-1,k+1,ii,4) - q(j-1,k-1,ii,4) )
 
          end if
  229   continue
      end if
  230 continue
c******************************************************************************
c
c   gradients at cell interfaces (needed for full ns terms). Note that 
c  we don't need derivatives in eta direction, because they have already
c  been accounted for in gfluxv.
c  Recall bci=0 => cell data, bci=1 => face data
c
c  t(7)       : d( c*c/(gm1*pr) )/d(xie)
c  t(8)       : d(u)/d(xie)
c  t(9)       : d(v)/d(xie)
c  t(10)      : d(w)/d(xie)
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
        do 330 j = 1,jdim
          jkv     = (j-1)*npl*kdim1 + (ipl-1)*kdim1
          if (j .eq. 1) then
c
c replace (j-1) data by boundary data
c for missing data at corners use one-sided differencing of known boundary values
c
            do 323 k = 1,kdim1

            t(jkv+k,7) = 0.25*gprgm1*
     .        (q(j,k,ii+1,5)/q(j,k,ii+1,1)  -
     .        (1.+bci(j,k,1))*qi0(j,k,5,1)/qi0(j,k,1,1) +
     .        bci(j,k,1)*q(j,k,ii,5)/q(j,k,ii,1) +
     .        2.*(1.+bcj(k,ii+1,1))*qj0(k,ii+1,5,1)/qj0(k,ii+1,1,1) -
     .        2.*bcj(k,ii+1,1)*q(j,k,ii+1,5)/q(j,k,ii+1,1) -
     .        2.*(1.+bcj(k,ii,1))*qj0(k,ii,5,1)/qj0(k,ii,1,1) +
     .        2.*bcj(k,ii,1)*q(j,k,ii,5)/q(j,k,ii,1) )

            t(jkv+k,8) = 0.25*(
     .        q(j,k,ii+1,2)  -
     .        (1.+bci(j,k,1))*qi0(j,k,2,1) +
     .        bci(j,k,1)*q(j,k,ii,2) +
     .        2.*(1.+bcj(k,ii+1,1))*qj0(k,ii+1,2,1) -
     .        2.*bcj(k,ii+1,1)*q(j,k,ii+1,2) -
     .        2.*(1.+bcj(k,ii,1))*qj0(k,ii,2,1) +
     .        2.*bcj(k,ii,1)*q(j,k,ii,2))

            t(jkv+k,9) = 0.25*(
     .        q(j,k,ii+1,3)  -
     .        (1.+bci(j,k,1))*qi0(j,k,3,1) +
     .        bci(j,k,1)*q(j,k,ii,3) +
     .        2.*(1.+bcj(k,ii+1,1))*qj0(k,ii+1,3,1) -
     .        2.*bcj(k,ii+1,1)*q(j,k,ii+1,3) -
     .        2.*(1.+bcj(k,ii,1))*qj0(k,ii,3,1) +
     .        2.*bcj(k,ii,1)*q(j,k,ii,3))

            t(jkv+k,10)= 0.25*(
     .        q(j,k,ii+1,4)  -
     .        (1.+bci(j,k,1))*qi0(j,k,4,1) +
     .        bci(j,k,1)*q(j,k,ii,4) +
     .        2.*(1.+bcj(k,ii+1,1))*qj0(k,ii+1,4,1) -
     .        2.*bcj(k,ii+1,1)*q(j,k,ii+1,4) -
     .        2.*(1.+bcj(k,ii,1))*qj0(k,ii,4,1) +
     .        2.*bcj(k,ii,1)*q(j,k,ii,4))

  323       continue
          else if (j .eq. jdim) then
c
c replace j data by boundary data
c for missing data at corners use one-sided differencing of known boundary values
c
            do 325 k = 1,kdim1

            t(jkv+k,7) = 0.25*gprgm1*(
     .        2.*(1.+bcj(k,ii+1,2))*qj0(k,ii+1,5,3)/qj0(k,ii+1,1,3) -
     .        2.*bcj(k,ii+1,2)*q(j-1,k,ii+1,5)/q(j-1,k,ii+1,1) -
     .        2.*(1.+bcj(k,ii,2))*qj0(k,ii,5,3)/qj0(k,ii,1,3) +
     .        2.*bcj(k,ii,2)*q(j-1,k,ii,5)/q(j-1,k,ii,1) +
     .        q(j-1,k,ii+1,5)/q(j-1,k,ii+1,1) -
     .        (1.+bci(j-1,k,1))*qi0(j-1,k,5,1)/qi0(j-1,k,1,1) +
     .        bci(j-1,k,1)*q(j-1,k,ii,5)/q(j-1,k,ii,1) )

            t(jkv+k,8) = 0.25*(
     .        2.*(1.+bcj(k,ii+1,2))*qj0(k,ii+1,2,3) -
     .        2.*bcj(k,ii+1,2)*q(j-1,k,ii+1,2) -
     .        2.*(1.+bcj(k,ii,2))*qj0(k,ii,2,3) +
     .        2.*bcj(k,ii,2)*q(j-1,k,ii,2) +
     .        q(j-1,k,ii+1,2) -
     .        (1.+bci(j-1,k,1))*qi0(j-1,k,2,1) +
     .        bci(j-1,k,1)*q(j-1,k,ii,2) )

            t(jkv+k,9) = 0.25*(
     .        2.*(1.+bcj(k,ii+1,2))*qj0(k,ii+1,3,3) -
     .        2.*bcj(k,ii+1,2)*q(j-1,k,ii+1,3) -
     .        2.*(1.+bcj(k,ii,2))*qj0(k,ii,3,3) +
     .        2.*bcj(k,ii,2)*q(j-1,k,ii,3) +
     .        q(j-1,k,ii+1,3) -
     .        (1.+bci(j-1,k,1))*qi0(j-1,k,3,1) +
     .        bci(j-1,k,1)*q(j-1,k,ii,3) )

            t(jkv+k,10)= 0.25*(
     .        2.*(1.+bcj(k,ii+1,2))*qj0(k,ii+1,4,3) -
     .        2.*bcj(k,ii+1,2)*q(j-1,k,ii+1,4) -
     .        2.*(1.+bcj(k,ii,2))*qj0(k,ii,4,3) +
     .        2.*bcj(k,ii,2)*q(j-1,k,ii,4) +
     .        q(j-1,k,ii+1,4) -
     .        (1.+bci(j-1,k,1))*qi0(j-1,k,4,1) +
     .        bci(j-1,k,1)*q(j-1,k,ii,4) )

  325       continue
          else
c
c-- general case for ii = 1
c
            do 327 k = 1,kdim1

            t(jkv+k,7) = 0.25*gprgm1*(
     .        q(j,k,ii+1,5)/q(j,k,ii+1,1) -
     .        (1.+bci(j,k,1))*qi0(j,k,5,1)/qi0(j,k,1,1) +
     .        bci(j,k,1)*q(j,k,ii,5)/q(j,k,ii,1) +
     .        q(j-1,k,ii+1,5)/q(j-1,k,ii+1,1) -
     .        (1.+bci(j-1,k,1))*qi0(j-1,k,5,1)/qi0(j-1,k,1,1) +
     .        bci(j-1,k,1)*q(j-1,k,ii,5)/q(j-1,k,ii,1) )

            t(jkv+k,8) = 0.25*(
     .        q(j,k,ii+1,2)  -
     .        (1.+bci(j,k,1))*qi0(j,k,2,1) +
     .        bci(j,k,1)*q(j,k,ii,2) +
     .        q(j-1,k,ii+1,2) -
     .        (1.+bci(j-1,k,1))*qi0(j-1,k,2,1) +
     .        bci(j-1,k,1)*q(j-1,k,ii,2) )

            t(jkv+k,9) = 0.25*(
     .        q(j,k,ii+1,3)  -
     .        (1.+bci(j,k,1))*qi0(j,k,3,1) +
     .        bci(j,k,1)*q(j,k,ii,3) +
     .        q(j-1,k,ii+1,3) -
     .        (1.+bci(j-1,k,1))*qi0(j-1,k,3,1) +
     .        bci(j-1,k,1)*q(j-1,k,ii,3) )

            t(jkv+k,10)= 0.25*(
     .        q(j,k,ii+1,4)  -
     .        (1.+bci(j,k,1))*qi0(j,k,4,1) +
     .        bci(j,k,1)*q(j,k,ii,4) +
     .        q(j-1,k,ii+1,4) -
     .        (1.+bci(j-1,k,1))*qi0(j-1,k,4,1) +
     .        bci(j-1,k,1)*q(j-1,k,ii,4) )

  327       continue
          end if
  330   continue

      else if (ii .eq. idim1) then
c
c replace (ii+1) data by boundary data
c
        do 340 j = 1,jdim
          jkv     = (j-1)*npl*kdim1 + (ipl-1)*kdim1
          if (j .eq. 1) then
c
c replace (j-1) data by boundary data
c for missing data at corners use one-sided differencing of known boundary values
c
            do 333 k = 1,kdim1

            t(jkv+k,7) = 0.25*gprgm1*(
     .        (1.+bci(j,k,2))*qi0(j,k,5,3)/qi0(j,k,1,3) -
     .        bci(j,k,2)*q(j,k,ii,5)/q(j,k,ii,1) -
     .        q(j,k,ii-1,5)/q(j,k,ii-1,1) +
     .        2.*(1.+bcj(k,ii,1))*qj0(k,ii,5,1)/qj0(k,ii,1,1) -
     .        2.*bcj(k,ii,1)*q(j,k,ii,5)/q(j,k,ii,1) -
     .        2.*(1.+bcj(k,ii-1,1))*qj0(k,ii-1,5,1)/qj0(k,ii-1,1,1) +
     .        2.*bcj(k,ii-1,1)*q(j,k,ii-1,5)/q(j,k,ii-1,1) )

            t(jkv+k,8) = 0.25*(
     .        (1.+bci(j,k,2))*qi0(j,k,2,3) -
     .        bci(j,k,2)*q(j,k,ii,2) -
     .        q(j,k,ii-1,2) +
     .        2.*(1.+bcj(k,ii,1))*qj0(k,ii,2,1) -
     .        2.*bcj(k,ii,1)*q(j,k,ii,2) -
     .        2.*(1.+bcj(k,ii-1,1))*qj0(k,ii-1,2,1) +
     .        2.*bcj(k,ii-1,1)*q(j,k,ii-1,2) )

            t(jkv+k,9) = 0.25*(
     .        (1.+bci(j,k,2))*qi0(j,k,3,3) -
     .        bci(j,k,2)*q(j,k,ii,3) -
     .        q(j,k,ii-1,3) +
     .        2.*(1.+bcj(k,ii,1))*qj0(k,ii,3,1) -
     .        2.*bcj(k,ii,1)*q(j,k,ii,3) -
     .        2.*(1.+bcj(k,ii-1,1))*qj0(k,ii-1,3,1) +
     .        2.*bcj(k,ii-1,1)*q(j,k,ii-1,3) )

            t(jkv+k,10)= 0.25*(
     .        (1.+bci(j,k,2))*qi0(j,k,4,3) -
     .        bci(j,k,2)*q(j,k,ii,4) -
     .        q(j,k,ii-1,4) +
     .        2.*(1.+bcj(k,ii,1))*qj0(k,ii,4,1) -
     .        2.*bcj(k,ii,1)*q(j,k,ii,4) -
     .        2.*(1.+bcj(k,ii-1,1))*qj0(k,ii-1,4,1) +
     .        2.*bcj(k,ii-1,1)*q(j,k,ii-1,4) )

  333       continue
          else if (j .eq. jdim) then
c
c replace j data by boundary data
c for missing data at corners use one-sided differencing of known boundary values
c
            do 335 k = 1,kdim1

            t(jkv+k,7) = 0.25*gprgm1*(
     .        2.*(1.+bcj(k,ii,2))*qj0(k,ii,5,3)/qj0(k,ii,1,3) -
     .        2.*bcj(k,ii,2)*q(j-1,k,ii,5)/q(j-1,k,ii,1) -
     .        2.*(1.+bcj(k,ii-1,2))*qj0(k,ii-1,5,3)/qj0(k,ii-1,1,3) +
     .        2.*bcj(k,ii-1,2)*q(j-1,k,ii-1,5)/q(j-1,k,ii-1,1) +
     .        (1.+bci(j-1,k,2))*qi0(j-1,k,5,3)/qi0(j-1,k,1,3) -
     .        bci(j-1,k,2)*q(j-1,k,ii,5)/q(j-1,k,ii,1) -
     .        q(j-1,k,ii-1,5)/q(j-1,k,ii-1,1) )

            t(jkv+k,8) = 0.25*(
     .        2.*(1.+bcj(k,ii,2))*qj0(k,ii,2,3) -
     .        2.*bcj(k,ii,2)*q(j-1,k,ii,2) -
     .        2.*(1.+bcj(k,ii-1,2))*qj0(k,ii-1,2,3) +
     .        2.*bcj(k,ii-1,2)*q(j-1,k,ii-1,2) +
     .        (1.+bci(j-1,k,2))*qi0(j-1,k,2,3) -
     .        bci(j-1,k,2)*q(j-1,k,ii,2) -
     .        q(j-1,k,ii-1,4) )

            t(jkv+k,9) = 0.25*(
     .        2.*(1.+bcj(k,ii,2))*qj0(k,ii,3,3) -
     .        2.*bcj(k,ii,2)*q(j-1,k,ii,3) -
     .        2.*(1.+bcj(k,ii-1,2))*qj0(k,ii-1,3,3) +
     .        2.*bcj(k,ii-1,2)*q(j-1,k,ii-1,3) +
     .        (1.+bci(j-1,k,2))*qi0(j-1,k,3,3) -
     .        bci(j-1,k,2)*q(j-1,k,ii,3) -
     .        q(j-1,k,ii-1,3) )

            t(jkv+k,10) = 0.25*(
     .        2.*(1.+bcj(k,ii,2))*qj0(k,ii,4,3) -
     .        2.*bcj(k,ii,2)*q(j-1,k,ii,4) -
     .        2.*(1.+bcj(k,ii-1,2))*qj0(k,ii-1,4,3) +
     .        2.*bcj(k,ii-1,2)*q(j-1,k,ii-1,4) +
     .        (1.+bci(j-1,k,2))*qi0(j-1,k,4,3) -
     .        bci(j-1,k,2)*q(j-1,k,ii,4) -
     .        q(j-1,k,ii-1,4) )

  335       continue
          else
c
c-- general case for ii = idim1
c
            do 337 k = 1,kdim1

            t(jkv+k,7) = 0.25*gprgm1*(
     .        (1.+bci(j,k,2))*qi0(j,k,5,3)/qi0(j,k,1,3) -
     .        bci(j,k,2)*q(j,k,ii,5)/q(j,k,ii,1) -
     .        q(j,k,ii-1,5)/q(j,k,ii-1,1) +
     .        (1.+bci(j-1,k,2))*qi0(j-1,k,5,3)/qi0(j-1,k,1,3) -
     .        bci(j-1,k,2)*q(j-1,k,ii,5)/q(j-1,k,ii,1) -
     .        q(j-1,k,ii-1,5)/ q(j-1,k,ii-1,1) )

            t(jkv+k,8) = 0.25*(
     .        (1.+bci(j,k,2))*qi0(j,k,2,3) -
     .        bci(j,k,2)*q(j,k,ii,2) -
     .        q(j,k,ii-1,2) +
     .        (1.+bci(j-1,k,2))*qi0(j-1,k,2,3) -
     .        bci(j-1,k,2)*q(j-1,k,ii,2) -
     .        q(j-1,k,ii-1,2) )

            t(jkv+k,9) = 0.25*(
     .        (1.+bci(j,k,2))*qi0(j,k,3,3) -
     .        bci(j,k,2)*q(j,k,ii,3) -
     .        q(j,k,ii-1,3) +
     .        (1.+bci(j-1,k,2))*qi0(j-1,k,3,3) -
     .        bci(j-1,k,2)*q(j-1,k,ii,3) -
     .        q(j-1,k,ii-1,3) )

            t(jkv+k,10)= 0.25*(
     .        (1.+bci(j,k,2))*qi0(j,k,4,3) -
     .        bci(j,k,2)*q(j,k,ii,4) -
     .        q(j,k,ii-1,4) +
     .        (1.+bci(j-1,k,2))*qi0(j-1,k,4,3) -
     .        bci(j-1,k,2)*q(j-1,k,ii,4) -
     .        q(j-1,k,ii-1,4) )

  337       continue
          end if
  340   continue
      else
c
c-- general case (ii .ne 1. and ii .ne. idim1)
c

        do 350 j = 1,jdim
          jkv     = (j-1)*npl*kdim1 + (ipl-1)*kdim1
          if (j .eq. 1) then
c
c replace (j-1) data by boundary data
c
            do 343 k = 1,kdim1

            t(jkv+k,7) = 0.25*gprgm1*(
     .        q(j,k,ii+1,5)/q(j,k,ii+1,1) -
     .        q(j,k,ii-1,5)/q(j,k,ii-1,1) +
     .        (1.+bcj(k,ii+1,1))*qj0(k,ii+1,5,1)/qj0(k,ii+1,1,1) -
     .        bcj(k,ii+1,1)*q(j,k,ii+1,5)/q(j,k,ii+1,1) -
     .        (1.+bcj(k,ii-1,1))*qj0(k,ii-1,5,1)/qj0(k,ii-1,1,1) +
     .        bcj(k,ii-1,1)*q(j,k,ii-1,5)/q(j,k,ii-1,1) )

            t(jkv+k,8) = 0.25*(
     .        q(j,k,ii+1,2) -q(j,k,ii-1,2) +
     .        (1.+bcj(k,ii+1,1))*qj0(k,ii+1,2,1) -
     .        bcj(k,ii+1,1)*q(j,k,ii+1,2) -
     .        (1.+bcj(k,ii-1,1))*qj0(k,ii-1,2,1) +
     .        bcj(k,ii-1,1)*q(j,k,ii-1,2) )

            t(jkv+k,9) = 0.25*(
     .        q(j,k,ii+1,3) -q(j,k,ii-1,3) +
     .        (1.+bcj(k,ii+1,1))*qj0(k,ii+1,3,1) -
     .        bcj(k,ii+1,1)*q(j,k,ii+1,3) -
     .        (1.+bcj(k,ii-1,1))*qj0(k,ii-1,3,1) +
     .        bcj(k,ii-1,1)*q(j,k,ii-1,3) )

            t(jkv+k,10)= 0.25*(
     .        q(j,k,ii+1,4) -q(j,k,ii-1,4) +
     .        (1.+bcj(k,ii+1,1))*qj0(k,ii+1,4,1) -
     .        bcj(k,ii+1,1)*q(j,k,ii+1,4) -
     .        (1.+bcj(k,ii-1,1))*qj0(k,ii-1,4,1) +
     .        bcj(k,ii-1,1)*q(j,k,ii-1,4) )

  343       continue
          else if (j .eq. jdim) then
c
c replace j data by boundary data
c
            do 345 k = 1,kdim1

            t(jkv+k,7)= 0.25*gprgm1*(
     .        (1.+bcj(k,ii+1,2))*qj0(k,ii+1,5,3)/qj0(k,ii+1,1,3) -
     .        bcj(k,ii+1,2)*q(j-1,k,ii+1,5)/q(j-1,k,ii+1,1) -
     .        (1.+bcj(k,ii-1,2))*qj0(k,ii-1,5,3)/qj0(k,ii-1,1,3) +
     .        bcj(k,ii-1,2)*q(j-1,k,ii-1,5)/q(j-1,k,ii-1,1) +
     .        q(j-1,k,ii+1,5)/q(j-1,k,ii+1,1) -
     .        q(j-1,k,ii-1,5)/q(j-1,k,ii-1,1) )

            t(jkv+k,8) = 0.25*(
     .        (1.+bcj(k,ii+1,2))*qj0(k,ii+1,2,3) -
     .        bcj(k,ii+1,2)*q(j-1,k,ii+1,2) -
     .        (1.+bcj(k,ii-1,2))*qj0(k,ii-1,2,3) +
     .        bcj(k,ii-1,2)*q(j-1,k,ii-1,2) +
     .        q(j-1,k,ii+1,2)-q(j-1,k,ii-1,2) )

            t(jkv+k,9) = 0.25*(
     .        (1.+bcj(k,ii+1,2))*qj0(k,ii+1,3,3) -
     .        bcj(k,ii+1,2)*q(j-1,k,ii+1,3) -
     .        (1.+bcj(k,ii-1,2))*qj0(k,ii-1,3,3) +
     .        bcj(k,ii-1,2)*q(j-1,k,ii-1,3) +
     .        q(j-1,k,ii+1,3)-q(j-1,k,ii-1,3) )

            t(jkv+k,10)= 0.25*(
     .        (1.+bcj(k,ii+1,2))*qj0(k,ii+1,4,3) -
     .        bcj(k,ii+1,2)*q(j-1,k,ii+1,4) -
     .        (1.+bcj(k,ii-1,2))*qj0(k,ii-1,4,3) +
     .        bcj(k,ii-1,2)*q(j-1,k,ii-1,4) +
     .        q(j-1,k,ii+1,4)-q(j-1,k,ii-1,4) )

  345       continue
          else
            do 347 k = 1,kdim1
c
c-- general case
c
            t(jkv+k,7)= 0.25*gprgm1*(
     .        q(j,k,ii+1,5)/q(j,k,ii+1,1) -
     .        q(j,k,ii-1,5)/q(j,k,ii-1,1) +
     .        q(j-1,k,ii+1,5)/q(j-1,k,ii+1,1) -
     .        q(j-1,k,ii-1,5)/q(j-1,k,ii-1,1) )

            t(jkv+k,8) = 0.25*(
     .        q(j,k,ii+1,2)  -q(j,k,ii-1,2) +
     .        q(j-1,k,ii+1,2)-q(j-1,k,ii-1,2) )

            t(jkv+k,9) = 0.25*(
     .        q(j,k,ii+1,3)  -q(j,k,ii-1,3) +
     .        q(j-1,k,ii+1,3)-q(j-1,k,ii-1,3) )

            t(jkv+k,10)= 0.25*(
     .        q(j,k,ii+1,4)  -q(j,k,ii-1,4) +
     .        q(j-1,k,ii+1,4)-q(j-1,k,ii-1,4) )

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
      t(izz+ks-1,11) =0.5e0*(t(izz+ks-1,18)+t(izz,18))
      t(izz+ks-1,12) =0.5e0*(t(izz+ks-1,19)+t(izz,19))
      t(izz+ks-1,13) =0.5e0*(t(izz+ks-1,20)+t(izz,20))
 1024 continue
c     interfaces at j=0 and j=jdim
cdir$ ivdep
      do 1025 izz=1,kv
      ab=1.+wj0(izz,15)
      bb=1.-wj0(izz,15)
      t(izz+n,11)=0.5*(ab*wj0(izz,18)+bb*t(izz+n-kv,18))
      t(izz+n,12)=0.5*(ab*wj0(izz,19)+bb*t(izz+n-kv,19))
      t(izz+n,13)=0.5*(ab*wj0(izz,20)+bb*t(izz+n-kv,20))
      ab=1.+wj0(izz,5)
      bb=1.-wj0(izz,5)
      t(izz,11) = 0.5*(ab*wj0(izz,8) +bb*t(izz,18))
      t(izz,12) = 0.5*(ab*wj0(izz,9) +bb*t(izz,19))
      t(izz,13) = 0.5*(ab*wj0(izz,10) +bb*t(izz,20))
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
c-- form ux, uy, uz, vx, vy, vz, wx, wy, wz, tx, ty and tz
c   excluding contributions in j/eta direction (already accounted for in 
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
c--  Note lambda = -2/3 mu, lambda + 2 mu = 4/3 mu 
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
c  (-)viscous flux =  Gv(j-1/2)-Gv(j+1/2)
cdir$ ivdep
      do 1027 izz=1,n
      t(izz,2)       = -t(izz+ks-1,16)+t(izz,16)
      t(izz,3)       = -t(izz+ks-1,17)+t(izz,17)
      t(izz,4)       = -t(izz+ks-1,18)+t(izz,18)
      t(izz,5)       = -t(izz+ks-1,19)+t(izz,19)
 1027 continue
c
c******************************************************************************
c
c  calculate residuals
c
      do 400 ipl=1,npl
      ii            = i+ipl-1
      do 400 j=1,jdim1
      jkv           = (j-1)*kdim1*npl + (ipl-1)*kdim1
c   for 2-D, t(3) should be identically zero
      if (i2d .eq. 1) then
        do k=1,kdim1
          t(k+jkv,3)=0.
        enddo
      end if
      do 400 l=2,5
cdir$ ivdep
      do 1030 k=1,kdim1
      res(j,k,ii,l) = res(j,k,ii,l)+t(k+jkv,l)
 1030 continue
  400 continue
      return
      end 
