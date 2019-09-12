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
      subroutine ffluxv1(k,npl,jdim,kdim,idim,res,q,qj0,qk0,qi0,
     .                  sj,sk,si,vol,t,nvtq,wi0,vist3d,vi0,
     .                  bcj,bck,bci,voli0,
     .                  nou,bou,nbuf,ibufdim,iadv)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Calculate right-hand residual contributions in the
c     I-direction due to the viscous terms missing in the thin-layer
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
      dimension si(jdim,kdim,idim,5),vol(jdim,kdim,idim-1)
      dimension sk(jdim,kdim,idim-1,5),sj(jdim,kdim,idim-1,5)

      dimension q(jdim,kdim,idim,5)
      dimension qj0(kdim,idim-1,5,4), qk0(jdim,idim-1,5,4),
     .          qi0(jdim,kdim,5,4)
      dimension res(jdim,kdim,idim-1,5),vist3d(jdim,kdim,idim),
     .          vi0(jdim,kdim,1,4),voli0(jdim,kdim,4)
      dimension t(nvtq,32),wi0(npl*(jdim-1),22)
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
      common /constit/ i_nonlin,c_nonlin,snonlin_lim,i_tauijs,i_qcr2000,
     .                 i_qcr2013,i_qcr2013v
      REAL :: coef_eddy  ! used to turn off eddy viscosity for stress calculation

      coef_eddy = 1.0
      if(ivisc(1)>=70) coef_eddy=0.
      if(i_tauijs .eq. 1) coef_eddy=0.
c
c******************************************************************************
c     written mid-2005   V. Venkatakrishnan & C. Rumsey
c******************************************************************************
c
      idim1  =  idim-1
      jdim1  =  jdim-1
      kdim1  =  kdim-1
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
      gprgm1 = gpr/gm1
c
      if (isklton.gt.0 .and. k.eq.1 .and. iadv.ge.0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1444)
      end if
 1444 format(52h   computing cross-derivative viscous fluxes, I-dir.)
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
      t1          =  si(j,kk,i,4)*t(j+jiv,15) 
      t(j+jiv,25) =  si(j,kk,i,1)*t1
      t(j+jiv,26) =  si(j,kk,i,2)*t1
      t(j+jiv,27) =  si(j,kk,i,3)*t1
 1014 continue
  200 continue
c
c t(28)  - t(30) : Components of grad (eta)
c grad(eta) = 1/4 [(sj(j,k,i,l)*sj(j,k,i,4) + sj(j,k,i-1,l)*sj(j,k,i-1,4) +
c sj(j+1,k,i,l)*sj(j+1,k,i,4) + sj(j+1,k,i-1,l)*sj(j+1,k,i-1,4)] /V_int,
c i=1,idim , j=1,jdim1, k=1,npl, l =1,3.  At the boundaries in i,
c since we don't have sj's, we assume that them to be equal to the interior one.
c
c Note: the code for i=1 is obtained from the general case by substituting
c i for i-1, and likewise the code for i=idim by substituting i-1 for i.
c
      do 210 kpl=1,npl
      kk  = k+kpl-1
      do 210 i=1,idim
      jiv  = (i-1)*npl*jdim1 + (kpl-1)*jdim1 
      if (i.eq.1) then
       do 204 j = 1,jdim1
         t(jiv+j,28) = 0.25*t(j+jiv,15)*(sj(j,kk,i,1)*sj(j,kk,i,4)
     .              + sj(j,kk,i,1)*sj(j,kk,i,4)
     .              + sj(j+1,kk,i,1)*sj(j+1,kk,i,4)
     .              + sj(j+1,kk,i,1)*sj(j+1,kk,i,4))
         t(jiv+j,29) = 0.25*t(j+jiv,15)*(sj(j,kk,i,2)*sj(j,kk,i,4)
     .              + sj(j,kk,i,2)*sj(j,kk,i,4)
     .              + sj(j+1,kk,i,2)*sj(j+1,kk,i,4)
     .              + sj(j+1,kk,i,2)*sj(j+1,kk,i,4))
         t(jiv+j,30) = 0.25*t(j+jiv,15)*(sj(j,kk,i,3)*sj(j,kk,i,4)
     .              + sj(j,kk,i,3)*sj(j,kk,i,4)
     .              + sj(j+1,kk,i,3)*sj(j+1,kk,i,4)
     .              + sj(j+1,kk,i,3)*sj(j+1,kk,i,4))
  204  continue
      else if (i.eq.idim) then
       do 206 j = 1,jdim1
         t(jiv+j,28) = 0.25*t(j+jiv,15)*(sj(j,kk,i-1,1)*sj(j,kk,i-1,4)
     .              + sj(j,kk,i-1,1)*sj(j,kk,i-1,4)
     .              + sj(j+1,kk,i-1,1)*sj(j+1,kk,i-1,4)
     .              + sj(j+1,kk,i-1,1)*sj(j+1,kk,i-1,4))
         t(jiv+j,29) = 0.25*t(j+jiv,15)*(sj(j,kk,i-1,2)*sj(j,kk,i-1,4)
     .              + sj(j,kk,i-1,2)*sj(j,kk,i-1,4)
     .              + sj(j+1,kk,i-1,2)*sj(j+1,kk,i-1,4)
     .              + sj(j+1,kk,i-1,2)*sj(j+1,kk,i-1,4))
         t(jiv+j,30) = 0.25*t(j+jiv,15)*(sj(j,kk,i-1,3)*sj(j,kk,i-1,4)
     .              + sj(j,kk,i-1,3)*sj(j,kk,i-1,4)
     .              + sj(j+1,kk,i-1,3)*sj(j+1,kk,i-1,4)
     .              + sj(j+1,kk,i-1,3)*sj(j+1,kk,i-1,4))
  206  continue
      else
c
c-- general case
c
       do 208 j = 1,jdim1
         t(jiv+j,28) = 0.25*t(j+jiv,15)*(sj(j,kk,i,1)*sj(j,kk,i,4)
     .              + sj(j,kk,i-1,1)*sj(j,kk,i-1,4)
     .              + sj(j+1,kk,i,1)*sj(j+1,kk,i,4)
     .              + sj(j+1,kk,i-1,1)*sj(j+1,kk,i-1,4))
         t(jiv+j,29) = 0.25*t(j+jiv,15)*(sj(j,kk,i,2)*sj(j,kk,i,4)
     .              + sj(j,kk,i-1,2)*sj(j,kk,i-1,4)
     .              + sj(j+1,kk,i,2)*sj(j+1,kk,i,4)
     .              + sj(j+1,kk,i-1,2)*sj(j+1,kk,i-1,4))
         t(jiv+j,30) = 0.25*t(j+jiv,15)*(sj(j,kk,i,3)*sj(j,kk,i,4)
     .              + sj(j,kk,i-1,3)*sj(j,kk,i-1,4)
     .              + sj(j+1,kk,i,3)*sj(j+1,kk,i,4)
     .              + sj(j+1,kk,i-1,3)*sj(j+1,kk,i-1,4))
  208  continue
      end if
  210 continue
c
c t(31),t(32),t(2)  : Components of grad (zeta)
c grad(zeta) = 1/4 [(sk(j,k,i,l)*sk(j,k,i,4) + sk(j,k,i-1,l)*sk(j,k,i-1,4) +
c sk(j,k+1,i,l)*sk(j,k+1,i,4) + sk(j,k+1,i-1,l)*sk(j,k+1,i-1,4)] /V_int,
c i=1,idim , j=1,jdim1, k=1,npl,l = 1,3.  At the boundaries in i,
c since we don't have sk's, we assume that them to be equal to the interior one.
c
c Note: the code for i=1 is obtained from the general case by substituting
c i for i-1, and likewise the code for i=idim by substituting i-1 for i.
c
      do 220 kpl=1,npl
      kk  = k+kpl-1
      do 220 i=1,idim
      jiv  = (i-1)*npl*jdim1 + (kpl-1)*jdim1 
      if (i.eq.1) then
       do 214 j = 1,jdim1
        t(jiv+j,31) = 0.25*t(j+jiv,15)*(sk(j,kk,i,1)*sk(j,kk,i,4)
     .              + sk(j,kk,i,1)*sk(j,kk,i,4)
     .              + sk(j,kk+1,i,1)*sk(j,kk+1,i,4)
     .              + sk(j,kk+1,i,1)*sk(j,kk+1,i,4))
        t(jiv+j,32) = 0.25*t(j+jiv,15)*(sk(j,kk,i,2)*sk(j,kk,i,4)
     .              + sk(j,kk,i,2)*sk(j,kk,i,4)
     .              + sk(j,kk+1,i,2)*sk(j,kk+1,i,4)
     .              + sk(j,kk+1,i,2)*sk(j,kk+1,i,4))
        t(jiv+j,2)  = 0.25*t(j+jiv,15)*(sk(j,kk,i,3)*sk(j,kk,i,4)
     .              + sk(j,kk,i,3)*sk(j,kk,i,4)
     .              + sk(j,kk+1,i,3)*sk(j,kk+1,i,4)
     .              + sk(j,kk+1,i,3)*sk(j,kk+1,i,4))
  214  continue
      else if (i.eq.idim) then
       do 216 j = 1,jdim1
        t(jiv+j,31) = 0.25*t(j+jiv,15)*(sk(j,kk,i-1,1)*sk(j,kk,i-1,4)
     .              + sk(j,kk,i-1,1)*sk(j,kk,i-1,4)
     .              + sk(j,kk+1,i-1,1)*sk(j,kk+1,i-1,4)
     .              + sk(j,kk+1,i-1,1)*sk(j,kk+1,i-1,4))
        t(jiv+j,32) = 0.25*t(j+jiv,15)*(sk(j,kk,i-1,2)*sk(j,kk,i-1,4)
     .              + sk(j,kk,i-1,2)*sk(j,kk,i-1,4)
     .              + sk(j,kk+1,i-1,2)*sk(j,kk+1,i-1,4)
     .              + sk(j,kk+1,i-1,2)*sk(j,kk+1,i-1,4))
        t(jiv+j,2)  = 0.25*t(j+jiv,15)*(sk(j,kk,i-1,3)*sk(j,kk,i-1,4)
     .              + sk(j,kk,i-1,3)*sk(j,kk,i-1,4)
     .              + sk(j,kk+1,i-1,3)*sk(j,kk+1,i-1,4)
     .              + sk(j,kk+1,i-1,3)*sk(j,kk+1,i-1,4))
  216  continue
      else
c
c-- general case
c
       do 218 j = 1,jdim1
        t(jiv+j,31) = 0.25*t(j+jiv,15)*(sk(j,kk,i,1)*sk(j,kk,i,4)
     .              + sk(j,kk,i-1,1)*sk(j,kk,i-1,4)
     .              + sk(j,kk+1,i,1)*sk(j,kk+1,i,4)
     .              + sk(j,kk+1,i-1,1)*sk(j,kk+1,i-1,4))
        t(jiv+j,32) = 0.25*t(j+jiv,15)*(sk(j,kk,i,2)*sk(j,kk,i,4)
     .              + sk(j,kk,i-1,2)*sk(j,kk,i-1,4)
     .              + sk(j,kk+1,i,2)*sk(j,kk+1,i,4)
     .              + sk(j,kk+1,i-1,2)*sk(j,kk+1,i-1,4))
        t(jiv+j,2)  = 0.25*t(j+jiv,15)*(sk(j,kk,i,3)*sk(j,kk,i,4)
     .              + sk(j,kk,i-1,3)*sk(j,kk,i-1,4)
     .              + sk(j,kk+1,i,3)*sk(j,kk+1,i,4)
     .              + sk(j,kk+1,i-1,3)*sk(j,kk+1,i-1,4))
  218  continue
      end if
  220 continue
c
c
c******************************************************************************
c
c   gradients at cell interfaces (needed for full ns terms). Note that 
c  we don't need derivatives in xie direction, because they have already
c  been accounted for in ffluxv.
c  Recall bcj=0 => cell data, bcj=1 => face data
c
c  t(3)       : d( c*c/(gm1*pr) )/d(eta)
c  t(4)       : d(u)/d(eta)
c  t(5)       : d(v)/d(eta)
c  t(6)       : d(w)/d(eta)
c
c     interior interfaces
c
      do 230 kpl = 1,npl
      kk = k + kpl -1
      do 230 i = 1,idim
      jiv  = (i-1)*npl*jdim1 + (kpl-1)*jdim1 
c
      if (i.eq.1) then
c
c replace (i-1) data by boundary data
c for missing data at corners use one-sided differencing of known boundary values
c
        do 225 j = 1,jdim1
          if (j. eq. 1) then
c
c replace (j-1) data by boundary data
c
           t(jiv+j,3) = 0.25*gprgm1*(
     .       q(j+1,kk,i,5)/q(j+1,kk,i,1) -
     .       (1.+bcj(kk,i,1))*qj0(kk,i,5,1)/qj0(kk,i,1,1) +
     .       bcj(kk,i,1)*q(j,kk,i,5)/q(j,kk,i,1) +
     .       2.*(1.+bci(j+1,kk,1))*qi0(j+1,kk,5,1)/qi0(j+1,kk,1,1) -
     .       2.*bci(j+1,kk,1)*q(j+1,kk,i,5)/q(j+1,kk,i,1) -
     .       2.*(1.+bci(j,kk,1))*qi0(j,kk,5,1)/qi0(j,kk,1,1) +
     .       2.*bci(j,kk,1)*q(j,kk,i,5)/q(j,kk,i,1) )

           t(jiv+j,4) = 0.25*(
     .       q(j+1,kk,i,2) -
     .       (1.+bcj(kk,i,1))*qj0(kk,i,2,1) +
     .       bcj(kk,i,1)*q(j,kk,i,2) +
     .       2.*(1.+bci(j+1,kk,1))*qi0(j+1,kk,2,1) -
     .       2.*bci(j+1,kk,1)*q(j+1,kk,i,2) -
     .       2.*(1.+bci(j,kk,1))*qi0(j,kk,2,1) +
     .       2.*bci(j,kk,1)*q(j,kk,i,2) )

           t(jiv+j,5) = 0.25*(
     .       q(j+1,kk,i,3) -
     .       (1.+bcj(kk,i,1))*qj0(kk,i,3,1) +
     .       bcj(kk,i,1)*q(j,kk,i,3) +
     .       2.*(1.+bci(j+1,kk,1))*qi0(j+1,kk,3,1) -
     .       2.*bci(j+1,kk,1)*q(j+1,kk,i,3) -
     .       2.*(1.+bci(j,kk,1))*qi0(j,kk,3,1) +
     .       2.*bci(j,kk,1)*q(j,kk,i,3) )

           t(jiv+j,6) = 0.25*(
     .       q(j+1,kk,i,4) -
     .       (1.+bcj(kk,i,1))*qj0(kk,i,4,1) +
     .       bcj(kk,i,1)*q(j,kk,i,4) +
     .       2.*(1.+bci(j+1,kk,1))*qi0(j+1,kk,4,1) -
     .       2.*bci(j+1,kk,1)*q(j+1,kk,i,4) -
     .       2.*(1.+bci(j,kk,1))*qi0(j,kk,4,1) +
     .       2.*bci(j,kk,1)*q(j,kk,i,4) )

          else if (j .eq. jdim1) then
c
c replace (j+1) data by boundary data 
c
           t(jiv+j,3) = 0.25*gprgm1*(
     .       (1.+bcj(kk,i,2))*qj0(kk,i,5,3)/qj0(kk,i,1,3) -
     .       bcj(kk,i,2)*q(j,kk,i,5)/q(j,kk,i,1) - 
     .       q(j-1,kk,i,5)/q(j-1,kk,i,1)  +
     .       2.*(1.+bci(j,kk,1))*qi0(j,kk,5,1)/qi0(j,kk,1,1) -
     .       2.*bci(j,kk,1)*q(j,kk,i,5)/q(j,kk,i,1) -
     .       2.*(1.+bci(j-1,kk,1))*qi0(j-1,kk,5,1)/qi0(j-1,kk,1,1) +
     .       2.*bci(j-1,kk,1)*q(j-1,kk,i,5)/q(j-1,kk,i,1) )

           t(jiv+j,4) = 0.25*(
     .       (1.+bcj(kk,i,2))*qj0(kk,i,2,3) -
     .       bcj(kk,i,2)*q(j,kk,i,2) - 
     .       q(j-1,kk,i,2)  +
     .       2.*(1.+bci(j,kk,1))*qi0(j,kk,2,1) -
     .       2.*bci(j,kk,1)*q(j,kk,i,2) -
     .       2.*(1.+bci(j-1,kk,1))*qi0(j-1,kk,2,1) +
     .       2.*bci(j-1,kk,1)*q(j-1,kk,i,2) )

           t(jiv+j,5) = 0.25*(
     .       (1.+bcj(kk,i,2))*qj0(kk,i,3,3) -
     .       bcj(kk,i,2)*q(j,kk,i,3) - 
     .       q(j-1,kk,i,3)  +
     .       2.*(1.+bci(j,kk,1))*qi0(j,kk,3,1) -
     .       2.*bci(j,kk,1)*q(j,kk,i,3) -
     .       2.*(1.+bci(j-1,kk,1))*qi0(j-1,kk,3,1) +
     .       2.*bci(j-1,kk,1)*q(j-1,kk,i,3) )

           t(jiv+j,6) = 0.25*(
     .       (1.+bcj(kk,i,2))*qj0(kk,i,4,3) -
     .       bcj(kk,i,2)*q(j,kk,i,4) - 
     .       q(j-1,kk,i,4)  +
     .       2.*(1.+bci(j,kk,1))*qi0(j,kk,4,1) -
     .       2.*bci(j,kk,1)*q(j,kk,i,4) -
     .       2.*(1.+bci(j-1,kk,1))*qi0(j-1,kk,4,1) +
     .       2.*bci(j-1,kk,1)*q(j-1,kk,i,4) )

          else
c
c-- general case for i = 1
c
            t(jiv+j,3) = 0.25*gprgm1*(
     .       q(j+1,kk,i,5)/q(j+1,kk,i,1) -
     .       q(j-1,kk,i,5)/q(j-1,kk,i,1) +
     .       (1.+bci(j+1,kk,1))*qi0(j+1,kk,5,1)/qi0(j+1,kk,1,1) -
     .       bci(j+1,kk,1)*q(j+1,kk,i,5)/q(j+1,kk,i,1) -
     .       (1.+bci(j-1,kk,1))*qi0(j-1,kk,5,1)/qi0(j-1,kk,1,1) +
     .       bci(j-1,kk,1)*q(j-1,kk,i,5)/q(j-1,kk,i,1) )

            t(jiv+j,4) = 0.25*(
     .       q(j+1,kk,i,2) -
     .       q(j-1,kk,i,2) +
     .       (1.+bci(j+1,kk,1))*qi0(j+1,kk,2,1) -
     .       bci(j+1,kk,1)*q(j+1,kk,i,2) -
     .       (1.+bci(j-1,kk,1))*qi0(j-1,kk,2,1) +
     .       bci(j-1,kk,1)*q(j-1,kk,i,2) )

            t(jiv+j,5) = 0.25*(
     .       q(j+1,kk,i,3) -
     .       q(j-1,kk,i,3) +
     .       (1.+bci(j+1,kk,1))*qi0(j+1,kk,3,1) -
     .       bci(j+1,kk,1)*q(j+1,kk,i,3) -
     .       (1.+bci(j-1,kk,1))*qi0(j-1,kk,3,1) +
     .       bci(j-1,kk,1)*q(j-1,kk,i,3) )

            t(jiv+j,6) = 0.25*(
     .       q(j+1,kk,i,4) -
     .       q(j-1,kk,i,4) +
     .       (1.+bci(j+1,kk,1))*qi0(j+1,kk,4,1) -
     .       bci(j+1,kk,1)*q(j+1,kk,i,4) -
     .       (1.+bci(j-1,kk,1))*qi0(j-1,kk,4,1) +
     .       bci(j-1,kk,1)*q(j-1,kk,i,4) )

          end if
  225   continue
c
      else if (i.eq.idim) then
c
c replace i data by boundary data
c for missing data at corners use one-sided differencing of known boundary values
c
        do 227 j = 1,jdim1
          if (j. eq. 1) then
c
c replace (j-1) data by boundary data
c
           t(jiv+j,3) = 0.25*gprgm1*(
     .       2.*(1.+bci(j+1,kk,2))*qi0(j+1,kk,5,3)/qi0(j+1,kk,1,3) -
     .       2.*bci(j+1,kk,2)*q(j+1,kk,i-1,5)/q(j+1,kk,i-1,1) -
     .       2.*(1.+bci(j,kk,2))*qi0(j,kk,5,3)/qi0(j,kk,1,3) +
     .       2.*bci(j,kk,2)*q(j,kk,i-1,5)/q(j,kk,i-1,1) +
     .       q(j+1,kk,i-1,5)/q(j+1,kk,i-1,1) -
     .       (1.+bcj(kk,i-1,1))*qj0(kk,i-1,5,1)/qj0(kk,i-1,1,1) +
     .       bcj(kk,i-1,1)*q(j,kk,i-1,5)/q(j,kk,i-1,1) )

           t(jiv+j,4) = 0.25*(
     .       2.*(1.+bci(j+1,kk,2))*qi0(j+1,kk,2,3) -
     .       2.*bci(j+1,kk,2)*q(j+1,kk,i-1,2) -
     .       2.*(1.+bci(j,kk,2))*qi0(j,kk,2,3) +
     .       2.*bci(j,kk,2)*q(j,kk,i-1,2) +
     .       q(j+1,kk,i-1,2) -
     .       (1.+bcj(kk,i-1,1))*qj0(kk,i-1,2,1) +
     .       bcj(kk,i-1,1)*q(j,kk,i-1,2) )

           t(jiv+j,5) = 0.25*(
     .       2.*(1.+bci(j+1,kk,2))*qi0(j+1,kk,3,3) -
     .       2.*bci(j+1,kk,2)*q(j+1,kk,i-1,3) -
     .       2.*(1.+bci(j,kk,2))*qi0(j,kk,3,3) +
     .       2.*bci(j,kk,2)*q(j,kk,i-1,3) +
     .       q(j+1,kk,i-1,3) -
     .       (1.+bcj(kk,i-1,1))*qj0(kk,i-1,3,1) +
     .       bcj(kk,i-1,1)*q(j,kk,i-1,3) )

           t(jiv+j,6) = 0.25*(
     .       2.*(1.+bci(j+1,kk,2))*qi0(j+1,kk,4,3) -
     .       2.*bci(j+1,kk,2)*q(j+1,kk,i-1,4) -
     .       2.*(1.+bci(j,kk,2))*qi0(j,kk,4,3) +
     .       2.*bci(j,kk,2)*q(j,kk,i-1,4) +
     .       q(j+1,kk,i-1,4) -
     .       (1.+bcj(kk,i-1,1))*qj0(kk,i-1,4,1) +
     .       bcj(kk,i-1,1)*q(j,kk,i-1,4) )

          else if (j .eq. jdim1) then
c
c replace (j+1) data by boundary data
c
           t(jiv+j,3) = 0.25*gprgm1*(
     .       2.*(1.+bci(j,kk,2))*qi0(j,kk,5,3)/qi0(j,kk,1,3) -
     .       2.*bci(j,kk,2)*q(j,kk,i-1,5)/q(j,kk,i-1,1) -
     .       2.*(1.+bci(j-1,kk,2))*qi0(j-1,kk,5,3)/qi0(j-1,kk,1,3) +
     .       2.*bci(j-1,kk,2)*q(j-1,kk,i-1,5)/q(j-1,kk,i-1,1) +
     .       (1.+bcj(kk,i-1,2))*qj0(kk,i-1,5,3)/qj0(kk,i-1,1,3) -
     .       bcj(kk,i-1,2)*q(j,kk,i-1,5)/q(j,kk,i-1,1) - 
     .       q(j-1,kk,i-1,5)/q(j-1,kk,i-1,1) )

           t(jiv+j,4) = 0.25*(
     .       2.*(1.+bci(j,kk,2))*qi0(j,kk,2,3) -
     .       2.*bci(j,kk,2)*q(j,kk,i-1,2) -
     .       2.*(1.+bci(j-1,kk,2))*qi0(j-1,kk,2,3) +
     .       2.*bci(j-1,kk,2)*q(j-1,kk,i-1,2) +
     .       (1.+bcj(kk,i-1,2))*qj0(kk,i-1,2,3) -
     .       bcj(kk,i-1,2)*q(j,kk,i-1,2) - 
     .       q(j-1,kk,i-1,2) )

           t(jiv+j,5) = 0.25*(
     .       2.*(1.+bci(j,kk,2))*qi0(j,kk,3,3) -
     .       2.*bci(j,kk,2)*q(j,kk,i-1,3) -
     .       2.*(1.+bci(j-1,kk,2))*qi0(j-1,kk,3,3) +
     .       2.*bci(j-1,kk,2)*q(j-1,kk,i-1,3) +
     .       (1.+bcj(kk,i-1,2))*qj0(kk,i-1,3,3) -
     .       bcj(kk,i-1,2)*q(j,kk,i-1,3) - 
     .       q(j-1,kk,i-1,3) )

           t(jiv+j,6) = 0.25*(
     .       2.*(1.+bci(j,kk,2))*qi0(j,kk,4,3) -
     .       2.*bci(j,kk,2)*q(j,kk,i-1,4) -
     .       2.*(1.+bci(j-1,kk,2))*qi0(j-1,kk,4,3) +
     .       2.*bci(j-1,kk,2)*q(j-1,kk,i-1,4) +
     .       (1.+bcj(kk,i-1,2))*qj0(kk,i-1,4,3) -
     .       bcj(kk,i-1,2)*q(j,kk,i-1,4) - 
     .       q(j-1,kk,i-1,4) )
          else
c
c-- general case for i = idim
c
           t(jiv+j,3) = 0.25*gprgm1*(
     .       (1.+bci(j+1,kk,2))*qi0(j+1,kk,5,3)/qi0(j+1,kk,1,3) -
     .       bci(j+1,kk,2)*q(j+1,kk,i-1,5)/q(j+1,kk,i-1,1) -
     .       (1.+bci(j-1,kk,2))*qi0(j-1,kk,5,3)/qi0(j-1,kk,1,3) +
     .       bci(j-1,kk,2)*q(j-1,kk,i-1,5)/q(j-1,kk,i-1,1) +
     .       q(j+1,kk,i-1,5)/q(j+1,kk,i-1,1) -
     .       q(j-1,kk,i-1,5)/q(j-1,kk,i-1,1) )

           t(jiv+j,4) = 0.25*(
     .       (1.+bci(j+1,kk,2))*qi0(j+1,kk,2,3) -
     .       bci(j+1,kk,2)*q(j+1,kk,i-1,2) -
     .       (1.+bci(j-1,kk,2))*qi0(j-1,kk,2,3) +
     .       bci(j-1,kk,2)*q(j-1,kk,i-1,2) +
     .       q(j+1,kk,i-1,2) -
     .       q(j-1,kk,i-1,2) )

           t(jiv+j,5) = 0.25*(
     .       (1.+bci(j+1,kk,2))*qi0(j+1,kk,3,3) -
     .       bci(j+1,kk,2)*q(j+1,kk,i-1,3) -
     .       (1.+bci(j-1,kk,2))*qi0(j-1,kk,3,3) +
     .       bci(j-1,kk,2)*q(j-1,kk,i-1,3) +
     .       q(j+1,kk,i-1,3) -
     .       q(j-1,kk,i-1,3) )

           t(jiv+j,6) = 0.25*(
     .       (1.+bci(j+1,kk,2))*qi0(j+1,kk,4,3) -
     .       bci(j+1,kk,2)*q(j+1,kk,i-1,4) -
     .       (1.+bci(j-1,kk,2))*qi0(j-1,kk,4,3) +
     .       bci(j-1,kk,2)*q(j-1,kk,i-1,4) +
     .       q(j+1,kk,i-1,4) -
     .       q(j-1,kk,i-1,4) )

          end if
  227   continue
      else
c
c-- general case for i between 1 and idim
c
        do 229 j = 1,jdim1
          if (j. eq. 1) then
c
c replace (j-1) data by boundary data
c
           t(jiv+j,3) = 0.25*gprgm1*(
     .       q(j+1,kk,i,5)/q(j+1,kk,i,1) -
     .       (1.+bcj(kk,i,1))*qj0(kk,i,5,1)/qj0(kk,i,1,1) +
     .       bcj(kk,i,1)*q(j,kk,i,5)/q(j,kk,i,1) +
     .       q(j+1,kk,i-1,5)/q(j+1,kk,i-1,1) -
     .       (1.+bcj(kk,i-1,1))*qj0(kk,i-1,5,1)/qj0(kk,i-1,1,1) +
     .       bcj(kk,i-1,1)*q(j,kk,i-1,5)/q(j,kk,i-1,1) )

           t(jiv+j,4) = 0.25*(
     .       q(j+1,kk,i,2) -
     .       (1.+bcj(kk,i,1))*qj0(kk,i,2,1) +
     .       bcj(kk,i,1)*q(j,kk,i,2) +
     .       q(j+1,kk,i-1,2) -
     .       (1.+bcj(kk,i-1,1))*qj0(kk,i-1,2,1) +
     .       bcj(kk,i-1,1)*q(j,kk,i-1,2) )

           t(jiv+j,5) = 0.25*(
     .       q(j+1,kk,i,3) -
     .       (1.+bcj(kk,i,1))*qj0(kk,i,3,1) +
     .       bcj(kk,i,1)*q(j,kk,i,3) +
     .       q(j+1,kk,i-1,3) -
     .       (1.+bcj(kk,i-1,1))*qj0(kk,i-1,3,1) +
     .       bcj(kk,i-1,1)*q(j,kk,i-1,3) )

           t(jiv+j,6) = 0.25*(
     .       q(j+1,kk,i,4) -
     .       (1.+bcj(kk,i,1))*qj0(kk,i,4,1) +
     .       bcj(kk,i,1)*q(j,kk,i,4) +
     .       q(j+1,kk,i-1,4) -
     .       (1.+bcj(kk,i-1,1))*qj0(kk,i-1,4,1) +
     .       bcj(kk,i-1,1)*q(j,kk,i-1,4) )

          else if (j .eq. jdim1) then
c
c replace (j+1) data by boundary data 
c
           t(jiv+j,3) = 0.25*gprgm1*(
     .       (1.+bcj(kk,i,2))*qj0(kk,i,5,3)/qj0(kk,i,1,3) -
     .       bcj(kk,i,2)*q(j,kk,i,5)/q(j,kk,i,1) - 
     .       q(j-1,kk,i,5)/q(j-1,kk,i,1)  +
     .       (1.+bcj(kk,i-1,2))*qj0(kk,i-1,5,3)/qj0(kk,i-1,1,3) -
     .       bcj(kk,i-1,2)*q(j-1,kk,i,5)/q(j-1,kk,i,1) - 
     .       q(j-1,kk,i-1,5)/q(j-1,kk,i-1,1) )

           t(jiv+j,4) = 0.25*(
     .       (1.+bcj(kk,i,2))*qj0(kk,i,2,3) -
     .       bcj(kk,i,2)*q(j,kk,i,2) - 
     .       q(j-1,kk,i,2)  +
     .       (1.+bcj(kk,i-1,2))*qj0(kk,i-1,2,3) -
     .       bcj(kk,i-1,2)*q(j-1,kk,i,2) - 
     .       q(j-1,kk,i-1,2) )

           t(jiv+j,5) = 0.25*(
     .       (1.+bcj(kk,i,2))*qj0(kk,i,3,3) -
     .       bcj(kk,i,2)*q(j,kk,i,3) - 
     .       q(j-1,kk,i,3)  +
     .       (1.+bcj(kk,i-1,2))*qj0(kk,i-1,3,3) -
     .       bcj(kk,i-1,2)*q(j-1,kk,i,3) - 
     .       q(j-1,kk,i-1,3) )

           t(jiv+j,6) = 0.25*(
     .       (1.+bcj(kk,i,2))*qj0(kk,i,4,3) -
     .       bcj(kk,i,2)*q(j,kk,i,4) - 
     .       q(j-1,kk,i,4)  +
     .       (1.+bcj(kk,i-1,2))*qj0(kk,i-1,4,3) -
     .       bcj(kk,i-1,2)*q(j-1,kk,i,4) - 
     .       q(j-1,kk,i-1,4) )
          else
c
c-- general case
c
            t(jiv+j,3) = 0.25*gprgm1*(
     .        q(j+1,kk,i,5)/q(j+1,kk,i,1) -
     .        q(j-1,kk,i,5)/q(j-1,kk,i,1) +
     .        q(j+1,kk,i-1,5)/q(j+1,kk,i-1,1) -
     .        q(j-1,kk,i-1,5)/q(j-1,kk,i-1,1) )

            t(jiv+j,4) = 0.25*(
     .        q(j+1,kk,i,2)   - q(j-1,kk,i,2) +
     .        q(j+1,kk,i-1,2) - q(j-1,kk,i-1,2) )

            t(jiv+j,5) = 0.25*(
     .        q(j+1,kk,i,3)   - q(j-1,kk,i,3) +
     .        q(j+1,kk,i-1,3) - q(j-1,kk,i-1,3) )

            t(jiv+j,6) = 0.25*(
     .        q(j+1,kk,i,4)   - q(j-1,kk,i,4) +
     .        q(j+1,kk,i-1,4) - q(j-1,kk,i-1,4) )

          end if
  229   continue
      end if
  230 continue
c
c******************************************************************************
c
c   gradients at cell interfaces (needed for full ns terms). Note that 
c  we don't need derivatives in xie direction, because they have already
c  been accounted for in ffluxv.
c  Recall bck=0 => cell data, bck=1 => face data
c
c  t(7)       : d( c*c/(gm1*pr) )/d(zeta)
c  t(8)       : d(u)/d(zeta)
c  t(9)       : d(v)/d(zeta)
c  t(10)      : d(w)/d(zeta)
c
c     interior interfaces
c
      do 380 kpl = 1,npl
      kk = k + kpl -1
      if (kk .eq. 1) then
c
c replace (kk-1) data by boundary data
c
        do 330 i = 1,idim
          jiv     = (i-1)*npl*jdim1 + (kpl-1)*jdim1 
c
          if (i .eq. 1) then
c
c replace (i-1) data by boundary data
c for missing data at corners use one-sided differencing of known boundary values
c
            do 323 j = 1,jdim1

            t(jiv+j,7) = 0.25*gprgm1*(
     .        q(j,kk+1,i,5)/q(j,kk+1,i,1) -
     .        (1.+bck(j,i,1))*qk0(j,i,5,1)/qk0(j,i,1,1) +
     .        bck(j,i,1)*q(j,kk,i,5)/q(j,kk,i,1) +
     .        2.*(1.+bci(j,kk+1,1))*qi0(j,kk+1,5,1)/qi0(j,kk+1,1,1) -
     .        2.*bci(j,kk+1,1)*q(j,kk+1,i,5)/q(j,kk+1,i,1) -
     .        2.*(1.+bci(j,kk,1))*qi0(j,kk,5,1)/qi0(j,kk,1,1) +
     .        2.*bci(j,kk,1)*q(j,kk,i,5)/q(j,kk,i,1) )

            t(jiv+j,8) = 0.25*(
     .        q(j,kk+1,i,2) -
     .        (1.+bck(j,i,1))*qk0(j,i,2,1) +
     .        bck(j,i,1)*q(j,kk,i,2) +
     .        2.*(1.+bci(j,kk+1,1))*qi0(j,kk+1,2,1) -
     .        2.*bci(j,kk+1,1)*q(j,kk+1,i,2) -
     .        2.*(1.+bci(j,kk,1))*qi0(j,kk,2,1) +
     .        2.*bci(j,kk,1)*q(j,kk,i,2) )

            t(jiv+j,9) = 0.25*(
     .        q(j,kk+1,i,3) -
     .        (1.+bck(j,i,1))*qk0(j,i,3,1) +
     .        bck(j,i,1)*q(j,kk,i,3) +
     .        2.*(1.+bci(j,kk+1,1))*qi0(j,kk+1,3,1) -
     .        2.*bci(j,kk+1,1)*q(j,kk+1,i,3) -
     .        2.*(1.+bci(j,kk,1))*qi0(j,kk,3,1) +
     .        2.*bci(j,kk,1)*q(j,kk,i,3) )

            t(jiv+j,10) = 0.25*(
     .        q(j,kk+1,i,4) -
     .        (1.+bck(j,i,1))*qk0(j,i,4,1) +
     .        bck(j,i,1)*q(j,kk,i,4) +
     .        2.*(1.+bci(j,kk+1,1))*qi0(j,kk+1,4,1) -
     .        2.*bci(j,kk+1,1)*q(j,kk+1,i,4) -
     .        2.*(1.+bci(j,kk,1))*qi0(j,kk,4,1) +
     .        2.*bci(j,kk,1)*q(j,kk,i,4) )

  323       continue
          else if (i .eq. idim) then
c
c replace i data by boundary data
c for missing data at corners use one-sided differencing of known boundary values
c
            do 325 j = 1,jdim1

            t(jiv+j,7) = 0.25*gprgm1*(
     .        2.*(1.+bci(j,kk+1,2))*qi0(j,kk+1,5,3)/qi0(j,kk+1,1,3) -
     .        2.*bci(j,kk+1,2)*q(j,kk+1,i-1,5)/q(j,kk+1,i-1,1) -
     .        2.*(1.+bci(j,kk,2))*qi0(j,kk,5,3)/qi0(j,kk,1,3) +
     .        2.*bci(j,kk,2)*q(j,kk,i-1,5)/q(j,kk,i-1,1) +
     .        q(j,kk+1,i-1,5)/q(j,kk+1,i-1,1) -
     .        (1.+bck(j,i-1,1))*qk0(j,i-1,5,1)/qk0(j,i-1,1,1) +
     .        bck(j,i-1,1)*q(j,kk,i-1,5)/q(j,kk,i-1,1) )

            t(jiv+j,8) = 0.25*(
     .        2.*(1.+bci(j,kk+1,2))*qi0(j,kk+1,2,3) -
     .        2.*bci(j,kk+1,2)*q(j,kk+1,i-1,2) -
     .        2.*(1.+bci(j,kk,2))*qi0(j,kk,2,3) +
     .        2.*bci(j,kk,2)*q(j,kk,i-1,2) +
     .        q(j,kk+1,i-1,2) -
     .        (1.+bck(j,i-1,1))*qk0(j,i-1,2,1) +
     .        bck(j,i-1,1)*q(j,kk,i-1,2) )

            t(jiv+j,9) = 0.25*(
     .        2.*(1.+bci(j,kk+1,2))*qi0(j,kk+1,3,3) -
     .        2.*bci(j,kk+1,2)*q(j,kk+1,i-1,3) -
     .        2.*(1.+bci(j,kk,2))*qi0(j,kk,3,3) +
     .        2.*bci(j,kk,2)*q(j,kk,i-1,3) +
     .        q(j,kk+1,i-1,3) -
     .        (1.+bck(j,i-1,1))*qk0(j,i-1,3,1) +
     .        bck(j,i-1,1)*q(j,kk,i-1,3) )

            t(jiv+j,10) = 0.25*(
     .        2.*(1.+bci(j,kk+1,2))*qi0(j,kk+1,4,3) -
     .        2.*bci(j,kk+1,2)*q(j,kk+1,i-1,4) -
     .        2.*(1.+bci(j,kk,2))*qi0(j,kk,4,3) +
     .        2.*bci(j,kk,2)*q(j,kk,i-1,4) +
     .        q(j,kk+1,i-1,4) -
     .        (1.+bck(j,i-1,1))*qk0(j,i-1,4,1) +
     .        bck(j,i-1,1)*q(j,kk,i-1,4) )

  325       continue
          else
c
c-- general case for kk = 1
c
            do 327 j = 1,jdim1

            t(jiv+j,7) = 0.25*gprgm1*(
     .        q(j,kk+1,i,5)/q(j,kk+1,i,1) -
     .        (1.+bck(j,i,1))*qk0(j,i,5,1)/qk0(j,i,1,1) +
     .        bck(j,i,1)*q(j,kk,i,5)/q(j,kk,i,1) +
     .        q(j,kk+1,i-1,5)/q(j,kk+1,i-1,1) -
     .        (1.+bck(j,i-1,1))*qk0(j,i-1,5,1)/qk0(j,i-1,1,1) +
     .        bck(j,i-1,1)*q(j,kk,i-1,5)/q(j,kk,i-1,1) )

            t(jiv+j,8) = 0.25*(
     .        q(j,kk+1,i,2) -
     .        (1.+bck(j,i,1))*qk0(j,i,2,1) +
     .        bck(j,i,1)*q(j,kk,i,2) +
     .        q(j,kk+1,i-1,2) -
     .        (1.+bck(j,i-1,1))*qk0(j,i-1,2,1) +
     .        bck(j,i-1,1)*q(j,kk,i-1,2) )

            t(jiv+j,9) = 0.25*(
     .        q(j,kk+1,i,3) -
     .        (1.+bck(j,i,1))*qk0(j,i,3,1) +
     .        bck(j,i,1)*q(j,kk,i,3) +
     .        q(j,kk+1,i-1,3) -
     .        (1.+bck(j,i-1,1))*qk0(j,i-1,3,1) +
     .        bck(j,i-1,1)*q(j,kk,i-1,3) )

            t(jiv+j,10) = 0.25*(
     .        q(j,kk+1,i,4) -
     .        (1.+bck(j,i,1))*qk0(j,i,4,1) +
     .        bck(j,i,1)*q(j,kk,i,4) +
     .        q(j,kk+1,i-1,4) -
     .        (1.+bck(j,i-1,1))*qk0(j,i-1,4,1) +
     .        bck(j,i-1,1)*q(j,kk,i-1,4) )

  327       continue
          end if
  330   continue

      else if (kk .eq. kdim1) then
c
c replace (kk+1) data by boundary data
c
        do 340 i = 1,idim
          jiv     = (i-1)*npl*jdim1 + (kpl-1)*jdim1 
          if (i .eq. 1) then
c
c replace (i-1) data by boundary data
c for missing data at corners use one-sided differencing of known boundary values
c
            do 333 j = 1,jdim1

            t(jiv+j,7) = 0.25*gprgm1*(
     .        (1.+bck(j,i,2))*qk0(j,i,5,3)/qk0(j,i,1,3) -
     .        bck(j,i,2)*q(j,kk,i,5)/q(j,kk,i,1) -
     .        q(j,kk-1,i,5)/q(j,kk-1,i,1) +
     .        2.*(1.+bci(j,kk,1))*qi0(j,kk,5,1)/qi0(j,kk,1,1) -
     .        2.*bci(j,kk,1)*q(j,kk,i,5)/q(j,kk,i,1) -
     .        2.*(1.+bci(j,kk-1,1))*qi0(j,kk-1,5,1)/qi0(j,kk-1,1,1) +
     .        2.*bci(j,kk-1,1)*q(j,kk-1,i,5)/q(j,kk-1,i,1) )

            t(jiv+j,8) = 0.25*(
     .        (1.+bck(j,i,2))*qk0(j,i,2,3) -
     .        bck(j,i,2)*q(j,kk,i,2) -
     .        q(j,kk-1,i,2) +
     .        2.*(1.+bci(j,kk,1))*qi0(j,kk,2,1) -
     .        2.*bci(j,kk,1)*q(j,kk,i,2) -
     .        2.*(1.+bci(j,kk-1,1))*qi0(j,kk-1,2,1) +
     .        2.*bci(j,kk-1,1)*q(j,kk-1,i,2) )

            t(jiv+j,9) = 0.25*(
     .        (1.+bck(j,i,2))*qk0(j,i,3,3) -
     .        bck(j,i,2)*q(j,kk,i,3) -
     .        q(j,kk-1,i,3) +
     .        2.*(1.+bci(j,kk,1))*qi0(j,kk,3,1) -
     .        2.*bci(j,kk,1)*q(j,kk,i,3) -
     .        2.*(1.+bci(j,kk-1,1))*qi0(j,kk-1,3,1) +
     .        2.*bci(j,kk-1,1)*q(j,kk-1,i,3) )

            t(jiv+j,10) = 0.25*(
     .        (1.+bck(j,i,2))*qk0(j,i,4,3) -
     .        bck(j,i,2)*q(j,kk,i,4) -
     .        q(j,kk-1,i,4) +
     .        2.*(1.+bci(j,kk,1))*qi0(j,kk,4,1) -
     .        2.*bci(j,kk,1)*q(j,kk,i,4) -
     .        2.*(1.+bci(j,kk-1,1))*qi0(j,kk-1,4,1) +
     .        2.*bci(j,kk-1,1)*q(j,kk-1,i,4) )

  333       continue
          else if (i .eq. idim) then
c
c replace i data by boundary data
c for missing data at corners use one-sided differencing of known boundary values
c
            do 335 j = 1,jdim1

            t(jiv+j,7) = 0.25*gprgm1*(
     .        2.*(1.+bci(j,kk,2))*qi0(j,kk,5,3)/qi0(j,kk,1,3) -
     .        2.*bci(j,kk,2)*q(j,kk,i-1,5)/q(j,kk,i-1,1) -
     .        2.*(1.+bci(j,kk-1,2))*qi0(j,kk-1,5,3)/qi0(j,kk-1,1,3) +
     .        2.*bci(j,kk-1,2)*q(j,kk-1,i-1,5)/q(j,kk-1,i-1,1) +
     .        (1.+bck(j,i-1,2))*qk0(j,i-1,5,3)/qk0(j,i-1,1,3) -
     .        bck(j,i-1,2)*q(j,kk,i-1,5)/q(j,kk,i-1,1) -
     .        q(j,kk-1,i-1,5)/q(j,kk-1,i-1,1) )

            t(jiv+j,8) = 0.25*(
     .        2.*(1.+bci(j,kk,2))*qi0(j,kk,2,3) -
     .        2.*bci(j,kk,2)*q(j,kk,i-1,2) -
     .        2.*(1.+bci(j,kk-1,2))*qi0(j,kk-1,2,3) +
     .        2.*bci(j,kk-1,2)*q(j,kk-1,i-1,2) +
     .        (1.+bck(j,i-1,2))*qk0(j,i-1,2,3) -
     .        bck(j,i-1,2)*q(j,kk,i-1,2) -
     .        q(j,kk-1,i-1,2) )

            t(jiv+j,9) = 0.25*(
     .        2.*(1.+bci(j,kk,2))*qi0(j,kk,3,3) -
     .        2.*bci(j,kk,2)*q(j,kk,i-1,3) -
     .        2.*(1.+bci(j,kk-1,2))*qi0(j,kk-1,3,3) +
     .        2.*bci(j,kk-1,2)*q(j,kk-1,i-1,3) +
     .        (1.+bck(j,i-1,2))*qk0(j,i-1,3,3) -
     .        bck(j,i-1,2)*q(j,kk,i-1,3) -
     .        q(j,kk-1,i-1,3) )

            t(jiv+j,10) = 0.25*(
     .        2.*(1.+bci(j,kk,2))*qi0(j,kk,4,3) -
     .        2.*bci(j,kk,2)*q(j,kk,i-1,4) -
     .        2.*(1.+bci(j,kk-1,2))*qi0(j,kk-1,4,3) +
     .        2.*bci(j,kk-1,2)*q(j,kk-1,i-1,4) +
     .        (1.+bck(j,i-1,2))*qk0(j,i-1,4,3) -
     .        bck(j,i-1,2)*q(j,kk,i-1,4) -
     .        q(j,kk-1,i-1,4) )

  335       continue
          else
c
c-- general case for kk = kdim1
c
            do 337 j = 1,jdim1

            t(jiv+j,7) = 0.25*gprgm1*(
     .        (1.+bck(j,i,2))*qk0(j,i,5,3)/qk0(j,i,1,3) -
     .        bck(j,i,2)*q(j,kk,i,5)/q(j,kk,i,1) -
     .        q(j,kk-1,i,5)/q(j,kk-1,i,1) +
     .        (1.+bck(j,i-1,2))*qk0(j,i-1,5,3)/qk0(j,i-1,1,3) -
     .        bck(j,i-1,2)*q(j,kk,i-1,5)/q(j,kk,i-1,1) -
     .        q(j,kk-1,i-1,5)/q(j,kk-1,i-1,1) )

            t(jiv+j,8) = 0.25*(
     .        (1.+bck(j,i,2))*qk0(j,i,2,3) -
     .        bck(j,i,2)*q(j,kk,i,2) -
     .        q(j,kk-1,i,2) +
     .        (1.+bck(j,i-1,2))*qk0(j,i-1,2,3) -
     .        bck(j,i-1,2)*q(j,kk,i-1,2) -
     .        q(j,kk-1,i-1,2) )

            t(jiv+j,9) = 0.25*(
     .        (1.+bck(j,i,2))*qk0(j,i,3,3) -
     .        bck(j,i,2)*q(j,kk,i,3) -
     .        q(j,kk-1,i,3) +
     .        (1.+bck(j,i-1,2))*qk0(j,i-1,3,3) -
     .        bck(j,i-1,2)*q(j,kk,i-1,3) -
     .        q(j,kk-1,i-1,3) )

            t(jiv+j,10) = 0.25*(
     .        (1.+bck(j,i,2))*qk0(j,i,4,3) -
     .        bck(j,i,2)*q(j,kk,i,4) -
     .        q(j,kk-1,i,4) +
     .        (1.+bck(j,i-1,2))*qk0(j,i-1,4,3) -
     .        bck(j,i-1,2)*q(j,kk,i-1,4) -
     .        q(j,kk-1,i-1,4) )

  337       continue
          end if
  340   continue
      else
c
c-- general case (kk .ne 1. and kk .ne. kdkm1)
c
        do 350 i = 1,idim
          jiv     = (i-1)*npl*jdim1 + (kpl-1)*jdim1 
          if (i .eq. 1) then
c
c replace (i-1) data by boundary data
c
            do 343 j = 1,jdim1

            t(jiv+j,7)  = 0.25*gprgm1*(
     .        q(j,kk+1,i,5)/q(j,kk+1,i,1) -
     .        q(j,kk-1,i,5)/q(j,kk-1,i,1) +
     .        (1.+bci(j,kk+1,1))*qi0(j,kk+1,5,1)/qi0(j,kk+1,1,1) -
     .        bci(j,kk+1,1)*q(j,kk+1,i,5)/q(j,kk+1,i,1) -
     .        (1.+bci(j,kk-1,1))*qi0(j,kk-1,5,1)/qi0(j,kk-1,1,1) +
     .        bci(j,kk-1,1)*q(j,kk-1,i,5)/q(j,kk-1,i,1) )

            t(jiv+j,8)  = 0.25*(
     .        q(j,kk+1,i,2) -q(j,kk-1,i,2) +
     .        (1.+bci(j,kk+1,1))*qi0(j,kk+1,2,1) -
     .        bci(j,kk+1,1)*q(j,kk+1,i,2) -
     .        (1.+bci(j,kk-1,1))*qi0(j,kk-1,2,1) +
     .        bci(j,kk-1,1)*q(j,kk-1,i,2) )

            t(jiv+j,9)  = 0.25*(
     .        q(j,kk+1,i,3) -q(j,kk-1,i,3) +
     .        (1.+bci(j,kk+1,1))*qi0(j,kk+1,3,1) -
     .        bci(j,kk+1,1)*q(j,kk+1,i,3) -
     .        (1.+bci(j,kk-1,1))*qi0(j,kk-1,3,1) +
     .        bci(j,kk-1,1)*q(j,kk-1,i,3) )

            t(jiv+j,10) = 0.25*(
     .        q(j,kk+1,i,4) -q(j,kk-1,i,4) +
     .        (1.+bci(j,kk+1,1))*qi0(j,kk+1,4,1) -
     .        bci(j,kk+1,1)*q(j,kk+1,i,4) -
     .        (1.+bci(j,kk-1,1))*qi0(j,kk-1,4,1) +
     .        bci(j,kk-1,1)*q(j,kk-1,i,4) )

  343       continue
          else if (i .eq. idim) then
c
c replace i data by boundary data
c
            do 345 j = 1,jdim1

            t(jiv+j,7)  = 0.25*gprgm1*(
     .        (1.+bci(j,kk+1,2))*qi0(j,kk+1,5,3)/qi0(j,kk+1,1,3) -
     .        bci(j,kk+1,2)*q(j,kk+1,i-1,5)/q(j,kk+1,i-1,1) -
     .        (1.+bci(j,kk-1,2))*qi0(j,kk-1,5,3)/qi0(j,kk-1,1,3) +
     .        bci(j,kk-1,2)*q(j,kk-1,i-1,5)/q(j,kk-1,i-1,1) +
     .        q(j,kk+1,i-1,5)/q(j,kk+1,i-1,1) -
     .        q(j,kk-1,i-1,5)/q(j,kk-1,i-1,1) )

            t(jiv+j,8)  = 0.25*(
     .        (1.+bci(j,kk+1,2))*qi0(j,kk+1,2,3) -
     .        bci(j,kk+1,2)*q(j,kk+1,i-1,2) -
     .        (1.+bci(j,kk-1,2))*qi0(j,kk-1,2,3) +
     .        bci(j,kk-1,2)*q(j,kk-1,i-1,2) +
     .        q(j,kk+1,i-1,2) -q(j,kk-1,i-1,2) )

            t(jiv+j,9)  = 0.25*(
     .        (1.+bci(j,kk+1,2))*qi0(j,kk+1,3,3) -
     .        bci(j,kk+1,2)*q(j,kk+1,i-1,3) -
     .        (1.+bci(j,kk-1,2))*qi0(j,kk-1,3,3) +
     .        bci(j,kk-1,2)*q(j,kk-1,i-1,3) +
     .        q(j,kk+1,i-1,3) -q(j,kk-1,i-1,3) )

            t(jiv+j,10) = 0.25*(
     .        (1.+bci(j,kk+1,2))*qi0(j,kk+1,4,3) -
     .        bci(j,kk+1,2)*q(j,kk+1,i-1,4) -
     .        (1.+bci(j,kk-1,2))*qi0(j,kk-1,4,3) +
     .        bci(j,kk-1,2)*q(j,kk-1,i-1,4) +
     .        q(j,kk+1,i-1,4) -q(j,kk-1,i-1,4) )

  345       continue
          else
            do 347 j = 1,jdim1
c
c-- general case
c
            t(jiv+j,7)  = 0.25*gprgm1*(
     .        q(j,kk+1,i,5)/q(j,kk+1,i,1) -
     .        q(j,kk-1,i,5)/q(j,kk-1,i,1) +
     .        q(j,kk+1,i-1,5)/q(j,kk+1,i-1,1) -
     .        q(j,kk-1,i-1,5)/q(j,kk-1,i-1,1) )

            t(jiv+j,8)  = 0.25*(
     .        q(j,kk+1,i,2)  -q(j,kk-1,i,2) +
     .        q(j,kk+1,i-1,2)-q(j,kk-1,i-1,2) )

            t(jiv+j,9)  = 0.25*(
     .        q(j,kk+1,i,3)  -q(j,kk-1,i,3) +
     .        q(j,kk+1,i-1,3)-q(j,kk-1,i-1,3) )

            t(jiv+j,10) = 0.25*(
     .        q(j,kk+1,i,4)  -q(j,kk-1,i,4) +
     .        q(j,kk+1,i-1,4)-q(j,kk-1,i-1,4) )

  347       continue
          end if
  350   continue
      end if
  380 continue
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
      do 1024 j=1,nn
      t(j+js-1,11) = 0.5e0*(t(j+js-1,18)+t(j,18))
      t(j+js-1,12) = 0.5e0*(t(j+js-1,19)+t(j,19))
      t(j+js-1,13) = 0.5e0*(t(j+js-1,20)+t(j,20))
 1024 continue
c     interfaces at i=0 and i=idim
      do 1025 izz=1,jv
      ab=1.+wi0(izz,15)
      bb=1.-wi0(izz,15)
      t(izz+n,11)=0.5*(ab*wi0(izz,18)+bb*t(izz+n-jv,18))
      t(izz+n,12)=0.5*(ab*wi0(izz,19)+bb*t(izz+n-jv,19))
      t(izz+n,13)=0.5*(ab*wi0(izz,20)+bb*t(izz+n-jv,20))
      ab=1.+wi0(izz,5)
      bb=1.-wi0(izz,5)
      t(izz,11) = 0.5*(ab*wi0(izz,8) +bb*t(izz,18))
      t(izz,12) = 0.5*(ab*wi0(izz,9) +bb*t(izz,19))
      t(izz,13) = 0.5*(ab*wi0(izz,10) +bb*t(izz,20))
 1025 continue
c
c******************************************************************************
c  calculate fluxes
c
c  viscous terms at interfaces
cdir$ ivdep
      do 1026 j=1,l0
c
c-- form ux, uy, uz, vx, vy, vz, wx, wy, wz, tx, ty and tz
c   excluding contributions in i/xie direction (already accounted for in
c   thin-layer approximation. Note that t = c*c/(gm1*pr)
c
       ux = t(j,4)*t(j,28) + t(j,8)*t(j,31)
       uy = t(j,4)*t(j,29) + t(j,8)*t(j,32)
       uz = t(j,4)*t(j,30) + t(j,8)*t(j,2)

       vx = t(j,5)*t(j,28) + t(j,9)*t(j,31)
       vy = t(j,5)*t(j,29) + t(j,9)*t(j,32)
       vz = t(j,5)*t(j,30) + t(j,9)*t(j,2)

       wx = t(j,6)*t(j,28) + t(j,10)*t(j,31)
       wy = t(j,6)*t(j,29) + t(j,10)*t(j,32)
       wz = t(j,6)*t(j,30) + t(j,10)*t(j,2)

       tx = t(j,3)*t(j,28) + t(j,7)*t(j,31)
       ty = t(j,3)*t(j,29) + t(j,7)*t(j,32)
       tz = t(j,3)*t(j,30) + t(j,7)*t(j,2)
c
c--  Note lambda = -2/3 mu, lambda + 2 mu = 4/3 mu
c--  u-momentum flux
c
      t(j,16)      = t(j,14)*(  t(j,25)*(4./3.*ux-2./3.*(vy + wz))
     .                            + t(j,26)*(uy + vx)
     .                            + t(j,27)*(uz + wx))

c
c--  v-momentum flux
c
      t(j,17)      = t(j,14)*(  t(j,25)*(vx + uy)
     .                            + t(j,26)*(4./3.*vy-2./3.*(ux + wz))
     .                            + t(j,27)*(vz + wy))

c
c--  w-momentum flux
c
      t(j,18)      = t(j,14)*(  t(j,25)*(wx+uz)
     .                            + t(j,26)*(wy+vz)
     .                            + t(j,27)*(4./3.*wz-2./3.*(vy+ux)))
c
c--  energy flux
c
      t(j,19) = t(j,14)*(
     .              t(j,25)*((4./3.*ux - 2./3*(vy+wz))*t(j,11) +
     .                         (uy+vx)*t(j,12) + (uz+wx)*t(j,13)) +
     .              t(j,26)*((vx+uy)*t(j,11) + (vz+wy)*t(j,13) +
     .                         (4./3.*vy - 2./3*(ux+wz))*t(j,12)) +
     .              t(j,27)*((wx+uz)*t(j,11) + (wy+vz)*t(j,12) +
     .                       (4./3.*wz - 2./3*(ux+vy))*t(j,13)) +
     .          t(j,23)*(t(j,25)*tx + t(j,26)*ty + t(j,27)*tz))
c
 1026 continue
c
c  (-)viscous flux =  Fv(i-1/2)-Fv(i+1/2) (momentum eqs.)
c
      do 1027 j=1,n
      t(j,2)  = -t(j+js-1,16) + t(j,16)
      t(j,3)  = -t(j+js-1,17) + t(j,17)
      t(j,4)  = -t(j+js-1,18) + t(j,18)
      t(j,5)  = -t(j+js-1,19) + t(j,19)
 1027 continue
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
      return
      end
