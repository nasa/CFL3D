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
      subroutine csurf(jdim,kdim,idim,x,y,z,sk,sj,si,q,ics,
     .                 ub,vb,wb,vmuk,vmuj,vmui,vol, bcj,bck,bci,
     .                 blank,xtbj,xtbk,xtbi,iuns,ncs,icsinfo,
     .                 sx,sy,sz,stot,pav,ptav,tav,ttav,xmav,fmdot,
     .                 cfxp,cfyp,cfzp,cfdp,cflp,cftp,cfxv,cfyv,cfzv,
     .                 cfdv,cflv,cftv,cfxmom,cfymom,cfzmom,cfdmom,
     .                 cflmom,cftmom,cfxtot,cfytot,cfztot,cfdtot,
     .                 cfltot,cfttot,maxcs,qj0,qk0,qi0)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Integrate control surface mass flow and momemtum/forces
c               Original coding by R. Cedar, GE Aircraft Engines
c
c               Modified 3/01 by R. Biedron to use either cell-face
c               values or cell-center values at zonal boundaries, 
c               depending on bc flag (bcj/bck/bci).
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension icsinfo(maxcs,9) 
      dimension ub(jdim,kdim,idim),vb(jdim,kdim,idim),wb(jdim,kdim,idim)
      dimension vmuk(jdim-1,idim-1,2),vol(jdim,kdim,idim-1),
     .          vmuj(kdim-1,idim-1,2),vmui(jdim-1,kdim-1,2)
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
      dimension xtbj(kdim,idim-1,3,2),xtbk(jdim,idim-1,3,2),
     .          xtbi(jdim,kdim,3,2)

      dimension sk(jdim,kdim,idim-1,5),sj(jdim,kdim,idim-1,5),
     .          si(jdim,kdim,idim,5)
      dimension q(jdim,kdim,idim,5),blank(jdim,kdim,idim)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
      dimension qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4),
     .          qi0(jdim,kdim,5,4)
c
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /fsum/ sref,cref,bref,xmc,ymc,zmc
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .        tur10(7)
      common /reyue/ reue,tinf,ivisc(3)
c
c**************************************************************************
c
c     sign conventions for forces:
c       cfx.....x-component of force
c              positive for force in +x direction
c       cfy.....y-component of force
c              positive for force in +y direction
c       cfz.....z-component of force
c              positive for force in +z direction
c
c     Note:- The viscous components of the forces will always be positve
c            in the respective directions, but the pressure and momentum
c            terms will only be positive if the normal direction of the
c            surface (fnorm) is correctly set in the input file.
c
c**************************************************************************
c
      ist   = icsinfo(ics,2)
      ifn   = icsinfo(ics,3) - 1
      jst   = icsinfo(ics,4)
      jfn   = icsinfo(ics,5) - 1
      kst   = icsinfo(ics,6)
      kfn   = icsinfo(ics,7) - 1
      iwall = icsinfo(ics,8)
      fnorm = float(icsinfo(ics,9))
      if (icsinfo(ics,9).eq.0) fnorm = 1.0
c
      cpc   = 2.e0/(gamma*xmach*xmach)
      const = 4./(reue*xmach)
c
      cosa  = cos(alpha)
      sina  = sin(alpha)
      cosb  = cos(beta)
      sinb  = sin(beta)
c
      sx     = 0.e0
      sy     = 0.e0
      sz     = 0.e0
      pav    = 0.e0
      ptav   = 0.e0
      tav    = 0.e0
      ttav   = 0.e0
      xmav   = 0.e0
      fmdot  = 0.e0
      cfxp   = 0.e0
      cfyp   = 0.e0 
      cfzp   = 0.e0
      cfxv   = 0.e0
      cfyv   = 0.e0 
      cfzv   = 0.e0
      cfxmom = 0.e0
      cfymom = 0.e0
      cfzmom = 0.e0
c
c******************************************************************
c     forces on k=constant surfaces
c******************************************************************
c
      if (kst.eq.kfn+1) then
         k    = kst
         kc   = k
         kcm1 = k - 1
         m = 0
c
         if (k.eq.1) then
            k    = k
            kc   = 1
            kcm1 = 1
            m    = 1
            mm   = 1
         end if
c
         if (k.eq.kdim) then
            k    = kdim
            kc   = kdim - 1
            kcm1 = kdim - 1
            m    = 2
            mm   = 3
         end if
c
         if (k.eq.1 .and. iwall.eq.1) then
            kk   = 1     
         end if
c
         if (k.eq.kdim .and. iwall.eq.1) then
            kk   = 2   
         end if
c
         do 10 i=ist,ifn
         do 15 j=jst,jfn
c
c..      area integration
c
         sx    = sx + sk(j,k,i,1)*sk(j,k,i,4)
         sy    = sy + sk(j,k,i,2)*sk(j,k,i,4)
         sz    = sz + sk(j,k,i,3)*sk(j,k,i,4)
c
c..      pressure force integration
c
         if (m .ne. 0) then
            aa = 1. +  bck(j,i,m)
            bb = 1. -  bck(j,i,m)
            dcp  = -(0.5*(aa*qk0(j,i,5,mm) + bb*q(j,kc,i,5))/p0
     .             -1.e0)*cpc*sk(j,k,i,4)
         else
            dcp  = -(0.5*(q(j,kc,i,5)+q(j,kcm1,i,5))/p0
     .             -1.e0)*cpc*sk(j,k,i,4)
         end if
         cfxp =  cfxp - dcp*sk(j,k,i,1)
         cfyp =  cfyp - dcp*sk(j,k,i,2)
         cfzp =  cfzp - dcp*sk(j,k,i,3)
c
c..      mass and momentum integration
c
         if (iwall.eq.0) then
            if (m .ne. 0) then
               aa = 1. +  bck(j,i,m)
               bb = 1. -  bck(j,i,m)
               q1 = 0.5*(aa*qk0(j,i,1,mm)
     .                        +      bb*q(j,kc,i,1)) 
               q2 = 0.5*(aa*qk0(j,i,2,mm)
     .                        +      bb*q(j,kc,i,2)) / xmach
               q3 = 0.5*(aa*qk0(j,i,3,mm)
     .                        +      bb*q(j,kc,i,3)) / xmach
               q4 = 0.5*(aa*qk0(j,i,4,mm)
     .                        +      bb*q(j,kc,i,4)) / xmach
               p1 = 0.5*(aa*qk0(j,i,5,mm)
     .                        +      bb*q(j,kc,i,5)) * gamma
            else
               q1  = 0.5*( q(j,kc,i,1) + q(j,kcm1,i,1))
               q2  = 0.5*( q(j,kc,i,2) + q(j,kcm1,i,2)) / xmach
               q3  = 0.5*( q(j,kc,i,3) + q(j,kcm1,i,3)) / xmach
               q4  = 0.5*( q(j,kc,i,4) + q(j,kcm1,i,4)) / xmach
               p1  = 0.5*( q(j,kc,i,5) + q(j,kcm1,i,5)) * gamma
            end if
            t1  = p1/q1
            xm1 = sqrt(xmach**2*(q2**2+q3**2+q4**2)/t1)
            tt  = t1*(1.0+0.5*gm1*xm1*xm1)
            term1 = (1.0+0.5*gm1*xm1*xm1)**(gamma/gm1)
            pt    = p1*term1
            cmdot = (q2*sk(j,k,i,1) + q3*sk(j,k,i,2) +
     .               q4*sk(j,k,i,3))*q1*sk(j,k,i,4)
            fmdot = fmdot + cmdot
            cfxmom = cfxmom + 2.0*q2*cmdot
            cfymom = cfymom + 2.0*q3*cmdot
            cfzmom = cfzmom + 2.0*q4*cmdot
c
            pav   = pav   + p1*cmdot
            ptav  = ptav  + pt*cmdot
            tav   = tav   + t1*cmdot
            ttav  = ttav  + tt*cmdot
            xmav  = xmav  + xm1*cmdot
         end if
c
c..      skin friction integration
c
         if (iwall.eq.1 .and. ivisc(3).gt.0) then
            urel  = ub(j,kc,i)
            vrel  = vb(j,kc,i)
            wrel  = wb(j,kc,i)
            if (iuns.gt.0) then
               urel  = ub(j,kc,i) - xtbk(j,i,1,kk)
               vrel  = vb(j,kc,i) - xtbk(j,i,2,kk)
               wrel  = wb(j,kc,i) - xtbk(j,i,3,kk)
            end if
            tau   = vmuk(j,i,kk)*const/vol(j,kc,i)*sk(j,k,i,4)**2
            vnorm = (urel*sk(j,k,i,1)+vrel*sk(j,k,i,2)
     .              +wrel*sk(j,k,i,3))
            cfxv  = cfxv + tau*(urel-vnorm*sk(j,k,i,1))
            cfyv  = cfyv + tau*(vrel-vnorm*sk(j,k,i,2))
            cfzv  = cfzv + tau*(wrel-vnorm*sk(j,k,i,3))
         end if
c
   15    continue
c
   10    continue
c
         if (abs(real(fmdot)).lt.0.00001) then
            pav  = 0.e0
            ptav = 0.e0
            tav  = 0.e0
            ttav = 0.e0
            xmav = 0.e0
         else            
            pav   = pav  / fmdot
            ptav  = ptav / fmdot
            tav   = tav  / fmdot
            ttav  = ttav / fmdot
            xmav  = xmav / fmdot
         end if
c
         sx      = sx*fnorm
         sy      = sy*fnorm
         sz      = sz*fnorm
         fmdot   = fmdot*fnorm
         cfxp    = cfxp*fnorm/sref
         cfyp    = cfyp*fnorm/sref
         cfzp    = cfzp*fnorm/sref
         cfxv    = cfxv/sref
         cfyv    = cfyv/sref
         cfzv    = cfzv/sref
         cfxmom  = cfxmom*fnorm/sref
         cfymom  = cfymom*fnorm/sref
         cfzmom  = cfzmom*fnorm/sref
c
         stot  = sqrt(sx*sx+sy*sy+sz*sz)
         cfxtot = cfxp+cfxv+cfxmom
         cfytot = cfyp+cfyv+cfymom
         cfztot = cfzp+cfzv+cfzmom
      end if
c
c******************************************************************
c     forces on j=constant surfaces
c******************************************************************
c
c
      if (jst.eq.jfn+1) then
         j    = jst
         jc   = j
         jcm1 = j - 1
          m    = 0
c
         if (j.eq.1) then
            j    = j
            jc   = 1
            jcm1 = 1
            m    = 1
            mm   = 1
         end if
c
         if (j.eq.jdim) then
            j    = jdim
            jc   = jdim - 1
            jcm1 = jdim - 1
            m    = 2
            mm   = 3
         end if
c
         if (j.eq.1 .and. iwall.eq.1) then
            jj   = 1     
         end if
c
         if (j.eq.jdim .and. iwall.eq.1) then
            jj   = 2   
         end if
c
         do 20 i=ist,ifn
         do 25 k=kst,kfn
c
c..      area integration
c
         sx    = sx + sj(j,k,i,1)*sj(j,k,i,4)
         sy    = sy + sj(j,k,i,2)*sj(j,k,i,4)
         sz    = sz + sj(j,k,i,3)*sj(j,k,i,4)
c
c..      pressure force integration
c
         if (m .ne. 0) then
            aa = 1. +  bcj(k,i,m)
            bb = 1. -  bcj(k,i,m)
            dcp  = -(0.5*(aa*qj0(k,i,5,mm) + bb*q(jc,k,i,5))/p0
     .             -1.e0)*cpc*sj(j,k,i,4)
         else
            dcp  = -(0.5*(q(jc,k,i,5)+q(jcm1,k,i,5))/p0
     .             -1.e0)*cpc*sj(j,k,i,4)
         end if
         cfxp =  cfxp - dcp*sj(j,k,i,1)
         cfyp =  cfyp - dcp*sj(j,k,i,2)
         cfzp =  cfzp - dcp*sj(j,k,i,3)
c
c..      mass and momentum integration
c
         if (iwall.eq.0) then
            if (m .ne. 0) then
               aa = 1. +  bcj(k,i,m)
               bb = 1. -  bcj(k,i,m)
               q1 = 0.5*(aa*qj0(k,i,1,mm)
     .                        +      bb*q(jc,k,i,1))
               q2 = 0.5*(aa*qj0(k,i,2,mm)
     .                        +      bb*q(jc,k,i,2)) / xmach
               q3 = 0.5*(aa*qj0(k,i,3,mm)
     .                        +      bb*q(jc,k,i,3)) / xmach
               q4 = 0.5*(aa*qj0(k,i,4,mm)
     .                        +      bb*q(jc,k,i,4)) / xmach
               p1 = 0.5*(aa*qj0(k,i,5,mm)
     .                        +      bb*q(jc,k,i,5)) * gamma
            else
               q1  = 0.5*( q(jc,k,i,1) + q(jcm1,k,i,1))
               q2  = 0.5*( q(jc,k,i,2) + q(jcm1,k,i,2)) / xmach
               q3  = 0.5*( q(jc,k,i,3) + q(jcm1,k,i,3)) / xmach
               q4  = 0.5*( q(jc,k,i,4) + q(jcm1,k,i,4)) / xmach
               p1  = 0.5*( q(jc,k,i,5) + q(jcm1,k,i,5)) * gamma
            end if
            t1  = p1/q1
c
            xm1 = sqrt(xmach**2*(q2**2+q3**2+q4**2)/t1)
            tt  = t1*(1.0+0.5*gm1*xm1*xm1)
            term1 = (1.0+0.5*gm1*xm1*xm1)**(gamma/gm1)
            pt    = p1*term1
            cmdot = (q2*sj(j,k,i,1) + q3*sj(j,k,i,2) +
     .               q4*sj(j,k,i,3))*q1*sj(j,k,i,4)
            fmdot = fmdot + cmdot
            cfxmom = cfxmom + 2.0*q2*cmdot
            cfymom = cfymom + 2.0*q3*cmdot
            cfzmom = cfzmom + 2.0*q4*cmdot
c
            pav   = pav   + p1*cmdot
            ptav  = ptav  + pt*cmdot
            tav   = tav   + t1*cmdot
            ttav  = ttav  + tt*cmdot
            xmav  = xmav  + xm1*cmdot
         end if
c
c..      skin friction integration
c
         if (iwall.eq.1 .and. ivisc(2).gt.0) then
            urel  = ub(jc,k,i)
            vrel  = vb(jc,k,i)
            wrel  = wb(jc,k,i)
            if (iuns.gt.0) then
               urel  = ub(jc,k,i) - xtbj(k,i,1,jj)
               vrel  = vb(jc,k,i) - xtbj(k,i,2,jj)
               wrel  = wb(jc,k,i) - xtbj(k,i,3,jj)
            end if
            tau   = vmuj(k,i,jj)*const/vol(jc,k,i)*sj(j,k,i,4)**2
            vnorm = (urel*sj(j,k,i,1)+vrel*sj(j,k,i,2)
     .              +wrel*sj(j,k,i,3))
            cfxv  = cfxv + tau*(urel-vnorm*sj(j,k,i,1))
            cfyv  = cfyv + tau*(vrel-vnorm*sj(j,k,i,2))
            cfzv  = cfzv + tau*(wrel-vnorm*sj(j,k,i,3))
         end if
c
   25    continue
c
   20    continue
c
         if (abs(real(fmdot)).lt.0.00001) then
            pav  = 0.e0
            ptav = 0.e0
            tav  = 0.e0
            ttav = 0.e0
            xmav = 0.e0
         else            
            pav   = pav  / fmdot
            ptav  = ptav / fmdot
            tav   = tav  / fmdot
            ttav  = ttav / fmdot
            xmav  = xmav / fmdot
         end if
c
         sx      = sx*fnorm
         sy      = sy*fnorm
         sz      = sz*fnorm
         fmdot   = fmdot*fnorm
         cfxp    = cfxp*fnorm/sref
         cfyp    = cfyp*fnorm/sref
         cfzp    = cfzp*fnorm/sref
         cfxv    = cfxv/sref
         cfyv    = cfyv/sref
         cfzv    = cfzv/sref
         cfxmom  = cfxmom*fnorm/sref
         cfymom  = cfymom*fnorm/sref
         cfzmom  = cfzmom*fnorm/sref
c
         stot  = sqrt(sx*sx+sy*sy+sz*sz)
         cfxtot = cfxp+cfxv+cfxmom
         cfytot = cfyp+cfyv+cfymom
         cfztot = cfzp+cfzv+cfzmom
c
      end if
c
c******************************************************************
c     forces on i=constant surfaces
c******************************************************************
c
      if (ist.eq.ifn+1) then
         i    = ist
         ic   = i
         icm1 = i - 1
         ii   = 1
         m    = 0
c
         if (i.eq.1) then
            ic   = 1
            icm1 = 1
            if (iwall.eq.1) ii = 1
            m    = 1
            mm   = 1
         end if
c
         if (i.eq.idim) then
            ic   = idim - 1
            icm1 = idim - 1
            if (iwall.eq.1) ii = 2
            m    = 2
            mm   = 3
         end if
c
         do 30 j=jst,jfn
         do 35 k=kst,kfn
c
c..      area integration
c
         sx    = sx + si(j,k,i,1)*si(j,k,i,4)
         sy    = sy + si(j,k,i,2)*si(j,k,i,4)
         sz    = sz + si(j,k,i,3)*si(j,k,i,4)
c
c..      pressure force integration
c
         if (m .ne. 0) then
            aa = 1. +  bci(j,k,m)
            bb = 1. -  bci(j,k,m)
            dcp  = -(0.5*(aa*qi0(j,k,5,mm) + bb*q(j,k,ic,5))/p0
     .             -1.e0)*cpc*si(j,k,i,4)
         else
            dcp  = -(0.5*(q(j,k,ic,5)+q(j,k,icm1,5))/p0
     .             -1.e0)*cpc*si(j,k,i,4)
         end if
         cfxp =  cfxp - dcp*si(j,k,i,1)
         cfyp =  cfyp - dcp*si(j,k,i,2)
         cfzp =  cfzp - dcp*si(j,k,i,3)
c
c..      mass and momentum integration
c
         if (iwall.eq.0) then
            if (m .ne. 0) then
               aa = 1. +  bci(j,k,m)
               bb = 1. -  bci(j,k,m)
               q1 = 0.5*(aa*qi0(j,k,1,mm)
     .                        +      bb*q(j,k,ic,1))
               q2 = 0.5*(aa*qi0(j,k,2,mm)
     .                        +      bb*q(j,k,ic,2)) / xmach
               q3 = 0.5*(aa*qi0(j,k,3,mm)
     .                        +      bb*q(j,k,ic,3)) / xmach
               q4 = 0.5*(aa*qi0(j,k,4,mm)
     .                        +      bb*q(j,k,ic,4)) / xmach
               p1 = 0.5*(aa*qi0(j,k,5,mm)
     .                        +      bb*q(j,k,ic,5)) * gamma
            else
               q1  = 0.5*( q(j,k,ic,1) + q(j,k,icm1,1))
               q2  = 0.5*( q(j,k,ic,2) + q(j,k,icm1,2)) / xmach
               q3  = 0.5*( q(j,k,ic,3) + q(j,k,icm1,3)) / xmach
               q4  = 0.5*( q(j,k,ic,4) + q(j,k,icm1,4)) / xmach
               p1  = 0.5*( q(j,k,ic,5) + q(j,k,icm1,5)) * gamma
            end if
            t1  = p1/q1
            xm1 = sqrt(xmach**2*(q2**2+q3**2+q4**2)/t1)
            tt  = t1*(1.0+0.5*gm1*xm1*xm1)
            term1 = (1.0+0.5*gm1*xm1*xm1)**(gamma/gm1)
            pt    = p1*term1
            cmdot = (q2*si(j,k,i,1) + q3*si(j,k,i,2) +
     .               q4*si(j,k,i,3))*q1*si(j,k,i,4)
            fmdot = fmdot + cmdot
            cfxmom = cfxmom + 2.0*q2*cmdot
            cfymom = cfymom + 2.0*q3*cmdot
            cfzmom = cfzmom + 2.0*q4*cmdot
c
            pav   = pav   + p1*cmdot
            ptav  = ptav  + pt*cmdot
            tav   = tav   + t1*cmdot
            ttav  = ttav  + tt*cmdot
            xmav  = xmav  + xm1*cmdot
         end if
c
c..      skin friction integration
c
         if (iwall.eq.1 .and. ivisc(1).gt.0) then
            urel  = ub(j,k,ic)
            vrel  = vb(j,k,ic)
            wrel  = wb(j,k,ic)
            if (iuns.gt.0) then
               urel  = ub(j,k,ic) - xtbi(j,k,1,ii)
               vrel  = vb(j,k,ic) - xtbi(j,k,2,ii)
               wrel  = wb(j,k,ic) - xtbi(j,k,3,ii)
            end if
            tau   = vmui(j,k,ii)*const/vol(j,k,ic)*si(j,k,i,4)**2
            vnorm = (urel*si(j,k,i,1)+vrel*si(j,k,i,2)
     .              +wrel*si(j,k,i,3))
            cfxv  = cfxv + tau*(urel-vnorm*si(j,k,i,1))
            cfyv  = cfyv + tau*(vrel-vnorm*si(j,k,i,2))
            cfzv  = cfzv + tau*(wrel-vnorm*si(j,k,i,3))
         end if
c
   35    continue
c
   30    continue
c
         if (abs(real(fmdot)).lt.0.00001) then
            pav  = 0.e0
            ptav = 0.e0
            tav  = 0.e0
            ttav = 0.e0
            xmav = 0.e0
         else            
            pav   = pav  / fmdot
            ptav  = ptav / fmdot
            tav   = tav  / fmdot
            ttav  = ttav / fmdot
            xmav  = xmav / fmdot
         end if
c
         sx      = sx*fnorm
         sy      = sy*fnorm
         sz      = sz*fnorm
         fmdot   = fmdot*fnorm
         cfxp    = cfxp*fnorm/sref
         cfyp    = cfyp*fnorm/sref
         cfzp    = cfzp*fnorm/sref
         cfxv    = cfxv/sref
         cfyv    = cfyv/sref
         cfzv    = cfzv/sref
         cfxmom  = cfxmom*fnorm/sref
         cfymom  = cfymom*fnorm/sref
         cfzmom  = cfzmom*fnorm/sref
c
         stot  = sqrt(sx*sx+sy*sy+sz*sz)
         cfxtot = cfxp+cfxv+cfxmom
         cfytot = cfyp+cfyv+cfymom
         cfztot = cfzp+cfzv+cfzmom
c
      end if
c
c..   calculate total magnitude and lift and drag components
c
      cftp   = sqrt(cfxp**2+cfyp**2+cfzp**2)
      cftv   = sqrt(cfxv**2+cfyv**2+cfzv**2)
      cftmom = sqrt(cfxmom**2+cfymom**2+cfzmom**2)
      cfttot = sqrt(cfxtot**2+cfytot**2+cfztot**2)
c
      cfdp    = cfxp*cosa*cosb+cfyp*cosa*sinb+cfzp*sina
      cflp    =-cfxp*sina*cosb-cfyp*sina*sinb+cfzp*cosa
      cfdv    = cfxv*cosa*cosb+cfyv*cosa*sinb+cfzv*sina
      cflv    =-cfxv*sina*cosb-cfyv*sina*sinb+cfzv*cosa
      cfdmom  = cfxmom*cosa*cosb+cfymom*cosa*sinb+cfzmom*sina
      cflmom  =-cfxmom*sina*cosb-cfymom*sina*sinb+cfzmom*cosa
      cfdtot  = cfxtot*cosa*cosb+cfytot*cosa*sinb+cfztot*sina
      cfltot  =-cfxtot*sina*cosb-cfytot*sina*sinb+cfztot*cosa
c
      return
      end 
