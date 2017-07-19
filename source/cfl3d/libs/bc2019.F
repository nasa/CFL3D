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
      subroutine bc2019(jdim,kdim,idim,q,qj0,qk0,qi0,sj,sk,si,bcj,bck,
     .                  bci,ista,iend,jsta,jend,ksta,kend,nface,
     .                  tursav,tj0,tk0,ti0,vist3d,vj0,vk0,vi0,
     .                  mdim,ndim,bcdata,filname,iuns,
     .                  nou,bou,nbuf,ibufdim,myid,nummem)
c
c     $Id$
c
c***********************************************************************
c     Purpose: Set inflow boundary conditions (typically for nozzle,
c              duct or engine flows), given total pressure ratio
c              relative to freestream total, and total temperature
c              relative to freestream total.  This boundary condition
c              is taken from that used by OVERFLOW.
c              Also can input up to two turbulence quantities.
c
c     pte   = total to free stream total pressure ratio at 
c             inlet (Ptotal/Ptotal_inf)
c     tte   = total to free stream total temperature ratio at 
c             inlet (Ttotal/Ttotal_inf)
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
      character*80 filname
c
      dimension nou(nbuf)
      dimension q(jdim,kdim,idim,5), qi0(jdim,kdim,5,4),
     .          qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
      dimension sk(jdim,kdim,idim-1,5),si(jdim,kdim,idim,5),
     .          sj(jdim,kdim,idim-1,5)
      dimension bcdata(mdim,ndim,2,12)
      dimension tursav(jdim,kdim,idim,nummem),tj0(kdim,idim-1,nummem,4),
     .          tk0(jdim,idim-1,nummem,4),ti0(jdim,kdim,nummem,4),
     .          vj0(kdim,idim-1,1,4),vk0(jdim,idim-1,1,4),
     .          vi0(jdim,kdim,1,4),vist3d(jdim,kdim,idim)
c
      common /maxiv/ ivmx
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .        tur10(7)
      common /reyue/ reue,tinf,ivisc(3)
      common /sklton/ isklton
      common /conversion/ radtodeg
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
      jend1 = jend-1
      kend1 = kend-1
      iend1 = iend-1
c
c     this bc makes use of only one plane of data   
c
      ip    = 1
c
c            * * * * * * * * * * * * * * * * * * * * * *
c            * standard boundary condition bctype=2019 *
c            * * * * * * * * * * * * * * * * * * * * * *
c
c******************************************************************************
c      j=1 boundary       overflow-type-nozzle total BCs           type 2019
c******************************************************************************
c
      if (nface.eq.3) then
c
c     check to see if turbulence data is input (itrflg1 = 1) or
c     if freestream values are to be used (itrflg1 = 0); the check
c     assumes if the first point has been set, all points have been
c
      itrflg1 = 0
      if (real(bcdata(1,1,ip,3)) .gt. -1.e10) itrflg1 = 1
c
      do 300 i=ista,iend1
      ii = i-ista+1
      js = (i-ista)*(kend-ksta)+1
c
      do 800 k=ksta,kend1
      kk = k-ksta+1
c
c     bcdata(kk,ii,ip,1) is Pt/Ptinf
c     bcdata(kk,ii,ip,2) is Tt/Ttinf
c
c   like overflow (vrgas=1):
c   freestream
      finf=1.+0.5*(gamma-1.)*xmach**2
      p0inf=p0*finf**(gamma/gm1)
      t0inf=p0/rho0*finf
      if (bcdata(kk,ii,ip,1) .ne. 0.) p00=p0inf*bcdata(kk,ii,ip,1)
      if (bcdata(kk,ii,ip,2) .ne. 0.) t00=t0inf*bcdata(kk,ii,ip,2)
      gm1h=0.5*(gamma-1.)
      gdgm1=gamma/gm1
c   external conditions
      rho=qj0(k,i,1,1)
      u  =qj0(k,i,2,1)
      v  =qj0(k,i,3,1)
      w  =qj0(k,i,4,1)
      p  =qj0(k,i,5,1)
      v2 =0.5*(u**2+v**2+w**2)
      ei =p/(rho*gm1)
      c2 =gamma*gm1*ei
      t  =gm1*ei
      fact=1.+gm1*v2/c2
      if (bcdata(kk,ii,ip,2) .eq. 0.) t00=t*fact
      if (bcdata(kk,ii,ip,1) .eq. 0.) p00=p*fact**gdgm1
      t0gr=t00*gamma
c     need outward-pointing (- sign for ilo,jlo,klo):
      uun=-(u*sj(1,k,i,1)+v*sj(1,k,i,2)+w*sj(1,k,i,3)+sj(1,k,i,5))
      un  = -uun*sj(1,k,i,1)
      vn  = -uun*sj(1,k,i,2)
      vw  = -uun*sj(1,k,i,3)
      ut  = u-un
      vt  = v-vn
      wt  = w-wn
      uut2= (ut**2+vt**2+wt**2)
c   internal conditions
      rhox=q(1,k,i,1)
      ux  =q(1,k,i,2)
      vx  =q(1,k,i,3)
      wx  =q(1,k,i,4)
      px  =q(1,k,i,5)
      v2x =0.5*(ux**2+vx**2+wx**2)
      eix =px/(rhox*gm1)
      c2x =gamma*gm1*eix
      cx  =sqrt(c2x)
c     need outward-pointing (- sign for ilo,jlo,klo):
      uunx=-(ux*sj(1,k,i,1)+vx*sj(1,k,i,2)+wx*sj(1,k,i,3)+sj(1,k,i,5))
      r2  =uunx + cx/gm1h
c   flow angle
      if (v2 .eq. 0.) then
        tan2p1 = 1.
      else if (uun .eq. 0.) then
        tan2p1=1.e20
      else
        tan2p1 = uut2/uun**2+1.
        tan2p1 = ccmincr(tan2p1,1.e20)
      end if
c   find speed of sound
      aterm = 1.+tan2p1/gm1h
      hbtrmm= tan2p1*r2
      cterm = hbtrmm*gm1h*r2 - t0gr
      sqterm= hbtrmm**2 - aterm*cterm
      c     = (hbtrmm + sqrt(ccmaxcr(sqterm,0.)))/aterm
      c2    = c**2
c   final values for setting BC
      uunold = uun
      uun    = -abs(r2 - c/gm1h)
      v2old  = v2
      v2     = 0.5*uun**2*tan2p1
      fact   = 1. + gm1*v2/c2
      p      = p00*fact**(-gdgm1)
      rho    = gamma*p/c2
      if (v2old .eq. 0.) then
        vratio = uun
        u = 1.
        v = 0.
        w = 0.
      else if (uunold .eq. 0. .or. tan2p1 .eq. 1.e20) then
        vratio = sqrt(v2/v2old)
      else
        vratio = uun/uunold
      end if
c
      qj0(k,i,1,1) = rho
      qj0(k,i,1,2) = qj0(k,i,1,1)
      qj0(k,i,2,1) = u*vratio
      qj0(k,i,2,2) = qj0(k,i,2,1)
      qj0(k,i,3,1) = v*vratio
      qj0(k,i,3,2) = qj0(k,i,3,1)
      qj0(k,i,4,1) = w*vratio
      qj0(k,i,4,2) = qj0(k,i,4,1)
      qj0(k,i,5,1) = p
      qj0(k,i,5,2) = qj0(k,i,5,1)
c
      bcj(k,i,1)   = 0.0
c
  800 continue
  300 continue
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 191 i=ista,iend1
        do 191 k=ksta,kend1
          vj0(k,i,1,1) = vist3d(1,k,i)
          vj0(k,i,1,2) = 0.0
  191   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 101 i=ista,iend1
        ii = i-ista+1
        do 101 k=ksta,kend1
          kk=k-ksta+1
          ubar=-(qj0(k,i,2,1)*sj(1,k,i,1)+qj0(k,i,3,1)*sj(1,k,i,2)+
     +           qj0(k,i,4,1)*sj(1,k,i,3))
          t1 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(kk,ii,ip,2+l)
          if (real(ubar) .lt. 0.) then
             tj0(k,i,l,1) = t1
             tj0(k,i,l,2) = t1
          else
             tj0(k,i,l,1) = tursav(1,k,i,l)
             tj0(k,i,l,2) = tj0(k,i,l,1)
          end if
  101   continue
        enddo
      end if
      end if
c
      end if
c
c******************************************************************************
c      j=jdim boundary       overflow-type-nozzle total BCs           type 2019
c******************************************************************************
c
      if (nface.eq.4) then
c
c     check to see if turbulence data is input (itrflg1 = 1) or
c     if freestream values are to be used (itrflg1 = 0); the check
c     assumes if the first point has been set, all points have been
c
      itrflg1 = 0
      if (real(bcdata(1,1,ip,3)) .gt. -1.e10) itrflg1 = 1
c
      do 310 i=ista,iend1
      ii = i-ista+1
      js = (i-ista)*(kend-ksta)+1
c
      do 810 k=ksta,kend1
      kk = k-ksta+1
c
c     bcdata(kk,ii,ip,1) is Pt/Ptinf
c     bcdata(kk,ii,ip,2) is Tt/Ttinf
c
c   like overflow (vrgas=1):
c   freestream
      finf=1.+0.5*(gamma-1.)*xmach**2
      p0inf=p0*finf**(gamma/gm1)
      t0inf=p0/rho0*finf
      if (bcdata(kk,ii,ip,1) .ne. 0.) p00=p0inf*bcdata(kk,ii,ip,1)
      if (bcdata(kk,ii,ip,2) .ne. 0.) t00=t0inf*bcdata(kk,ii,ip,2)
      gm1h=0.5*(gamma-1.)
      gdgm1=gamma/gm1
c   external conditions
      rho=qj0(k,i,1,3)
      u  =qj0(k,i,2,3)
      v  =qj0(k,i,3,3)
      w  =qj0(k,i,4,3)
      p  =qj0(k,i,5,3)
      v2 =0.5*(u**2+v**2+w**2)
      ei =p/(rho*gm1)
      c2 =gamma*gm1*ei
      t  =gm1*ei
      fact=1.+gm1*v2/c2
      if (bcdata(kk,ii,ip,2) .eq. 0.) t00=t*fact
      if (bcdata(kk,ii,ip,1) .eq. 0.) p00=p*fact**gdgm1
      t0gr=t00*gamma
c     need outward-pointing (- sign for ilo,jlo,klo):
      uun=+(u*sj(jdim,k,i,1)+v*sj(jdim,k,i,2)+w*sj(jdim,k,i,3)+
     +     sj(jdim,k,i,5))
      un  = +uun*sj(jdim,k,i,1)
      vn  = +uun*sj(jdim,k,i,2)
      vw  = +uun*sj(jdim,k,i,3)
      ut  = u-un
      vt  = v-vn
      wt  = w-wn
      uut2= (ut**2+vt**2+wt**2)
c   internal conditions
      rhox=q(jdim1,k,i,1)
      ux  =q(jdim1,k,i,2)
      vx  =q(jdim1,k,i,3)
      wx  =q(jdim1,k,i,4)
      px  =q(jdim1,k,i,5)
      v2x =0.5*(ux**2+vx**2+wx**2)
      eix =px/(rhox*gm1)
      c2x =gamma*gm1*eix
      cx  =sqrt(c2x)
c     need outward-pointing (- sign for ilo,jlo,klo):
      uunx=+(ux*sj(jdim,k,i,1)+vx*sj(jdim,k,i,2)+wx*sj(jdim,k,i,3)+
     +      sj(jdim,k,i,5))
      r2  =uunx + cx/gm1h
c   flow angle
      if (v2 .eq. 0.) then
        tan2p1 = 1.
      else if (uun .eq. 0.) then
        tan2p1=1.e20
      else
        tan2p1 = uut2/uun**2+1.
        tan2p1 = ccmincr(tan2p1,1.e20)
      end if
c   find speed of sound
      aterm = 1.+tan2p1/gm1h
      hbtrmm= tan2p1*r2
      cterm = hbtrmm*gm1h*r2 - t0gr
      sqterm= hbtrmm**2 - aterm*cterm
      c     = (hbtrmm + sqrt(ccmaxcr(sqterm,0.)))/aterm
      c2    = c**2
c   final values for setting BC
      uunold = uun
      uun    = -abs(r2 - c/gm1h)
      v2old  = v2
      v2     = 0.5*uun**2*tan2p1
      fact   = 1. + gm1*v2/c2
      p      = p00*fact**(-gdgm1)
      rho    = gamma*p/c2
      if (v2old .eq. 0.) then
        vratio = uun
        u = 1.
        v = 0.
        w = 0.
      else if (uunold .eq. 0. .or. tan2p1 .eq. 1.e20) then
        vratio = sqrt(v2/v2old)
      else
        vratio = uun/uunold
      end if
c
      qj0(k,i,1,3) = rho
      qj0(k,i,1,4) = qj0(k,i,1,3)
      qj0(k,i,2,3) = u*vratio
      qj0(k,i,2,4) = qj0(k,i,2,3)
      qj0(k,i,3,3) = v*vratio
      qj0(k,i,3,4) = qj0(k,i,3,3)
      qj0(k,i,4,3) = w*vratio
      qj0(k,i,4,4) = qj0(k,i,4,3)
      qj0(k,i,5,3) = p
      qj0(k,i,5,4) = qj0(k,i,5,3)
c
      bcj(k,i,2)   = 0.0
c
  810 continue
  310 continue
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 291 i=ista,iend1
        do 291 k=ksta,kend1
          vj0(k,i,1,3) = vist3d(jdim-1,k,i)
          vj0(k,i,1,4) = 0.0
  291   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 201 i=ista,iend1
        ii = i-ista+1
        do 201 k=ksta,kend1
          kk = k-ksta+1
          ubar=qj0(k,i,2,3)*sj(jdim,k,i,1)+qj0(k,i,3,3)*sj(jdim,k,i,2)+
     +         qj0(k,i,4,3)*sj(jdim,k,i,3)
          t1 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(kk,ii,ip,2+l)
          if (real(ubar) .lt. 0.) then
             tj0(k,i,l,3) = t1
             tj0(k,i,l,4) = t1
          else
             tj0(k,i,l,3) = tursav(jdim-1,k,i,l)
             tj0(k,i,l,4) = tj0(k,i,l,3)
          end if
  201   continue
        enddo
      end if
      end if
c
      end if
c
c******************************************************************************
c      k=1 boundary       overflow-type-nozzle total BCs           type 2019
c******************************************************************************
c
      if (nface.eq.5) then
c
c     check to see if turbulence data is input (itrflg1 = 1) or
c     if freestream values are to be used (itrflg1 = 0); the check
c     assumes if the first point has been set, all points have been
c
      itrflg1 = 0
      if (real(bcdata(1,1,ip,3)) .gt. -1.e10) itrflg1 = 1
c
      do 320 i=ista,iend1
      ii = i-ista+1
      js = (i-ista)*(jend-jsta)+1
c
      do 820 j=jsta,jend1
      jj = j-jsta+1
c
c     bcdata(jj,ii,ip,1) is Pt/Ptinf
c     bcdata(jj,ii,ip,2) is Tt/Ttinf
c
c   like overflow (vrgas=1):
c   freestream
      finf=1.+0.5*(gamma-1.)*xmach**2
      p0inf=p0*finf**(gamma/gm1)
      t0inf=p0/rho0*finf
      if (bcdata(jj,ii,ip,1) .ne. 0.) p00=p0inf*bcdata(jj,ii,ip,1)
      if (bcdata(jj,ii,ip,2) .ne. 0.) t00=t0inf*bcdata(jj,ii,ip,2)
      gm1h=0.5*(gamma-1.)
      gdgm1=gamma/gm1
c   external conditions
      rho=qk0(j,i,1,1)
      u  =qk0(j,i,2,1)
      v  =qk0(j,i,3,1)
      w  =qk0(j,i,4,1)
      p  =qk0(j,i,5,1)
      v2 =0.5*(u**2+v**2+w**2)
      ei =p/(rho*gm1)
      c2 =gamma*gm1*ei
      t  =gm1*ei
      fact=1.+gm1*v2/c2
      if (bcdata(jj,ii,ip,2) .eq. 0.) t00=t*fact
      if (bcdata(jj,ii,ip,1) .eq. 0.) p00=p*fact**gdgm1
      t0gr=t00*gamma
c     need outward-pointing (- sign for ilo,jlo,klo):
      uun=-(u*sk(j,1,i,1)+v*sk(j,1,i,2)+w*sk(j,1,i,3)+sk(j,1,i,5))
      un  = -uun*sk(j,1,i,1)
      vn  = -uun*sk(j,1,i,2)
      vw  = -uun*sk(j,1,i,3)
      ut  = u-un
      vt  = v-vn
      wt  = w-wn
      uut2= (ut**2+vt**2+wt**2)
c   internal conditions
      rhox=q(j,1,i,1)
      ux  =q(j,1,i,2)
      vx  =q(j,1,i,3)
      wx  =q(j,1,i,4)
      px  =q(j,1,i,5)
      v2x =0.5*(ux**2+vx**2+wx**2)
      eix =px/(rhox*gm1)
      c2x =gamma*gm1*eix
      cx  =sqrt(c2x)
c     need outward-pointing (- sign for ilo,jlo,klo):
      uunx=-(ux*sk(j,1,i,1)+vx*sk(j,1,i,2)+wx*sk(j,1,i,3)+sk(j,1,i,5))
      r2  =uunx + cx/gm1h
c   flow angle
      if (v2 .eq. 0.) then
        tan2p1 = 1.
      else if (uun .eq. 0.) then
        tan2p1=1.e20
      else
        tan2p1 = uut2/uun**2+1.
        tan2p1 = ccmincr(tan2p1,1.e20)
      end if
c   find speed of sound
      aterm = 1.+tan2p1/gm1h
      hbtrmm= tan2p1*r2
      cterm = hbtrmm*gm1h*r2 - t0gr
      sqterm= hbtrmm**2 - aterm*cterm
      c     = (hbtrmm + sqrt(ccmaxcr(sqterm,0.)))/aterm
      c2    = c**2
c   final values for setting BC
      uunold = uun
      uun    = -abs(r2 - c/gm1h)
      v2old  = v2
      v2     = 0.5*uun**2*tan2p1
      fact   = 1. + gm1*v2/c2
      p      = p00*fact**(-gdgm1)
      rho    = gamma*p/c2
      if (v2old .eq. 0.) then
        vratio = uun
        u = 1.
        v = 0.
        w = 0.
      else if (uunold .eq. 0. .or. tan2p1 .eq. 1.e20) then
        vratio = sqrt(v2/v2old)
      else
        vratio = uun/uunold
      end if
c
      qk0(j,i,1,1) = rho
      qk0(j,i,1,2) = qk0(j,i,1,1)
      qk0(j,i,2,1) = u*vratio
      qk0(j,i,2,2) = qk0(j,i,2,1)
      qk0(j,i,3,1) = v*vratio
      qk0(j,i,3,2) = qk0(j,i,3,1)
      qk0(j,i,4,1) = w*vratio
      qk0(j,i,4,2) = qk0(j,i,4,1)
      qk0(j,i,5,1) = p
      qk0(j,i,5,2) = qk0(j,i,5,1)
c
      bck(j,i,1)   = 0.0
c
  820 continue
  320 continue
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 391 i=ista,iend1
        do 391 j=jsta,jend1
          vk0(j,i,1,1) = vist3d(j,1,i)
          vk0(j,i,1,2) = 0.0
  391   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 301 i=ista,iend1
        ii = i-ista+1
        do 301 j=jsta,jend1
          jj = j-jsta+1
          ubar=-(qk0(j,i,2,1)*sk(j,1,i,1)+qk0(j,i,3,1)*sk(j,1,i,2)+
     +           qk0(j,i,4,1)*sk(j,1,i,3))
          t1 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(jj,ii,ip,2+l)
          if (real(ubar) .lt. 0.) then
             tk0(j,i,l,1) = t1
             tk0(j,i,l,2) = t1
          else
             tk0(j,i,l,1) = tursav(j,1,i,l)
             tk0(j,i,l,2) = tk0(j,i,l,1)
          end if
  301   continue
        enddo
      end if
      end if
c
      end if
c
c******************************************************************************
c      k=kdim boundary       overflow-type-nozzle total BCs           type 2019
c******************************************************************************
c
      if (nface.eq.6) then
c
c     check to see if turbulence data is input (itrflg1 = 1) or
c     if freestream values are to be used (itrflg1 = 0); the check
c     assumes if the first point has been set, all points have been
c
      itrflg1 = 0
      if (real(bcdata(1,1,ip,3)) .gt. -1.e10) itrflg1 = 1
c
      do 330 i=ista,iend1
      ii = i-ista+1
      js = (i-ista)*(jend-jsta)+1
c
      do 830 j=jsta,jend1
      jj = j-jsta+1
c
c     bcdata(jj,ii,ip,1) is Pt/Ptinf
c     bcdata(jj,ii,ip,2) is Tt/Ttinf
c
c   like overflow (vrgas=1):
c   freestream
      finf=1.+0.5*(gamma-1.)*xmach**2
      p0inf=p0*finf**(gamma/gm1)
      t0inf=p0/rho0*finf
      if (bcdata(jj,ii,ip,1) .ne. 0.) p00=p0inf*bcdata(jj,ii,ip,1)
      if (bcdata(jj,ii,ip,2) .ne. 0.) t00=t0inf*bcdata(jj,ii,ip,2)
      gm1h=0.5*(gamma-1.)
      gdgm1=gamma/gm1
c   external conditions
      rho=qk0(j,i,1,3)
      u  =qk0(j,i,2,3)
      v  =qk0(j,i,3,3)
      w  =qk0(j,i,4,3)
      p  =qk0(j,i,5,3)
      v2 =0.5*(u**2+v**2+w**2)
      ei =p/(rho*gm1)
      c2 =gamma*gm1*ei
      t  =gm1*ei
      fact=1.+gm1*v2/c2
      if (bcdata(jj,ii,ip,2) .eq. 0.) t00=t*fact
      if (bcdata(jj,ii,ip,1) .eq. 0.) p00=p*fact**gdgm1
      t0gr=t00*gamma
c     need outward-pointing (- sign for ilo,jlo,klo):
      uun=+(u*sk(j,kdim,i,1)+v*sk(j,kdim,i,2)+w*sk(j,kdim,i,3)+
     +     sk(j,kdim,i,5))
      un  = +uun*sk(j,kdim,i,1)
      vn  = +uun*sk(j,kdim,i,2)
      vw  = +uun*sk(j,kdim,i,3)
      ut  = u-un
      vt  = v-vn
      wt  = w-wn
      uut2= (ut**2+vt**2+wt**2)
c   internal conditions
      rhox=q(j,kdim1,i,1)
      ux  =q(j,kdim1,i,2)
      vx  =q(j,kdim1,i,3)
      wx  =q(j,kdim1,i,4)
      px  =q(j,kdim1,i,5)
      v2x =0.5*(ux**2+vx**2+wx**2)
      eix =px/(rhox*gm1)
      c2x =gamma*gm1*eix
      cx  =sqrt(c2x)
c     need outward-pointing (- sign for ilo,jlo,klo):
      uunx=+(ux*sk(j,kdim,i,1)+vx*sk(j,kdim,i,2)+wx*sk(j,kdim,i,3)+
     +      sk(j,kdim,i,5))
      r2  =uunx + cx/gm1h
c   flow angle
      if (v2 .eq. 0.) then
        tan2p1 = 1.
      else if (uun .eq. 0.) then
        tan2p1=1.e20
      else
        tan2p1 = uut2/uun**2+1.
        tan2p1 = ccmincr(tan2p1,1.e20)
      end if
c   find speed of sound
      aterm = 1.+tan2p1/gm1h
      hbtrmm= tan2p1*r2
      cterm = hbtrmm*gm1h*r2 - t0gr
      sqterm= hbtrmm**2 - aterm*cterm
      c     = (hbtrmm + sqrt(ccmaxcr(sqterm,0.)))/aterm
      c2    = c**2
c   final values for setting BC
      uunold = uun
      uun    = -abs(r2 - c/gm1h)
      v2old  = v2
      v2     = 0.5*uun**2*tan2p1
      fact   = 1. + gm1*v2/c2
      p      = p00*fact**(-gdgm1)
      rho    = gamma*p/c2
      if (v2old .eq. 0.) then
        vratio = uun
        u = 1.
        v = 0.
        w = 0.
      else if (uunold .eq. 0. .or. tan2p1 .eq. 1.e20) then
        vratio = sqrt(v2/v2old)
      else
        vratio = uun/uunold
      end if
c
      qk0(j,i,1,3) = rho
      qk0(j,i,1,4) = qk0(j,i,1,3)
      qk0(j,i,2,3) = u*vratio
      qk0(j,i,2,4) = qk0(j,i,2,3)
      qk0(j,i,3,3) = v*vratio
      qk0(j,i,3,4) = qk0(j,i,3,3)
      qk0(j,i,4,3) = w*vratio
      qk0(j,i,4,4) = qk0(j,i,4,3)
      qk0(j,i,5,3) = p
      qk0(j,i,5,4) = qk0(j,i,5,3)
c
      bck(j,i,2)   = 0.0
c
  830 continue
  330 continue
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 491 i=ista,iend1
        do 491 j=jsta,jend1
          vk0(j,i,1,3) = vist3d(j,kdim-1,i)
          vk0(j,i,1,4) = 0.0
  491   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 401 i=ista,iend1
        ii = i-ista+1
        do 401 j=jsta,jend1
          jj = j-jsta+1
          ubar=qk0(j,i,2,3)*sk(j,kdim,i,1)+qk0(j,i,3,3)*sk(j,kdim,i,2)+
     +         qk0(j,i,4,3)*sk(j,kdim,i,3)
          t1 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(jj,ii,ip,2+l)
          if (real(ubar) .lt. 0.) then
             tk0(j,i,l,3) = t1
             tk0(j,i,l,4) = t1
          else
             tk0(j,i,l,3) = tursav(j,kdim-1,i,l)
             tk0(j,i,l,4) = tk0(j,i,l,3)
          end if
  401   continue
        enddo
      end if
      end if
c
      end if
c
c******************************************************************************
c      i=1 boundary       overflow-type-nozzle total BCs           type 2019
c******************************************************************************
c
      if (nface.eq.1) then
c
c     check to see if turbulence data is input (itrflg1 = 1) or
c     if freestream values are to be used (itrflg1 = 0); the check
c     assumes if the first point has been set, all points have been
c
      itrflg1 = 0
      if (real(bcdata(1,1,ip,3)) .gt. -1.e10) itrflg1 = 1
c
      do 340 k=ksta,kend1
      kk = k-ksta+1
      js = (k-ksta)*(jend-jsta)+1
c
      do 840 j=jsta,jend1
      jj = j-jsta+1
c
c     bcdata(jj,kk,ip,1) is Pt/Ptinf
c     bcdata(jj,kk,ip,2) is Tt/Ttinf
c
c   like overflow (vrgas=1):
c   freestream
      finf=1.+0.5*(gamma-1.)*xmach**2
      p0inf=p0*finf**(gamma/gm1)
      t0inf=p0/rho0*finf
      if (bcdata(jj,kk,ip,1) .ne. 0.) p00=p0inf*bcdata(jj,kk,ip,1)
      if (bcdata(jj,kk,ip,2) .ne. 0.) t00=t0inf*bcdata(jj,kk,ip,2)
      gm1h=0.5*(gamma-1.)
      gdgm1=gamma/gm1
c   external conditions
      rho=qi0(j,k,1,1)
      u  =qi0(j,k,2,1)
      v  =qi0(j,k,3,1)
      w  =qi0(j,k,4,1)
      p  =qi0(j,k,5,1)
      v2 =0.5*(u**2+v**2+w**2)
      ei =p/(rho*gm1)
      c2 =gamma*gm1*ei
      t  =gm1*ei
      fact=1.+gm1*v2/c2
      if (bcdata(jj,kk,ip,2) .eq. 0.) t00=t*fact
      if (bcdata(jj,kk,ip,1) .eq. 0.) p00=p*fact**gdgm1
      t0gr=t00*gamma
c     need outward-pointing (- sign for ilo,jlo,klo):
      uun=-(u*si(j,k,1,1)+v*si(j,k,1,2)+w*si(j,k,1,3)+si(j,k,1,5))
      un  = -uun*si(j,k,1,1)
      vn  = -uun*si(j,k,1,2)
      vw  = -uun*si(j,k,1,3)
      ut  = u-un
      vt  = v-vn
      wt  = w-wn
      uut2= (ut**2+vt**2+wt**2)
c   internal conditions
      rhox=q(j,k,1,1)
      ux  =q(j,k,1,2)
      vx  =q(j,k,1,3)
      wx  =q(j,k,1,4)
      px  =q(j,k,1,5)
      v2x =0.5*(ux**2+vx**2+wx**2)
      eix =px/(rhox*gm1)
      c2x =gamma*gm1*eix
      cx  =sqrt(c2x)
c     need outward-pointing (- sign for ilo,jlo,klo):
      uunx=-(ux*si(j,k,1,1)+vx*si(j,k,1,2)+wx*si(j,k,1,3)+si(j,k,1,5))
      r2  =uunx + cx/gm1h
c   flow angle
      if (v2 .eq. 0.) then
        tan2p1 = 1.
      else if (uun .eq. 0.) then
        tan2p1=1.e20
      else
        tan2p1 = uut2/uun**2+1.
        tan2p1 = ccmincr(tan2p1,1.e20)
      end if
c   find speed of sound
      aterm = 1.+tan2p1/gm1h
      hbtrmm= tan2p1*r2
      cterm = hbtrmm*gm1h*r2 - t0gr
      sqterm= hbtrmm**2 - aterm*cterm
      c     = (hbtrmm + sqrt(ccmaxcr(sqterm,0.)))/aterm
      c2    = c**2
c   final values for setting BC
      uunold = uun
      uun    = -abs(r2 - c/gm1h)
      v2old  = v2
      v2     = 0.5*uun**2*tan2p1
      fact   = 1. + gm1*v2/c2
      p      = p00*fact**(-gdgm1)
      rho    = gamma*p/c2
      if (v2old .eq. 0.) then
        vratio = uun
        u = 1.
        v = 0.
        w = 0.
      else if (uunold .eq. 0. .or. tan2p1 .eq. 1.e20) then
        vratio = sqrt(v2/v2old)
      else
        vratio = uun/uunold
      end if
c
      qi0(j,k,1,1) = rho
      qi0(j,k,1,2) = qi0(j,k,1,1)
      qi0(j,k,2,1) = u*vratio
      qi0(j,k,2,2) = qi0(j,k,2,1)
      qi0(j,k,3,1) = v*vratio
      qi0(j,k,3,2) = qi0(j,k,3,1)
      qi0(j,k,4,1) = w*vratio
      qi0(j,k,4,2) = qi0(j,k,4,1)
      qi0(j,k,5,1) = p
      qi0(j,k,5,2) = qi0(j,k,5,1)
c
      bci(j,k,1)   = 0.0
c
  840 continue
  340 continue
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 591 k=ksta,kend1
        do 591 j=jsta,jend1
          vi0(j,k,1,1) = vist3d(j,k,1)
          vi0(j,k,1,2) = 0.0
  591   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 501 k=ksta,kend1
        kk = k-ksta+1
        do 501 j=jsta,jend1
          jj = j-jsta+1
          ubar=-(qi0(j,k,2,1)*si(j,k,1,1)+qi0(j,k,3,1)*si(j,k,1,2)+
     +           qi0(j,k,4,1)*si(j,k,1,3))
          t1 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(jj,kk,ip,2+l)
          if (real(ubar) .lt. 0.) then
             ti0(j,k,l,1) = t1
             ti0(j,k,l,2) = t1
          else
             ti0(j,k,l,1) = tursav(j,k,1,l)
             ti0(j,k,l,2) = ti0(j,k,l,1)
          end if
  501   continue
        enddo
      end if
      end if
c
      end if
c
c******************************************************************************
c      i=idim boundary       overflow-type-nozzle total BCs           type 2019
c******************************************************************************
c
      if (nface.eq.2) then
c
c     check to see if turbulence data is input (itrflg1 = 1) or
c     if freestream values are to be used (itrflg1 = 0); the check
c     assumes if the first point has been set, all points have been
c
      itrflg1 = 0
      if (real(bcdata(1,1,ip,3)) .gt. -1.e10) itrflg1 = 1
c
      do 350 k=ksta,kend1
      kk = k-ksta+1
      js = (k-ksta)*(jend-jsta)+1
c
      do 850 j=jsta,jend1
      jj = j-jsta+1
c
c     bcdata(jj,kk,ip,1) is Pt/Ptinf
c     bcdata(jj,kk,ip,2) is Tt/Ttinf
c
c   like overflow (vrgas=1):
c   freestream
      finf=1.+0.5*(gamma-1.)*xmach**2
      p0inf=p0*finf**(gamma/gm1)
      t0inf=p0/rho0*finf
      if (bcdata(jj,kk,ip,1) .ne. 0.) p00=p0inf*bcdata(jj,kk,ip,1)
      if (bcdata(jj,kk,ip,2) .ne. 0.) t00=t0inf*bcdata(jj,kk,ip,2)
      gm1h=0.5*(gamma-1.)
      gdgm1=gamma/gm1
c   external conditions
      rho=qi0(j,k,1,3)
      u  =qi0(j,k,2,3)
      v  =qi0(j,k,3,3)
      w  =qi0(j,k,4,3)
      p  =qi0(j,k,5,3)
      v2 =0.5*(u**2+v**2+w**2)
      ei =p/(rho*gm1)
      c2 =gamma*gm1*ei
      t  =gm1*ei
      fact=1.+gm1*v2/c2
      if (bcdata(jj,kk,ip,2) .eq. 0.) t00=t*fact
      if (bcdata(jj,kk,ip,1) .eq. 0.) p00=p*fact**gdgm1
      t0gr=t00*gamma
c     need outward-pointing (- sign for ilo,jlo,klo):
      uun=+(u*si(j,k,idim,1)+v*si(j,k,idim,2)+w*si(j,k,idim,3)+
     +     si(j,k,idim,5))
      un  = +uun*si(j,k,idim,1)
      vn  = +uun*si(j,k,idim,2)
      vw  = +uun*si(j,k,idim,3)
      ut  = u-un
      vt  = v-vn
      wt  = w-wn
      uut2= (ut**2+vt**2+wt**2)
c   internal conditions
      rhox=q(j,k,idim1,1)
      ux  =q(j,k,idim1,2)
      vx  =q(j,k,idim1,3)
      wx  =q(j,k,idim1,4)
      px  =q(j,k,idim1,5)
      v2x =0.5*(ux**2+vx**2+wx**2)
      eix =px/(rhox*gm1)
      c2x =gamma*gm1*eix
      cx  =sqrt(c2x)
c     need outward-pointing (- sign for ilo,jlo,klo):
      uunx=+(ux*si(j,k,idim,1)+vx*si(j,k,idim,2)+wx*si(j,k,idim,3)+
     +      si(j,k,idim,5))
      r2  =uunx + cx/gm1h
c   flow angle
      if (v2 .eq. 0.) then
        tan2p1 = 1.
      else if (uun .eq. 0.) then
        tan2p1=1.e20
      else
        tan2p1 = uut2/uun**2+1.
        tan2p1 = ccmincr(tan2p1,1.e20)
      end if
c   find speed of sound
      aterm = 1.+tan2p1/gm1h
      hbtrmm= tan2p1*r2
      cterm = hbtrmm*gm1h*r2 - t0gr
      sqterm= hbtrmm**2 - aterm*cterm
      c     = (hbtrmm + sqrt(ccmaxcr(sqterm,0.)))/aterm
      c2    = c**2
c   final values for setting BC
      uunold = uun
      uun    = -abs(r2 - c/gm1h)
      v2old  = v2
      v2     = 0.5*uun**2*tan2p1
      fact   = 1. + gm1*v2/c2
      p      = p00*fact**(-gdgm1)
      rho    = gamma*p/c2
      if (v2old .eq. 0.) then
        vratio = uun
        u = 1.
        v = 0.
        w = 0.
      else if (uunold .eq. 0. .or. tan2p1 .eq. 1.e20) then
        vratio = sqrt(v2/v2old)
      else
        vratio = uun/uunold
      end if
c
      qi0(j,k,1,3) = rho
      qi0(j,k,1,4) = qi0(j,k,1,3)
      qi0(j,k,2,3) = u*vratio
      qi0(j,k,2,4) = qi0(j,k,2,3)
      qi0(j,k,3,3) = v*vratio
      qi0(j,k,3,4) = qi0(j,k,3,3)
      qi0(j,k,4,3) = w*vratio
      qi0(j,k,4,4) = qi0(j,k,4,3)
      qi0(j,k,5,3) = p
      qi0(j,k,5,4) = qi0(j,k,5,3)
c
      bci(j,k,2)   = 0.0
c
  850 continue
  350 continue
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 691 k=ksta,kend1
        do 691 j=jsta,jend1
          vi0(j,k,1,3) = vist3d(j,k,idim-1)
          vi0(j,k,1,4) = 0.0
  691   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 601 k=ksta,kend1
        kk = k-ksta+1
        do 601 j=jsta,jend1
          jj = j-jsta+1
          ubar=qi0(j,k,2,3)*si(j,k,idim,1)+qi0(j,k,3,3)*si(j,k,idim,2)+
     +         qi0(j,k,4,3)*si(j,k,idim,3)
          t1 = (1 - itrflg1)*tur10(l) + itrflg1*bcdata(jj,kk,ip,2+l)
          if (real(ubar) .lt. 0.) then
             ti0(j,k,l,3) = t1
             ti0(j,k,l,4) = t1
          else
             ti0(j,k,l,3) = tursav(j,k,idim-1,l)
             ti0(j,k,l,4) = ti0(j,k,l,3)
          end if
  601   continue
        enddo
      end if
      end if
c
      end if
c
      return
      end
