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
      subroutine force(jdim,kdim,idim,x,y,z,sk,sj,si,cl,cd,cz,cy,cx,
     .                 cmy,cmx,cmz,chd,swet,i00,ub,vb,wb,vmuk,vmuj,vmui,
     .                 vol,ifo,jfo,kfo,bcj,bck,bci,blank,nbl,
     .                 xtbj,xtbk,xtbi,iuns,qj0,qk0,qi0,nbci0,nbcj0,
     .                 nbck0,nbcidim,nbcjdim,nbckdim,ibcinfo,jbcinfo,
     .                 kbcinfo,nn,maxseg)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Integrate the forces on the body.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension ub(jdim,kdim,idim),vb(jdim,kdim,idim),wb(jdim,kdim,idim)
      dimension vmuk(jdim-1,idim-1,2),vol(jdim,kdim,idim-1),
     .          vmuj(kdim-1,idim-1,2),vmui(jdim-1,kdim-1,2)
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
      dimension xtbj(kdim,idim-1,3,2),xtbk(jdim,idim-1,3,2),
     .          xtbi(jdim,kdim,3,2)

      dimension sk(jdim,kdim,idim-1,5),sj(jdim,kdim,idim-1,5),
     .          si(jdim,kdim,idim,5)
      dimension blank(jdim,kdim,idim)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
      dimension qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4),
     .          qi0(jdim,kdim,5,4)
      dimension nbci0(nn),nbcidim(nn),nbcj0(nn),nbcjdim(nn),
     .          nbck0(nn),nbckdim(nn),ibcinfo(nn,maxseg,7,2),
     .          jbcinfo(nn,maxseg,7,2),kbcinfo(nn,maxseg,7,2)
c
      common /drag/ cdv,cdp
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /fsum/ sref,cref,bref,xmc,ymc,zmc
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /twod/ i2d
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .        tur10(7)
      common /reyue/ reue,tinf,ivisc(3)
      common /igrdtyp/ ip3dgrd,ialph
c
c**************************************************************************
c
c     for force computation on block nbl:
c
c     kbcinfo(nbl,nseg,6,1) > 0...add contribution from k=1,    segment nseg 
c                           <=0...do not add contribution
c
c     kbcinfo(nbl,nseg,6,2) > 0...add contribution from k=kdim, segment nseg 
c                           <=0...do not add contribution
c
c     jbcinfo(nbl,nseg,6,1) > 0...add contribution from j=1,    segment nseg
c                           <=0...do not add contribution
c
c     jbcinfo(nbl,nseg,6,2) > 0...add contribution from j=jdim, segment nseg
c                           <=0...do not add contribution
c
c     ibcinfo(nbl,nseg,6,1) > 0...add contribution from i=1,    segment nseg
c                           <=0...do not add contribution
c
c     ibcinfo(nbl,nseg,6,2) > 0...add contribution from i=idim, segment nseg
c                           <=0...do not add contribution
c
c     sign conventions for forces and moments:
c       cmx....moment about x-axis through the moment center (xmc,ymc,zmc)
c              positive value for ccw moment when viewed from +x axis
c       cmy....moment about y-axis through the moment center (xmc,ymc,zmc)
c              positive value for ccw moment when viewed from +y axis
c       cmz....moment about z-axis through the moment center (xmc,ymc,zmc)
c              positive value for ccw moment when viewed from +z axis
c       cx.....x-component of force
c              positive for force in +x direction
c       cy.....y-component of force
c              positive for force in +y direction
c       cz.....z-component of force
c              positive for force in +z direction
c
c**************************************************************************
c
      cpc   = 2.e0/(gamma*xmach*xmach)
      const = 4./(reue*xmach)
c
      cosa  = cos(alpha)
      sina  = sin(alpha)
      cosb  = cos(beta)
      sinb  = sin(beta)
c
      cl    = 0.e0 
      cd    = 0.e0 
      cz    = 0.e0
      cy    = 0.e0 
      cx    = 0.e0
      cmy   = 0.e0 
      cmx   = 0.e0
      cmz   = 0.e0
      chd   = 0.e0 
      swet  = 0.e0 
      cdp   = 0.e0
      cdv   = 0.e0
c
c******************************************************************
c     forces on k=constant surfaces
c******************************************************************
c
      do 1 kk=1,2
c
      if (kk.eq.1) then
         nseg = nbck0(nbl) 
         k    = 1
         kc   = 1
      else
         nseg = nbckdim(nbl)
         k    = kdim
         kc   = kdim - 1
      end if
c
      do 5 ns = 1,nseg
c
      if (kbcinfo(nbl,ns,6,kk) .gt. 0 ) then
c
      ist = kbcinfo(nbl,ns,2,kk)
      ifn = kbcinfo(nbl,ns,3,kk) - 1
      jst = kbcinfo(nbl,ns,4,kk) 
      jfn = kbcinfo(nbl,ns,5,kk) - 1
c
      sgn = 1.0
      if(kk .gt. 1) sgn = -1.0
c
      do 10 i=ist,ifn 
      cxl   = 0.e0
      cyl   = 0.e0
      czl   = 0.e0
      cmyl  = 0.e0
      cmxl  = 0.e0
      cmzl  = 0.e0
      chdl  = 0.e0 
      xas   = 0.e0
      yas   = 0.e0
      zas   = 0.e0
      swetl = 0.e0
      cdvl  = 0.e0
      cdpl  = 0.e0
c
      do 15 j=jst,jfn
c 
      xa  = .25e0*( x(j,k,i) + x(j+1,k,i) + x(j,k,i+1) + x(j+1,k,i+1) )
      ya  = .25e0*( y(j,k,i) + y(j+1,k,i) + y(j,k,i+1) + y(j+1,k,i+1) )
      za  = .25e0*( z(j,k,i) + z(j+1,k,i) + z(j,k,i+1) + z(j+1,k,i+1) )
c 
      dcp = -(qk0(j,i,5,kk+kk-1)/p0-1.e0)*cpc*sk(j,k,i,4)
      dcx = dcp*sk(j,k,i,1)*sgn 
      dcy = dcp*sk(j,k,i,2)*sgn 
      dcz = dcp*sk(j,k,i,3)*sgn 
c     pressure drag
      dcdpl = dcx*cosa*cosb-dcy*sinb+dcz*sina*cosb
c
      dcdvl = 0. 
      if (ivisc(3).gt.0 .and. (abs(kbcinfo(nbl,ns,1,kk)).eq.2004 .or.
     .    abs(kbcinfo(nbl,ns,1,kk)).eq.2014 .or.
     .    abs(kbcinfo(nbl,ns,1,kk)).eq.2024 .or.
     .    abs(kbcinfo(nbl,ns,1,kk)).eq.2034 .or.
     .    abs(kbcinfo(nbl,ns,1,kk)).eq.2016)) then
         dcxp  = dcx
         dcyp  = dcy
         dczp  = dcz
         urel  = ub(j,kc,i)
         vrel  = vb(j,kc,i)
         wrel  = wb(j,kc,i)
         if (iuns.gt.0 .and. kbcinfo(nbl,ns,1,kk) .gt. 0) then
            urel  = ub(j,kc,i) - xtbk(j,i,1,kk)
            vrel  = vb(j,kc,i) - xtbk(j,i,2,kk)
            wrel  = wb(j,kc,i) - xtbk(j,i,3,kk)
         end if
         tau   = vmuk(j,i,kk)*const/vol(j,kc,i)*sk(j,k,i,4)**2
         vnorm = (urel*sk(j,k,i,1)+vrel*sk(j,k,i,2)
     .           +wrel*sk(j,k,i,3))*sgn
         dcx   = tau*(urel-vnorm*sk(j,k,i,1)*sgn)
         dcz   = tau*(wrel-vnorm*sk(j,k,i,3)*sgn)
         dcy   = tau*(vrel-vnorm*sk(j,k,i,2)*sgn)
c        viscous drag
         dcdvl  = dcx*cosa*cosb-dcy*sinb+dcz*sina*cosb
         dcx = dcxp + dcx
         dcy = dcyp + dcy
         dcz = dczp + dcz
      end if
c
c     only use contributions from points with interface (solid surface)
c     boundary conditions, and only from those points not blanked out
c
      fact = bck(j,i,kk)*blank(j,kc,i)
      dcx  = dcx*fact
      dcy  = dcy*fact
      dcz  = dcz*fact
      dcdpl = dcdpl*fact
      dcdvl = dcdvl*fact
      dsw  = sk(j,k,i,4)*fact
c
      swetl = swetl+dsw 
      cxl   = cxl+dcx 
      cyl   = cyl+dcy 
      czl   = czl+dcz 
      cdpl  = cdpl + dcdpl
      cdvl  = cdvl + dcdvl
      cmyl  = cmyl-dcz*(xa-xmc)+dcx*(za-zmc)
      cmxl  = cmxl+dcz*(ya-ymc)-dcy*(za-zmc)
      cmzl  = cmzl+dcy*(xa-xmc)-dcx*(ya-ymc)
   15 continue
c
c     integrated values
      cds   = cxl*cosa*cosb-cyl*sinb+czl*sina*cosb
      cls   =-cxl*sina+czl*cosa
      cmys  = cmyl
      cmxs  = cmxl
      cmzs  = cmzl 
      swets = swetl 
      cl    = cl+cls 
      cd    = cd+cds 
      cz    = cz+czl
      cy    = cy+cyl
      cx    = cx+cxl 
      cmy   = cmy+cmys
      cmx   = cmx+cmxs
      cmz   = cmz+cmzs 
      chd   = 1.0
      swet  = swet+swets 
      cdp   = cdp + cdpl
      cdv   = cdv + cdvl
c
   10 continue
      end if
    5 continue
    1 continue
c
c******************************************************************
c     forces on j=constant surfaces
c******************************************************************
c
c
      do 6 jj=1,2
c
      if (jj.eq.1) then
         nseg = nbcj0(nbl)
         j    = 1
         jc   = 1
      else
         nseg = nbcjdim(nbl)
         j    = jdim
         jc   = jdim - 1
      end if
c
      do 35 ns = 1,nseg
c
      if (jbcinfo(nbl,ns,6,jj) .gt. 0 ) then
c
      ist = jbcinfo(nbl,ns,2,jj)
      ifn = jbcinfo(nbl,ns,3,jj) - 1
      kst = jbcinfo(nbl,ns,4,jj)
      kfn = jbcinfo(nbl,ns,5,jj) - 1
c
      sgn = 1.0
      if(jj .gt. 1) sgn = -1.0
c
      do 40 i=ist,ifn
      cxl   = 0.e0
      cyl   = 0.e0
      czl   = 0.e0
      cmyl  = 0.e0
      cmxl  = 0.e0
      cmzl  = 0.e0
      chdl  = 0.e0
      swetl = 0.e0
      cdvl  = 0.e0
      cdpl  = 0.e0
c
      do 45 k=kst,kfn
c
      xa  = .25e0*( x(j,k,i) + x(j,k+1,i) + x(j,k,i+1) + x(j,k+1,i+1) )
      ya  = .25e0*( y(j,k,i) + y(j,k+1,i) + y(j,k,i+1) + y(j,k+1,i+1) )
      za  = .25e0*( z(j,k,i) + z(j,k+1,i) + z(j,k,i+1) + z(j,k+1,i+1) )
c
      dcp = -(qj0(k,i,5,jj+jj-1)/p0-1.e0)*cpc*sj(j,k,i,4)
      dcx = dcp*sj(j,k,i,1)*sgn
      dcy = dcp*sj(j,k,i,2)*sgn
      dcz = dcp*sj(j,k,i,3)*sgn
c     pressure drag
      dcdpl = dcx*cosa*cosb-dcy*sinb+dcz*sina*cosb
c
      dcdvl = 0.
      if (ivisc(2).gt.0 .and. (abs(jbcinfo(nbl,ns,1,jj)).eq.2004 .or.
     .    abs(jbcinfo(nbl,ns,1,jj)).eq.2014 .or.
     .    abs(jbcinfo(nbl,ns,1,jj)).eq.2024 .or.
     .    abs(jbcinfo(nbl,ns,1,jj)).eq.2034 .or.
     .    abs(jbcinfo(nbl,ns,1,jj)).eq.2016)) then
         dcxp  = dcx
         dcyp  = dcy
         dczp  = dcz
         urel  = ub(jc,k,i)
         vrel  = vb(jc,k,i)
         wrel  = wb(jc,k,i)
         if (iuns.gt.0 .and. jbcinfo(nbl,ns,1,jj) .gt. 0) then
            urel  = ub(jc,k,i) - xtbj(k,i,1,jj)
            vrel  = vb(jc,k,i) - xtbj(k,i,2,jj)
            wrel  = wb(jc,k,i) - xtbj(k,i,3,jj)
         end if
         tau   = vmuj(k,i,jj)*const/vol(jc,k,i)*sj(j,k,i,4)**2
         vnorm = (urel*sj(j,k,i,1)+vrel*sj(j,k,i,2)
     .           +wrel*sj(j,k,i,3))*sgn
         dcx   = tau*(urel-vnorm*sj(j,k,i,1)*sgn)
         dcz   = tau*(wrel-vnorm*sj(j,k,i,3)*sgn)
         dcy   = tau*(vrel-vnorm*sj(j,k,i,2)*sgn)
c        viscous drag
         dcdvl  = dcx*cosa*cosb-dcy*sinb+dcz*sina*cosb
         dcx = dcxp + dcx
         dcy = dcyp + dcy
         dcz = dczp + dcz
      end if
c
c     only use contributions from points with interface (solid surface)
c     boundary conditions, and only from those points not blanked out
c
      fact = bcj(k,i,jj)*blank(jc,k,i)
      dcx  = dcx*fact
      dcy  = dcy*fact
      dcz  = dcz*fact
      dcdpl = dcdpl*fact
      dcdvl = dcdvl*fact
      dsw  = sj(j,k,i,4)*fact
c
      swetl = swetl+dsw
      cxl   = cxl+dcx
      cyl   = cyl+dcy
      czl   = czl+dcz
      cdpl  = cdpl + dcdpl
      cdvl  = cdvl + dcdvl
      cmyl  = cmyl-dcz*(xa-xmc)+dcx*(za-zmc)
      cmxl  = cmxl+dcz*(ya-ymc)-dcy*(za-zmc)
      cmzl  = cmzl+dcy*(xa-xmc)-dcx*(ya-ymc)
c
   45 continue
c
c     integrated values
      cds   = cxl*cosa*cosb-cyl*sinb+czl*sina*cosb
      cls   =-cxl*sina+czl*cosa
      cmys  = cmyl
      cmxs  = cmxl
      cmzs  = cmzl
      swets = swetl
      cl    = cl+cls
      cd    = cd+cds
      cz    = cz+czl
      cy    = cy+cyl
      cx    = cx+cxl
      cmy   = cmy+cmys
      cmx   = cmx+cmxs
      cmz   = cmz+cmzs
      chd   = 1.0
      swet  = swet+swets
      cdp   = cdp + cdpl
      cdv   = cdv + cdvl
c
   40 continue
      end if
   35 continue
    6 continue
c
c******************************************************************
c     forces on i=constant surfaces
c******************************************************************
c
      do 7 ii=1,2
c
      if (ii.eq.1) then
         nseg = nbci0(nbl)
         i    = 1
         ic   = 1
      else
         nseg = nbcidim(nbl)
         i    = idim
         ic   = idim - 1
      end if
c
      do 65 ns = 1,nseg
c
      if (ibcinfo(nbl,ns,6,ii) .gt. 0 ) then
c
      jst = ibcinfo(nbl,ns,2,ii)
      jfn = ibcinfo(nbl,ns,3,ii) - 1
      kst = ibcinfo(nbl,ns,4,ii)
      kfn = ibcinfo(nbl,ns,5,ii) - 1
c
      sgn = 1.0
      if(ii .gt. 1) sgn = -1.0
c
      do 70 j=jst,jfn
      cxl   = 0.e0
      cyl   = 0.e0
      czl   = 0.e0
      cmyl  = 0.e0
      cmxl  = 0.e0
      cmzl  = 0.e0
      chdl  = 0.e0
      swetl = 0.e0
      cdvl  = 0.e0
      cdpl  = 0.e0
c
      do 75 k=kst,kfn
c
      xa  = .25e0*( x(j,k,i) + x(j,k+1,i) + x(j+1,k,i) + x(j+1,k+1,i) )
      ya  = .25e0*( y(j,k,i) + y(j,k+1,i) + y(j+1,k,i) + y(j+1,k+1,i) )
      za  = .25e0*( z(j,k,i) + z(j,k+1,i) + z(j+1,k,i) + z(j+1,k+1,i) )
c
      dcp = -(qi0(j,k,5,ii+ii-1)/p0-1.e0)*cpc*si(j,k,i,4)
      dcx = dcp*si(j,k,i,1)*sgn
      dcy = dcp*si(j,k,i,2)*sgn
      dcz = dcp*si(j,k,i,3)*sgn
c     pressure drag
      dcdpl = dcx*cosa*cosb-dcy*sinb+dcz*sina*cosb
c
      dcdvl = 0.
      if (ivisc(1).gt.0 .and. (abs(ibcinfo(nbl,ns,1,ii)).eq.2004 .or.
     .    abs(ibcinfo(nbl,ns,1,ii)).eq.2014 .or.
     .    abs(ibcinfo(nbl,ns,1,ii)).eq.2024 .or.
     .    abs(ibcinfo(nbl,ns,1,ii)).eq.2034 .or.
     .    abs(ibcinfo(nbl,ns,1,ii)).eq.2016)) then
         dcxp  = dcx
         dcyp  = dcy
         dczp  = dcz
         urel  = ub(j,k,ic)
         vrel  = vb(j,k,ic)
         wrel  = wb(j,k,ic)
         if (iuns.gt.0 .and. ibcinfo(nbl,ns,1,ii) .gt. 0) then
            urel  = ub(j,k,ic) - xtbi(j,k,1,ii)
            vrel  = vb(j,k,ic) - xtbi(j,k,2,ii)
            wrel  = wb(j,k,ic) - xtbi(j,k,3,ii)
         end if
         tau   = vmui(j,k,ii)*const/vol(j,k,ic)*si(j,k,i,4)**2
         vnorm = (urel*si(j,k,i,1)+vrel*si(j,k,i,2)
     .           +wrel*si(j,k,i,3))*sgn
         dcx   = tau*(urel-vnorm*si(j,k,i,1)*sgn)
         dcz   = tau*(wrel-vnorm*si(j,k,i,3)*sgn)
         dcy   = tau*(vrel-vnorm*si(j,k,i,2)*sgn)
c        viscous drag
         dcdvl  = dcx*cosa*cosb-dcy*sinb+dcz*sina*cosb
         dcx = dcxp + dcx
         dcy = dcyp + dcy
         dcz = dczp + dcz
      end if
c
c     only use contributions from points with interface (solid surface)
c     boundary conditions, and only from those points not blanked out
c
      fact = bci(j,k,ii)*blank(j,k,ic)
      dcx  = dcx*fact
      dcy  = dcy*fact
      dcz  = dcz*fact
      dcdpl = dcdpl*fact
      dcdvl = dcdvl*fact
      dsw  = si(j,k,i,4)*fact
c
      swetl = swetl+dsw
      cxl   = cxl+dcx
      cyl   = cyl+dcy
      czl   = czl+dcz
      cdpl  = cdpl + dcdpl
      cdvl  = cdvl + dcdvl
      cmyl  = cmyl-dcz*(xa-xmc)+dcx*(za-zmc)
      cmxl  = cmxl+dcz*(ya-ymc)-dcy*(za-zmc)
      cmzl  = cmzl+dcy*(xa-xmc)-dcx*(ya-ymc)
   75 continue
c
c     integrated values
      cds   = cxl*cosa*cosb-cyl*sinb+czl*sina*cosb
      cls   =-cxl*sina+czl*cosa
      cmys  = cmyl
      cmxs  = cmxl
      cmzs  = cmzl
      swets = swetl
      cl    = cl+cls
      cd    = cd+cds
      cz    = cz+czl
      cy    = cy+cyl
      cx    = cx+cxl
      cmy   = cmy+cmys
      cmx   = cmx+cmxs
      cmz   = cmz+cmzs
      chd   = 1.0
      swet  = swet+swets
      cdp   = cdp + cdpl
      cdv   = cdv + cdvl
c
   70 continue
      end if
   65 continue
    7 continue
c
c     force cmx and cmy=0 if 2-D; then user doesn't have to worry about
c     inputting "correct" cell center location for the non-2-D-direction
      if(i2d .eq. 1) then
        cmx=0.
        cmz=0.
      end if
c
c     note: the forces above are computed in the standard cfl3d coordinate
c     system where z is "up". if ialph = 1, then must map these into a system 
c     with y "up" so that when output, they are consistent with the input
c     grid. If x,y,z is the standard cfl3d system and x',y',z' is the 
c     input grid system (with x'=x), then for ialph .ne. 0 the transform
c     z = y', y = -z' maps the input grid onto the standard cfl3d system
c     (see also subroutine rp3d). Conversely z' = -y, y' = z' maps the
c     cfl3d system back into the input system (this is what is used below,
c     where the ' has been dropped from the left hand side)
c
      if (ialph.gt.0) then
         temp = cmz
         cmz  = -cmy
         cmy  = temp
         temp = cz
         cz   = -cy
         cy   = temp
      end if
c
      return
      end 
