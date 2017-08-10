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
      subroutine bc2006(jdim,kdim,idim,q,qj0,qk0,qi0,sj,sk,si,bcj,bck,
     .                  bci,xtbj,xtbk,xtbi,atbj,atbk,atbi,ista,iend,
     .                  jsta,jend,ksta,kend,nface,tursav,tj0,tk0,
     .                  ti0,vist3d,vj0,vk0,vi0,mdim,ndim,bcdata,
     .                  filname,iuns,x,y,z,nblc,jdimc,kdimc,idimc,
     .                  qj0c,qk0c,qi0c,ichk,nou,bou,nbuf,ibufdim,myid,
     .                  nummem)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Set pressure via radial equilibrium condition;
c     extrapolate other quantities
c
c     The radial equilibrium equation is dp/dr = rho*vtheta**2/r
c     where vtheta is the circumferential velocity component. given 
c     the pressure at one radius, the radial equilibrium equation is
c     integrated to give the pressure at all other radii. Integration 
c     is carried out via the trapezoidal rule. 
c 
c     contents of bcdata array for this boundary condition:
c   
c     bcdata(m,n,ip,1) = block from which integration of equilibrium
c                     equation is to be continued (if > 0)
c     bcdata(m,n,ip,2) = p/pinf at the initial radial position
c     bcdata(m,n,ip,3) = indicates coordinate index/direction for 
c                     integration (see below)
c     bcdata(m,n,ip,4) = indicates which physical coordinate is held 
c                     constant for radial integration: 1=x, 2=y, 3=z
c
c     only the bcdata values on the boundaries m=1/mdim or
c     n=1/ndim are actually used - interior values are ignored
c
c     lijk = abs(int(bcdata(m,n,ip,3)), indicates computational
c            coordinate direction in which the radial equilibrium
c            equation dp/dr = rho*vtheta**2/r is to be integrated:
c            1=i, 2=j, 3=k
c     ldir = int(bcdata(m,n,ip,3)/abs(int(bcdata(m,n,ip,3)), indicates in
c            which direction integration is to be carried out:
c            from 1 to max (ldir>0) or max to 1 (ldir<0)
c
c     nblc = bcdata(m,n,ip,1), block from which radial integration is
c            to be continued, if non-zero. qi0c/qj0c/qk0c are the 
c            ghost-cell data in  block nblc, which has dimensions
c            idimc,jdimc,kdimc
c
c     a check is made to determine if the data in the ghost cells
c     of the adjacent block have not yet been set. this occurs in
c     the first call to the subroutine if the block the integration
c     is continued from has a larger block number than the current 
c     block. if not yet set, then the nearest interior point in the
c     current block is used instead (this check is needed only for 
c     cases where integration is continued from another block).
c
c                            
c     NOTES:  
c
c     this boundary condition assumes that one coordinate 
c     direction on the block face is essentially radial and
c     the other is essentially circumferential, and that 
c     the integration is being carried out in the radial 
c     direction. as there is no way to verify this in the code,
c     care must be exercised when using this bc to insure that
c     this restriction and the ones below are met.
c
c     must have bcdata(m,n,ip,3) = bcdata(1,1,ip,3) for all boundary 
c     indicies of m,n, i.e. direction of integration cannot
c     change from point to point
c
c     continuation of the radial equilibrium condition through
c     block boundaries is restricted to blocks that have the
c     SAME orientation. For example, if the equilibrium condition
c     is to be continued through a k-boundary from an adjacent
c     block, then i and j must run in the same direction in both
c     blocks. This also implies that the k-boundary must be k=1
c     in one block and k=kdim in the second (e.g. the boundary
c     cannot be k=1 in both blocks)
c
c     the segment must run the entire length of the block face in
c     the direction in which the integration is being carried
c     out.  e.g. if the  integration is being carried out in the 
c     k-direction, then must have ksta = 1 and kend = kdim (this
c     restriction applies only if the integration is being
c     continued from another block)
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
     .          qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4),
     .          x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
      dimension qi0c(jdimc,kdimc,5,4),qj0c(kdimc,idimc-1,5,4),
     .          qk0c(jdimc,idimc-1,5,4)
      dimension sk(jdim,kdim,idim-1,5),si(jdim,kdim,idim,5),
     .          sj(jdim,kdim,idim-1,5)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
      dimension xtbj(kdim,idim-1,3,2),xtbk(jdim,idim-1,3,2),
     .          xtbi(jdim,kdim,3,2),atbj(kdim,idim-1,3,2),
     .          atbk(jdim,idim-1,3,2),atbi(jdim,kdim,3,2)
      dimension bcdata(mdim,ndim,2,12)
      dimension tursav(jdim,kdim,idim,nummem),tj0(kdim,idim-1,nummem,4),
     .          tk0(jdim,idim-1,nummem,4),ti0(jdim,kdim,nummem,4),
     .          vj0(kdim,idim-1,1,4),vk0(jdim,idim-1,1,4),
     .          vi0(jdim,kdim,1,4),vist3d(jdim,kdim,idim)
c
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .        tur10(7)
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /reyue/ reue,tinf,ivisc(3)
      common /sklton/ isklton
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
      ipp    = 1
c
      lijk = abs(int(real(bcdata(1,1,1,3))))
      ldir = int(real(bcdata(1,1,1,3)))/abs(int(real(bcdata(1,1,1,3))))
c
      xfact = 1.0
      yfact = 1.0
      zfact = 1.0
      if (bcdata(1,1,1,4) .eq. 1) then
         xfact = 0.0
      else if (bcdata(1,1,1,4) .eq. 2) then
         yfact = 0.0
      else if (bcdata(1,1,1,4) .eq. 3) then
         zfact = 0.0
      else
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),*)' stopping...bcdata(4) = ',
     .   bcdata(1,1,1,4),' should be +1, +2, or +3' 
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
c
c            * * * * * * * * * * * * * * * * * * * * * *
c            * standard boundary condition bctype=2006 *
c            * * * * * * * * * * * * * * * * * * * * * *
c
c******************************************************************************
c      j=1 boundary   set pressure via radial equilibrium           bctype 2006
c******************************************************************************
c
      if (nface.eq.3) then
c
      j = 1
c
c     extrapolate all primatives except pressure
c
      do 100 l=1,4
      do 100 i=ista,iend1
      do 100 k=ksta,kend1
      qj0(k,i,l,1) = q(1,k,i,l)
      qj0(k,i,l,2) = qj0(k,i,l,1)
      bcj(k,i,1) = 0.0
  100 continue
c
      if (lijk .eq. 1) then
c
c        radial integration of dp/dr in i-direction
c        (k is the theta-direction)
c
         if (ldir .gt. 0) then
            i0   = ista
            is   = ista + 1
            ie   = iend
            ii0  = 1
            iinc = 1
            iic  = idimc-1
         else 
            i0   = iend
            is   = iend - 1
            ie   = ista 
            ii0  = iend - ista
            iinc = 0
            iic  = 1
         end if
c
         if (nblc .ne. 0) then
c
c           set data for continuation
c
            if (ichk .gt. 0) then
               do 110 k=ksta,kend1
               kk = k-ksta+1
               bcdata(kk,ii0,ipp,2) = qj0c(k,iic,5,1)/p0
  110          continue
            else
               do 120 k=ksta,kend1
               kk = k-ksta+1
               bcdata(kk,ii0,ipp,2) = q(j,k,ii0,5)/p0
  120          continue
            end if
         end if
c
         do 130 k=ksta,kend1
         kk   = k-ksta+1
         pm1  = bcdata(kk,ii0,ipp,2)*p0
         i    = i0
         kp   = k  + 1
         jp   = j  + 1
         xav0 = 0.25*(x(j,k,i)  + x(j,kp,i)
     .              + x(jp,k,i) + x(jp,kp,i))
         yav0 = 0.25*(y(j,k,i)  + y(j,kp,i)
     .              + y(jp,k,i) + y(jp,kp,i))
         zav0 = 0.25*(z(j,k,i)  + z(j,kp,i)
     .              + z(jp,k,i) + z(jp,kp,i))
         xav0 = xav0*xfact
         yav0 = yav0*yfact
         zav0 = zav0*zfact
         r0   = sqrt(xav0**2 + yav0**2 + zav0**2)
         do 130 i=is,ie,ldir
         im   = i  - iinc
         kp   = k  + 1
         jp   = j  + 1
         xav1 = 0.25*(x(j,k,i)  + x(j,kp,i)
     .              + x(jp,k,i) + x(jp,kp,i))
         yav1 = 0.25*(y(j,k,i)  + y(j,kp,i)
     .              + y(jp,k,i) + y(jp,kp,i))
         zav1 = 0.25*(z(j,k,i)  + z(j,kp,i)
     .              + z(jp,k,i) + z(jp,kp,i))
         xav1 = xav1*xfact
         yav1 = yav1*yfact
         zav1 = zav1*zfact
         r1   = sqrt(xav1**2 + yav1**2 + zav1**2)
         rho  = q(j,k,im,1 )
         dnxav = 0.5*(sk(j,k,im,1) + sk(j,k+1,im,1))
         dnyav = 0.5*(sk(j,k,im,2) + sk(j,k+1,im,2))
         dnzav = 0.5*(sk(j,k,im,3) + sk(j,k+1,im,3))
         qthet = q(j,k,im,2)*dnxav + q(j,k,im,3)*dnyav
     .         + q(j,k,im,4)*dnzav
         rav   = 0.5*(r1 + r0)
         dr    = r1 - r0
         dp    = dr*rho*qthet**2/rav
         qj0(k,im,5,1) = pm1 + dp
         qj0(k,im,5,2) = qj0(k,im,5,1)
         r0  = r1
         pm1 = qj0(k,im,5,1)
  130    continue
      end if
c
      if (lijk .eq. 3) then
c
c        radial integration of dp/dr in k-direction
c        (i is the theta-direction)
c
         if (ldir .gt. 0) then
            k0   = ksta
            ks   = ksta + 1
            ke   = kend
            kk0  = 1
            kinc = 1
            kkc  = kdimc-1
         else 
            k0   = kend
            ks   = kend - 1
            ke   = ksta 
            kk0  = kend - ksta
            kinc = 0
            kkc  = 1
         end if
c
         if (nblc .ne. 0) then
c
c           set data for continuation
c
            if (ichk .gt. 0) then
               do 140 i=ista,iend1
               ii = i-ista+1
               bcdata(kk0,ii,ipp,2) = qj0c(kkc,i,5,1)/p0
  140          continue
            else
               do 150 i=ista,iend1
               ii = i-ista+1
               bcdata(kk0,ii,ipp,2) = q(j,kk0,i,5)/p0
  150          continue
            end if
         end if
c
         do 160 i=ista,iend1
         ii   = i-ista+1
         pm1  = bcdata(kk0,ii,ipp,2)*p0
         k    = k0
         ip   = i  + 1
         jp   = j  + 1
         xav0 = 0.25*(x(j,k,i)  + x(j,k,ip)
     .              + x(jp,k,i) + x(jp,k,ip))
         yav0 = 0.25*(y(j,k,i)  + y(j,k,ip)
     .              + y(jp,k,i) + y(jp,k,ip))
         zav0 = 0.25*(z(j,k,i)  + z(j,k,ip)
     .              + z(jp,k,i) + z(jp,k,ip))
         xav0 = xav0*xfact
         yav0 = yav0*yfact
         zav0 = zav0*zfact
         r0   = sqrt(xav0**2 + yav0**2 + zav0**2)
         do 160 k=ks,ke,ldir
         km   = k - kinc
         ip   = i  + 1
         jp   = j  + 1
         xav1 = 0.25*(x(j,k,i)  + x(j,k,ip)
     .              + x(jp,k,i) + x(jp,k,ip))
         yav1 = 0.25*(y(j,k,i)  + y(j,k,ip)
     .              + y(jp,k,i) + y(jp,k,ip))
         zav1 = 0.25*(z(j,k,i)  + z(j,k,ip)
     .              + z(jp,k,i) + z(jp,k,ip))
         xav1 = xav1*xfact
         yav1 = yav1*yfact
         zav1 = zav1*zfact
         r1   = sqrt(xav1**2 + yav1**2 + zav1**2)
         rho  = q(j,km,i,1 )
         dnxav = 0.5*(si(j,km,i,1) + si(j,km,i+1,1))
         dnyav = 0.5*(si(j,km,i,2) + si(j,km,i+1,2))
         dnzav = 0.5*(si(j,km,i,3) + si(j,km,i+1,3))
         qthet = q(j,km,i,2)*dnxav + q(j,km,i,3)*dnyav
     .         + q(j,km,i,4)*dnzav
         rav   = 0.5*(r1 + r0)
         dr    = r1 - r0
         dp    = dr*rho*qthet**2/rav
         qj0(km,i,5,1) = pm1 + dp
         qj0(km,i,5,2) = qj0(km,i,5,1)
         r0  = r1
         pm1 = qj0(km,i,5,1)
  160    continue
      end if
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 170 i=ista,iend1
        do 170 k=ksta,kend1
          vj0(k,i,1,1) = vist3d(1,k,i)
          vj0(k,i,1,2) = vist3d(1,k,i)
  170   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 180 i=ista,iend1
        do 180 k=ksta,kend1
          tj0(k,i,l,1) = tursav(1,k,i,l)
          tj0(k,i,l,2) = tursav(1,k,i,l)
  180   continue
        enddo
      end if
      end if
      end if
c
c******************************************************************************
c      j=jdim boundary   set pressure via radial equilibrium        bctype 2006
c******************************************************************************
c 
      if (nface.eq.4) then
c
      j = jdim - 1
c
c     extrapolate all primatives except pressure
c
      do 200 l=1,4
      do 200 i=ista,iend1
      do 200 k=ksta,kend1
      qj0(k,i,l,3) = q(jdim1,k,i,l)
      qj0(k,i,l,4) = qj0(k,i,l,3)
      bcj(k,i,2) = 0.0
  200 continue
c
      if (lijk .eq. 1) then
c
c        radial integration of dp/dr in i-direction
c        (k is the theta-direction)
c
         if (ldir .gt. 0) then
            i0   = ista
            is   = ista + 1
            ie   = iend
            ii0  = 1
            iinc = 1
            iic  = idimc-1
         else 
            i0   = iend
            is   = iend - 1
            ie   = ista 
            ii0  = iend - ista
            iinc = 0
            iic  = 1
         end if
c
         if (nblc .ne. 0) then
c
c           set data for continuation
c
            if (ichk .gt. 0) then
               do 210 k=ksta,kend1
               kk = k-ksta+1
               bcdata(kk,ii0,ipp,2) = qj0c(k,iic,5,3)/p0
  210          continue
            else
               do 220 k=ksta,kend1
               kk = k-ksta+1
               bcdata(kk,ii0,ipp,2) = q(j,k,ii0,5)/p0
  220          continue
            end if
         end if
c
         do 230 k=ksta,kend1
         kk   = k-ksta+1
         pm1  = bcdata(kk,ii0,ipp,2)*p0
         i    = i0
         kp   = k  + 1
         jp   = j  + 1
         xav0 = 0.25*(x(j,k,i)  + x(j,kp,i)
     .              + x(jp,k,i) + x(jp,kp,i))
         yav0 = 0.25*(y(j,k,i)  + y(j,kp,i)
     .              + y(jp,k,i) + y(jp,kp,i))
         zav0 = 0.25*(z(j,k,i)  + z(j,kp,i)
     .              + z(jp,k,i) + z(jp,kp,i))
         xav0 = xav0*xfact
         yav0 = yav0*yfact
         zav0 = zav0*zfact
         r0   = sqrt(xav0**2 + yav0**2 + zav0**2)
         do 230 i=is,ie,ldir
         im   = i  - iinc
         kp   = k  + 1
         jp   = j  + 1
         xav1 = 0.25*(x(j,k,i)  + x(j,kp,i)
     .              + x(jp,k,i) + x(jp,kp,i))
         yav1 = 0.25*(y(j,k,i)  + y(j,kp,i)
     .              + y(jp,k,i) + y(jp,kp,i))
         zav1 = 0.25*(z(j,k,i)  + z(j,kp,i)
     .              + z(jp,k,i) + z(jp,kp,i))
         xav1 = xav1*xfact
         yav1 = yav1*yfact
         zav1 = zav1*zfact
         r1   = sqrt(xav1**2 + yav1**2 + zav1**2)
         rho  = q(j,k,im,1 )
         dnxav = 0.5*(sk(j,k,im,1) + sk(j,k+1,im,1))
         dnyav = 0.5*(sk(j,k,im,2) + sk(j,k+1,im,2))
         dnzav = 0.5*(sk(j,k,im,3) + sk(j,k+1,im,3))
         qthet = q(j,k,im,2)*dnxav + q(j,k,im,3)*dnyav
     .         + q(j,k,im,4)*dnzav
         rav   = 0.5*(r1 + r0)
         dr    = r1 - r0
         dp    = dr*rho*qthet**2/rav
         qj0(k,im,5,3) = pm1 + dp
         qj0(k,im,5,4) = qj0(k,im,5,3)
         r0  = r1
         pm1 = qj0(k,im,5,3)
  230    continue
      end if
c
      if (lijk .eq. 3) then
c
c        radial integration of dp/dr in k-direction
c        (i is the theta-direction)
c
         if (ldir .gt. 0) then
            k0   = ksta
            ks   = ksta + 1
            ke   = kend
            kk0  = 1
            kinc = 1
            kkc  = kdimc-1
         else 
            k0   = kend
            ks   = kend - 1
            ke   = ksta 
            kk0  = kend - ksta
            kinc = 0
            kkc  = 1
         end if
c
         if (nblc .ne. 0) then
c
c           set data for continuation
c
            if (ichk .gt. 0) then
               do 240 i=ista,iend1
               ii = i-ista+1
               bcdata(kk0,ii,ipp,2) = qj0c(kkc,i,5,3)/p0
  240          continue
            else
               do 250 i=ista,iend1
               ii = i-ista+1
               bcdata(kk0,ii,ipp,2) = q(j,kk0,i,5)/p0
  250          continue
            end if
         end if
c
         do 260 i=ista,iend1
         ii   = i-ista+1
         pm1  = bcdata(kk0,ii,ipp,2)*p0
         k    = k0
         ip   = i  + 1
         jp   = j  + 1
         xav0 = 0.25*(x(j,k,i)  + x(j,k,ip)
     .              + x(jp,k,i) + x(jp,k,ip))
         yav0 = 0.25*(y(j,k,i)  + y(j,k,ip)
     .              + y(jp,k,i) + y(jp,k,ip))
         zav0 = 0.25*(z(j,k,i)  + z(j,k,ip)
     .              + z(jp,k,i) + z(jp,k,ip))
         xav0 = xav0*xfact
         yav0 = yav0*yfact
         zav0 = zav0*zfact
         do 260 k=ks,ke,ldir
         km   = k - kinc
         ip   = i  + 1
         jp   = j  + 1
         xav1 = 0.25*(x(j,k,i)  + x(j,k,ip)
     .              + x(jp,k,i) + x(jp,k,ip))
         yav1 = 0.25*(y(j,k,i)  + y(j,k,ip)
     .              + y(jp,k,i) + y(jp,k,ip))
         zav1 = 0.25*(z(j,k,i)  + z(j,k,ip)
     .              + z(jp,k,i) + z(jp,k,ip))
         xav1 = xav1*xfact
         yav1 = yav1*yfact
         zav1 = zav1*zfact
         r1   = sqrt(xav1**2 + yav1**2 + zav1**2)
         rho  = q(j,km,i,1 )
         dnxav = 0.5*(si(j,km,i,1) + si(j,km,i+1,1))
         dnyav = 0.5*(si(j,km,i,2) + si(j,km,i+1,2))
         dnzav = 0.5*(si(j,km,i,3) + si(j,km,i+1,3))
         qthet = q(j,km,i,2)*dnxav + q(j,km,i,3)*dnyav
     .         + q(j,km,i,4)*dnzav
         rav   = 0.5*(r1 + r0)
         dr    = r1 - r0
         dp    = dr*rho*qthet**2/rav
         qj0(km,i,5,3) = pm1 + dp
         qj0(km,i,5,4) = qj0(km,i,5,3)
         r0  = r1
         pm1 = qj0(km,i,5,3)
  260    continue
      end if
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 270 i=ista,iend1
        do 270 k=ksta,kend1
          vj0(k,i,1,3) = vist3d(jdim1,k,i)
          vj0(k,i,1,4) = vist3d(jdim1,k,i)
  270   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 280 i=ista,iend1
        do 280 k=ksta,kend1
          tj0(k,i,l,3) = tursav(jdim1,k,i,l)
          tj0(k,i,l,4) = tursav(jdim1,k,i,l)
  280   continue
        enddo
      end if
      end if
      end if
c
c******************************************************************************
c      k=1 boundary   set pressure via radial equilibrium           bctype 2006
c******************************************************************************
c
      if (nface.eq.5) then
c
      k = 1
c
c     extrapolate all primatives except pressure
c
      do 300 l=1,4
      do 300 i=ista,iend1
      do 300 j=jsta,jend1
      qk0(j,i,l,1) = q(j,1,i,l)
      qk0(j,i,l,2) = qk0(j,i,l,1)
      bck(j,i,1) = 0.0
  300 continue
c
      if (lijk .eq. 1) then
c
c        radial integration of dp/dr in i-direction
c        (j is the theta-direction)
c
         if (ldir .gt. 0) then
            i0   = ista
            is   = ista + 1
            ie   = iend
            ii0  = 1
            iinc = 1
            iic  = idimc-1
         else 
            i0   = iend
            is   = iend - 1
            ie   = ista 
            ii0  = iend - ista
            iinc = 0
            iic  = 1
         end if
c
         if (nblc .ne. 0) then
c
c           set data for continuation
c
            if (ichk .gt. 0) then
               do 310 j=jsta,jend1
               jj = j-jsta+1
               bcdata(jj,ii0,ipp,2) = qk0c(j,iic,5,1)
  310          continue
            else
               do 320 j=jsta,jend1
               jj = j-jsta+1
               bcdata(jj,ii0,ipp,2) = q(j,k,ii0,5)/p0
  320          continue
            end if
         end if
c
         do 330 j=jsta,jend1
         jj   = j-jsta+1
         pm1  = bcdata(jj,ii0,ipp,2)*p0
         i    = i0
         kp   = k  + 1
         jp   = j  + 1
         xav0 = 0.25*(x(j,k,i)  + x(j,kp,i)
     .              + x(jp,k,i) + x(jp,kp,i))
         yav0 = 0.25*(y(j,k,i)  + y(j,kp,i)
     .              + y(jp,k,i) + y(jp,kp,i))
         zav0 = 0.25*(z(j,k,i)  + z(j,kp,i)
     .              + z(jp,k,i) + z(jp,kp,i))
         xav0 = xav0*xfact
         yav0 = yav0*yfact
         zav0 = zav0*zfact
         r0   = sqrt(xav0**2 + yav0**2 + zav0**2)
         do 330 i=is,ie,ldir
         im   = i  - iinc
         kp   = k  + 1
         jp   = j  + 1
         xav1 = 0.25*(x(j,k,i)  + x(j,kp,i)
     .              + x(jp,k,i) + x(jp,kp,i))
         yav1 = 0.25*(y(j,k,i)  + y(j,kp,i)
     .              + y(jp,k,i) + y(jp,kp,i))
         zav1 = 0.25*(z(j,k,i)  + z(j,kp,i)
     .              + z(jp,k,i) + z(jp,kp,i))
         xav1 = xav1*xfact
         yav1 = yav1*yfact
         zav1 = zav1*zfact
         r1   = sqrt(xav1**2 + yav1**2 + zav1**2)
         rho  = q(j,k,im,1 )
         dnxav = 0.5*(sj(j,k,im,1) + sj(j+1,k,im,1))
         dnyav = 0.5*(sj(j,k,im,2) + sj(j+1,k,im,2))
         dnzav = 0.5*(sj(j,k,im,3) + sj(j+1,k,im,3))
         qthet = q(j,k,im,2)*dnxav + q(j,k,im,3)*dnyav
     .         + q(j,k,im,4)*dnzav
         rav   = 0.5*(r1 + r0)
         dr    = r1 - r0
         dp    = dr*rho*qthet**2/rav
         qk0(j,im,5,1) = pm1 + dp
         qk0(j,im,5,2) = qk0(j,im,5,1)
         r0  = r1
         pm1 = qk0(j,im,5,1)
  330    continue
      end if
c
      if (lijk .eq. 2) then
c
c        radial integration of dp/dr in j-direction
c        (i is the theta-direction)
c
         if (ldir .gt. 0) then
            j0   = jsta
            js   = jsta + 1
            je   = jend
            jj0  = 1
            jinc = 1
            jjc  = jdimc-1
         else
            j0   = jend
            js   = jend - 1
            je   = jsta
            jj0  = jend - jsta
            jinc = 0
            jjc  = 1
         end if
c
         if (nblc .ne. 0) then
c
c           set data for continuation
c
            if (ichk .gt. 0) then
               do 340 i=ista,iend1
               ii = i-ista+1
               bcdata(jj0,ii,ipp,2) = qk0c(jjc,i,5,1)/p0
  340          continue
            else
               do 350 i=ista,iend1
               ii = i-ista+1
               bcdata(jj0,ii,ipp,2) = q(jj0,k,i,5)/p0
  350          continue
            end if
         end if
c
         do 360 i=ista,iend1
         ii   = i-ista+1
         pm1  = bcdata(jj0,ii,ipp,2)*p0
         j    = j0
         ip   = i  + 1
         kp   = k  + 1
         xav0 = 0.25*(x(j,k,i)  + x(j,kp,i)
     .              + x(j,k,ip) + x(j,kp,ip))
         yav0 = 0.25*(y(j,k,i)  + y(j,kp,i)
     .              + y(j,k,ip) + y(j,kp,ip))
         zav0 = 0.25*(z(j,k,i)  + z(j,kp,i)
     .              + z(j,k,ip) + z(j,kp,ip))
         xav0 = xav0*xfact
         yav0 = yav0*yfact
         zav0 = zav0*zfact
         r0   = sqrt(xav0**2 + yav0**2 + zav0**2)
         do 360 j=js,je,ldir
         jm   = j  - jinc
         kp   = k  + 1
         ip   = i  + 1
         xav1 = 0.25*(x(j,k,i)  + x(j,kp,i)
     .              + x(j,k,ip) + x(j,kp,ip))
         yav1 = 0.25*(y(j,k,i)  + y(j,kp,i)
     .              + y(j,k,ip) + y(j,kp,ip))
         zav1 = 0.25*(z(j,k,i)  + z(j,kp,i)
     .              + z(j,k,ip) + z(j,kp,ip))
         xav1 = xav1*xfact
         yav1 = yav1*yfact
         zav1 = zav1*zfact
         r1   = sqrt(xav1**2 + yav1**2 + zav1**2)
         rho  = q(jm,k,i,1 )
         dnxav = 0.5*(si(jm,k,i,1) + si(jm,k,i+1,1))
         dnyav = 0.5*(si(jm,k,i,2) + si(jm,k,i+1,2))
         dnzav = 0.5*(si(jm,k,i,3) + si(jm,k,i+1,3))
         qthet = q(jm,k,i,2)*dnxav + q(jm,k,i,3)*dnyav
     .         + q(jm,k,i,4)*dnzav
         rav   = 0.5*(r1 + r0)
         dr    = r1 - r0
         dp    = dr*rho*qthet**2/rav
         qk0(jm,i,5,1) = pm1 + dp
         qk0(jm,i,5,2) = qk0(jm,i,5,1)
         r0  = r1
         pm1 = qk0(jm,i,5,1)
  360    continue
      end if
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 370 i=ista,iend1
        do 370 j=jsta,jend1
          vk0(j,i,1,1) = vist3d(j,1,i)
          vk0(j,i,1,2) = vist3d(j,1,i)
  370   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 380 i=ista,iend1
        do 380 j=jsta,jend1
          tk0(j,i,l,1) = tursav(j,1,i,l)
          tk0(j,i,l,2) = tursav(j,1,i,l)
  380   continue
        enddo
      end if
      end if
      end if
c
c******************************************************************************
c      k=kdim boundary   set pressure via radial equilibrium        bctype 2006
c******************************************************************************
c
      if (nface.eq.6) then
c
      k = kdim - 1
c
c     extrapolate all primatives except pressure
c
      do 400 l=1,4
      do 400 i=ista,iend1
      do 400 j=jsta,jend1
      qk0(j,i,l,3) = q(j,kdim1,i,l)
      qk0(j,i,l,4) = qk0(j,i,l,3)
      bck(j,i,2) = 0.0
  400 continue
c
      if (lijk .eq. 1) then
c
c        radial integration of dp/dr in i-direction
c        (j is the theta-direction)
c
         if (ldir .gt. 0) then
            i0   = ista
            is   = ista + 1
            ie   = iend
            ii0  = 1
            iinc = 1
            iic  = idimc-1
         else 
            i0   = iend
            is   = iend - 1
            ie   = ista 
            ii0  = iend - ista
            iinc = 0
            iic  = 1
         end if
c
         if (nblc .ne. 0) then
c
c           set data for continuation
c
            if (ichk .gt. 0) then
               do 410 j=jsta,jend1
               jj = j-jsta+1
               bcdata(jj,ii0,ipp,2) = qk0c(j,iic,5,3)/p0
  410          continue
            else
               do 420 j=jsta,jend1
               jj = j-jsta+1
               bcdata(jj,ii0,ipp,2) = q(j,k,ii0,5)/p0
  420          continue
            end if
         end if
c
         do 430 j=jsta,jend1
         jj   = j-jsta+1
         pm1  = bcdata(jj,ii0,ipp,2)*p0
         i    = i0
         kp   = k  + 1
         jp   = j  + 1
         xav0 = 0.25*(x(j,k,i)  + x(j,kp,i)
     .              + x(jp,k,i) + x(jp,kp,i))
         yav0 = 0.25*(y(j,k,i)  + y(j,kp,i)
     .              + y(jp,k,i) + y(jp,kp,i))
         zav0 = 0.25*(z(j,k,i)  + z(j,kp,i)
     .              + z(jp,k,i) + z(jp,kp,i))
         xav0 = xav0*xfact
         yav0 = yav0*yfact
         zav0 = zav0*zfact
         r0   = sqrt(xav0**2 + yav0**2 + zav0**2)
         do 430 i=is,ie,ldir
         im   = i  - iinc
         kp   = k  + 1
         jp   = j  + 1
         xav1 = 0.25*(x(j,k,i)  + x(j,kp,i)
     .              + x(jp,k,i) + x(jp,kp,i))
         yav1 = 0.25*(y(j,k,i)  + y(j,kp,i)
     .              + y(jp,k,i) + y(jp,kp,i))
         zav1 = 0.25*(z(j,k,i)  + z(j,kp,i)
     .              + z(jp,k,i) + z(jp,kp,i))
         xav1 = xav1*xfact
         yav1 = yav1*yfact
         zav1 = zav1*zfact
         r1   = sqrt(xav1**2 + yav1**2 + zav1**2)
         rho  = q(j,k,im,1 )
         dnxav = 0.5*(sj(j,k,im,1) + sj(j+1,k,im,1))
         dnyav = 0.5*(sj(j,k,im,2) + sj(j+1,k,im,2))
         dnzav = 0.5*(sj(j,k,im,3) + sj(j+1,k,im,3))
         qthet = q(j,k,im,2)*dnxav + q(j,k,im,3)*dnyav
     .         + q(j,k,im,4)*dnzav
         rav   = 0.5*(r1 + r0)
         dr    = r1 - r0
         dp    = dr*rho*qthet**2/rav
         qk0(j,im,5,3) = pm1 + dp
         qk0(j,im,5,4) = qk0(j,im,5,3)
         r0  = r1
         pm1 = qk0(j,im,5,3)
  430    continue
      end if
c
      if (lijk .eq. 2) then
c
c        radial integration of dp/dr in j-direction
c        (i is the theta-direction)
c
         if (ldir .gt. 0) then
            j0   = jsta
            js   = jsta + 1
            je   = jend
            jj0  = 1
            jinc = 1
            jjc  = jdimc-1
         else
            j0   = jend
            js   = jend - 1
            je   = jsta
            jj0  = jend - jsta
            jinc = 0
            jjc  = 1
         end if
c
         if (nblc .ne. 0) then
c
c           set data for continuation
c
            if (ichk .gt. 0) then
               do 440 i=ista,iend1
               ii = i-ista+1
               bcdata(jj0,ii,ipp,2) = qk0c(jjc,i,5,3)/p0
  440          continue
            else
               do 450 i=ista,iend1
               ii = i-ista+1
               bcdata(jj0,ii,ipp,2) = q(jj0,k,i,5)/p0
  450          continue
            end if
         end if
c
         do 460 i=ista,iend1
         ii   = i-ista+1
         pm1  = bcdata(jj0,ii,ipp,2)*p0
         j    = j0
         ip   = i  + 1
         kp   = k  + 1
         xav0 = 0.25*(x(j,k,i)  + x(j,kp,i)
     .              + x(j,k,ip) + x(j,kp,ip))
         yav0 = 0.25*(y(j,k,i)  + y(j,kp,i)
     .              + y(j,k,ip) + y(j,kp,ip))
         zav0 = 0.25*(z(j,k,i)  + z(j,kp,i)
     .              + z(j,k,ip) + z(j,kp,ip))
         xav0 = xav0*xfact
         yav0 = yav0*yfact
         zav0 = zav0*zfact
         r0   = sqrt(xav0**2 + yav0**2 + zav0**2)
         do 460 j=js,je,ldir
         jm   = j  - jinc
         kp   = k  + 1
         ip   = i  + 1
         xav1 = 0.25*(x(j,k,i)  + x(j,kp,i)
     .              + x(j,k,ip) + x(j,kp,ip))
         yav1 = 0.25*(y(j,k,i)  + y(j,kp,i)
     .              + y(j,k,ip) + y(j,kp,ip))
         zav1 = 0.25*(z(j,k,i)  + z(j,kp,i)
     .              + z(j,k,ip) + z(j,kp,ip))
         xav1 = xav1*xfact
         yav1 = yav1*yfact
         zav1 = zav1*zfact
         r1   = sqrt(xav1**2 + yav1**2 + zav1**2)
         rho  = q(jm,k,i,1 )
         dnxav = 0.5*(si(jm,k,i,1) + si(jm,k,i+1,1))
         dnyav = 0.5*(si(jm,k,i,2) + si(jm,k,i+1,2))
         dnzav = 0.5*(si(jm,k,i,3) + si(jm,k,i+1,3))
         qthet = q(jm,k,i,2)*dnxav + q(jm,k,i,3)*dnyav
     .         + q(jm,k,i,4)*dnzav
         rav   = 0.5*(r1 + r0)
         dr    = r1 - r0
         dp    = dr*rho*qthet**2/rav
         qk0(jm,i,5,3) = pm1 + dp
         qk0(jm,i,5,4) = qk0(jm,i,5,3)
         r0  = r1
         pm1 = qk0(jm,i,5,1)
  460    continue
      end if
c
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 470 i=ista,iend1
        do 470 j=jsta,jend1
          vk0(j,i,1,3) = vist3d(j,kdim1,i)
          vk0(j,i,1,4) = vist3d(j,kdim1,i)
  470   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 480 i=ista,iend1
        do 480 j=jsta,jend1
          tk0(j,i,l,3) = tursav(j,kdim1,i,l)
          tk0(j,i,l,4) = tursav(j,kdim1,i,l)
  480   continue
        enddo
      end if
      end if
      end if
c
c******************************************************************************
c      i=1 boundary   set pressure via radial equilibrium           bctype 2006
c******************************************************************************
c
      if (nface.eq.1) then
c
      i = 1
c
c     extrapolate all primatives except pressure
c
      do 500 l=1,4
      do 500 k=ksta,kend1
      do 500 j=jsta,jend1
      qi0(j,k,l,1) = q(j,k,1,l)
      qi0(j,k,l,2) = qi0(j,k,l,1)
      bci(j,k,1) = 0.0
  500 continue
c
      if (lijk .eq. 2) then
c
c        radial integration of dp/dr in j-direction
c        (k is the theta-direction)
c
         if (ldir .gt. 0) then
            j0   = jsta
            js   = jsta + 1
            je   = jend
            jj0  = 1
            jinc = 1
            jjc  = jdimc-1
         else
            j0   = jend
            js   = jend - 1
            je   = jsta
            jj0  = jend - jsta
            jinc = 0
            jjc  = 1
         end if
c
         if (nblc .ne. 0) then
c
c           set data for continuation
c
            if (ichk .gt. 0) then
               do 510 k=ksta,kend1
               kk = k-ksta+1
               bcdata(jj0,kk,ipp,2) = qi0c(jjc,k,5,1)/p0
  510          continue
            else
               do 520 k=ksta,kend1
               kk = k-ksta+1
               bcdata(jj0,kk,ipp,2) = q(jj0,k,i,5)/p0
  520          continue
            end if
         end if
c
         do 530 k=ksta,kend1
         kk   = k-ksta+1
         pm1  = bcdata(jj0,kk,ipp,2)*p0
         j    = j0
         ip   = i  + 1
         kp   = k  + 1
         xav0 = 0.25*(x(j,k,i)  + x(j,kp,i)
     .              + x(j,k,ip) + x(j,kp,ip))
         yav0 = 0.25*(y(j,k,i)  + y(j,kp,i)
     .              + y(j,k,ip) + y(j,kp,ip))
         zav0 = 0.25*(z(j,k,i)  + z(j,kp,i)
     .              + z(j,k,ip) + z(j,kp,ip))
         xav0 = xav0*xfact
         yav0 = yav0*yfact
         zav0 = zav0*zfact
         r0   = sqrt(xav0**2 + yav0**2 + zav0**2)
         do 530 j=js,je,ldir
         jm   = j  - jinc
         kp   = k  + 1
         ip   = i  + 1
         xav1 = 0.25*(x(j,k,i)  + x(j,kp,i)
     .              + x(j,k,ip) + x(j,kp,ip))
         yav1 = 0.25*(y(j,k,i)  + y(j,kp,i)
     .              + y(j,k,ip) + y(j,kp,ip))
         zav1 = 0.25*(z(j,k,i)  + z(j,kp,i)
     .              + z(j,k,ip) + z(j,kp,ip))
         xav1 = xav1*xfact
         yav1 = yav1*yfact
         zav1 = zav1*zfact
         r1   = sqrt(xav1**2 + yav1**2 + zav1**2)
         rho  = q(jm,k,i,1 )
         dnxav = 0.5*(sk(jm,k,i,1) + sk(jm,k+1,i,1))
         dnyav = 0.5*(sk(jm,k,i,2) + sk(jm,k+1,i,2))
         dnzav = 0.5*(sk(jm,k,i,3) + sk(jm,k+1,i,3))
         qthet = q(jm,k,i,2)*dnxav + q(jm,k,i,3)*dnyav
     .         + q(jm,k,i,4)*dnzav
         rav   = 0.5*(r1 + r0)
         dr    = r1 - r0
         dp    = dr*rho*qthet**2/rav
         qi0(jm,k,5,1) = pm1 + dp
         qi0(jm,k,5,2) = qi0(jm,k,5,1)
         r0  = r1
         pm1 = qi0(jm,k,5,3)
  530    continue
      end if
c
      if (lijk .eq. 3) then
c
c        radial integration of dp/dr in k-direction
c        (j is the theta-direction)
c
         if (ldir .gt. 0) then
            k0   = ksta
            ks   = ksta + 1
            ke   = kend
            kk0  = 1
            kinc = 1
            kkc  = kdimc-1
         else
            k0   = kend
            ks   = kend - 1
            ke   = ksta
            kk0  = kend - ksta
            kinc = 0
            kkc  = 1
         end if
c
         if (nblc .ne. 0) then
c
c           set data for continuation
c
            if (ichk .gt. 0) then
               do 540 j=jsta,jend1
               jj = j-jsta+1
               bcdata(jj,kk0,ipp,2) = qi0c(j,kkc,5,1)/p0
  540          continue
            else
               do 550 j=jsta,jend1
               jj = j-jsta+1
               bcdata(jj,kk0,ipp,2) = q(j,kk0,i,5)/p0
  550          continue
            end if
         end if
c
         do 560 j=jsta,jend1
         jj   = j-jsta+1
         pm1  = bcdata(jj,kk0,ipp,2)*p0
         k    = k0
         ip   = i  + 1
         jp   = j  + 1
         xav0 = 0.25*(x(j,k,i)  + x(j,k,ip)
     .              + x(jp,k,i) + x(jp,k,ip))
         yav0 = 0.25*(y(j,k,i)  + y(j,k,ip)
     .              + y(jp,k,i) + y(jp,k,ip))
         zav0 = 0.25*(z(j,k,i)  + z(j,k,ip)
     .              + z(jp,k,i) + z(jp,k,ip))
         xav0 = xav0*xfact
         yav0 = yav0*yfact
         zav0 = zav0*zfact
         r0   = sqrt(xav0**2 + yav0**2 + zav0**2)
         do 560 k=ks,ke,ldir
         km   = k - kinc
         ip   = i  + 1
         jp   = j  + 1
         xav1 = 0.25*(x(j,k,i)  + x(j,k,ip)
     .              + x(jp,k,i) + x(jp,k,ip))
         yav1 = 0.25*(y(j,k,i)  + y(j,k,ip)
     .              + y(jp,k,i) + y(jp,k,ip))
         zav1 = 0.25*(z(j,k,i)  + z(j,k,ip)
     .              + z(jp,k,i) + z(jp,k,ip))
         xav1 = xav1*xfact
         yav1 = yav1*yfact
         zav1 = zav1*zfact
         r1   = sqrt(xav1**2 + yav1**2 + zav1**2)
         rho  = q(j,km,i,1 )
         dnxav = 0.5*(sj(j,km,i,1) + sj(j+1,km,i,1))
         dnyav = 0.5*(sj(j,km,i,2) + sj(j+1,km,i,2))
         dnzav = 0.5*(sj(j,km,i,3) + sj(j+1,km,i,3))
         qthet = q(j,km,i,2)*dnxav + q(j,km,i,3)*dnyav
     .         + q(j,km,i,4)*dnzav
         rav   = 0.5*(r1 + r0)
         dr    = r1 - r0
         dp    = dr*rho*qthet**2/rav
         qi0(j,km,5,1) = pm1 + dp
         qi0(j,km,5,2) = qi0(j,km,5,1)
         r0  = r1
         pm1 = qi0(j,km,5,1)
  560    continue
      end if
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 570 k=ksta,kend1
        do 570 j=jsta,jend1
          vi0(j,k,1,1) = vist3d(j,k,1)
          vi0(j,k,1,2) = vist3d(j,k,1)
  570   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 580 k=ksta,kend1
        do 580 j=jsta,jend1
          ti0(j,k,l,1) = tursav(j,k,1,l)
          ti0(j,k,l,2) = tursav(j,k,1,l)
  580   continue
        enddo
      end if
      end if
      end if
c
c******************************************************************************
c      i=idim boundary   set pressure via radial equilibrium        bctype 2006
c******************************************************************************
c
      if (nface.eq.2) then
c
      i = idim - 1
c
c     extrapolate all primatives except pressure
c
      do 600 l=1,4
      do 600 k=ksta,kend1
      do 600 j=jsta,jend1
      qi0(j,k,l,3) = q(j,k,idim1,l)
      qi0(j,k,l,4) = qi0(j,k,l,3)
      bci(j,k,2) = 0.0
  600 continue
c
      if (lijk .eq. 2) then
c
c        radial integration of dp/dr in j-direction
c        (k is theta-direction)
c
         if (ldir .gt. 0) then
            j0   = jsta
            js   = jsta + 1
            je   = jend
            jj0  = 1
            jinc = 1
            jjc  = jdimc-1
         else
            j0   = jend
            js   = jend - 1
            je   = jsta
            jj0  = jend - jsta
            jinc = 0
            jjc  = 1
         end if
c
         if (nblc .ne. 0) then
c
c           set data for continuation
c
            if (ichk .gt. 0) then
               do 610 k=ksta,kend1
               kk = k-ksta+1
               bcdata(jj0,kk,ipp,2) = qi0c(jjc,k,5,3)/p0
  610          continue
            else
               do 620 k=ksta,kend1
               kk = k-ksta+1
               bcdata(jj0,kk,ipp,2) = q(jj0,k,i,5)/p0
  620          continue
            end if
         end if
c
         do 630 k=ksta,kend1
         kk   = k-ksta+1
         pm1  = bcdata(jj0,kk,ipp,2)*p0
         j    = j0
         ip   = i  + 1
         kp   = k  + 1
         xav0 = 0.25*(x(j,k,i)  + x(j,kp,i)
     .              + x(j,k,ip) + x(j,kp,ip))
         yav0 = 0.25*(y(j,k,i)  + y(j,kp,i)
     .              + y(j,k,ip) + y(j,kp,ip))
         zav0 = 0.25*(z(j,k,i)  + z(j,kp,i)
     .              + z(j,k,ip) + z(j,kp,ip))
         xav0 = xav0*xfact
         yav0 = yav0*yfact
         zav0 = zav0*zfact
         r0   = sqrt(xav0**2 + yav0**2 + zav0**2)
         do 630 j=js,je,ldir
         jm   = j  - jinc
         kp   = k  + 1
         ip   = i  + 1
         xav1 = 0.25*(x(j,k,i)  + x(j,kp,i)
     .              + x(j,k,ip) + x(j,kp,ip))
         yav1 = 0.25*(y(j,k,i)  + y(j,kp,i)
     .              + y(j,k,ip) + y(j,kp,ip))
         zav1 = 0.25*(z(j,k,i)  + z(j,kp,i)
     .              + z(j,k,ip) + z(j,kp,ip))
         xav1 = xav1*xfact
         yav1 = yav1*yfact
         zav1 = zav1*zfact
         r1   = sqrt(xav1**2 + yav1**2 + zav1**2)
         rho  = q(jm,k,i,1 )
         dnxav = 0.5*(sk(jm,k,i,1) + sk(jm,k+1,i,1))
         dnyav = 0.5*(sk(jm,k,i,2) + sk(jm,k+1,i,2))
         dnzav = 0.5*(sk(jm,k,i,3) + sk(jm,k+1,i,3))
         qthet = q(jm,k,i,2)*dnxav + q(jm,k,i,3)*dnyav
     .         + q(jm,k,i,4)*dnzav
         rav   = 0.5*(r1 + r0)
         dr    = r1 - r0
         dp    = dr*rho*qthet**2/rav
         qi0(jm,k,5,3) = pm1 + dp
         qi0(jm,k,5,4) = qi0(jm,k,5,3)
         r0  = r1
         pm1 = qi0(jm,k,5,3)
  630    continue
      end if
c
      if (lijk .eq. 3) then
c
c        radial integration of dp/dr in k-direction
c        (j is the theta-direction)
c
         if (ldir .gt. 0) then
            k0   = ksta
            ks   = ksta + 1
            ke   = kend
            kk0  = 1
            kinc = 1
            kkc  = kdimc-1
         else
            k0   = kend
            ks   = kend - 1
            ke   = ksta
            kk0  = kend - ksta
            kinc = 0
            kkc  = 1
         end if
c
         if (nblc .ne. 0) then
c
c           set data for continuation
c
            if (ichk .gt. 0) then
               do 640 j=jsta,jend1
               jj = j-jsta+1
               bcdata(jj,kk0,ipp,2) = qi0c(j,kkc,5,3)/p0
  640          continue
            else
               do 650 j=jsta,jend1
               jj = j-jsta+1
               bcdata(jj,kk0,ipp,2) = q(j,kk0,i,5)/p0
  650          continue
            end if
         end if
c
         do 660 j=jsta,jend1
         jj   = j-jsta+1
         pm1  = bcdata(jj,kk0,ipp,2)*p0
         k    = k0
         ip   = i  + 1
         jp   = j  + 1
         xav0 = 0.25*(x(j,k,i)  + x(j,k,ip)
     .              + x(jp,k,i) + x(jp,k,ip))
         yav0 = 0.25*(y(j,k,i)  + y(j,k,ip)
     .              + y(jp,k,i) + y(jp,k,ip))
         zav0 = 0.25*(z(j,k,i)  + z(j,k,ip)
     .              + z(jp,k,i) + z(jp,k,ip))
         xav0 = xav0*xfact
         yav0 = yav0*yfact
         zav0 = zav0*zfact
         r0   = sqrt(xav0**2 + yav0**2 + zav0**2)
         do 660 k=ks,ke,ldir
         km   = k - kinc
         ip   = i  + 1
         jp   = j  + 1
         xav1 = 0.25*(x(j,k,i)  + x(j,k,ip)
     .              + x(jp,k,i) + x(jp,k,ip))
         yav1 = 0.25*(y(j,k,i)  + y(j,k,ip)
     .              + y(jp,k,i) + y(jp,k,ip))
         zav1 = 0.25*(z(j,k,i)  + z(j,k,ip)
     .              + z(jp,k,i) + z(jp,k,ip))
         xav1 = xav1*xfact
         yav1 = yav1*yfact
         zav1 = zav1*zfact
         r1   = sqrt(xav1**2 + yav1**2 + zav1**2)
         rho  = q(j,km,i,1 )
         dnxav = 0.5*(sj(j,km,i,1) + sj(j+1,km,i,1))
         dnyav = 0.5*(sj(j,km,i,2) + sj(j+1,km,i,2))
         dnzav = 0.5*(sj(j,km,i,3) + sj(j+1,km,i,3))
         qthet = q(j,km,i,2)*dnxav + q(j,km,i,3)*dnyav 
     .         + q(j,km,i,4)*dnzav
         rav   = 0.5*(r1 + r0)
         dr    = r1 - r0
         dp    = dr*rho*qthet**2/rav
         qi0(j,km,5,3) = pm1 + dp
         qi0(j,km,5,4) = qi0(j,km,5,3)
         r0  = r1
         pm1 = qi0(j,km,5,3)
  660    continue
      end if
c
      if (ivisc(3).ge.2 .or. ivisc(2).ge.2 .or. ivisc(1).ge.2) then
        do 670 k=ksta,kend1
        do 670 j=jsta,jend1
          vi0(j,k,1,3) = vist3d(j,k,idim1)
          vi0(j,k,1,4) = vist3d(j,k,idim1)
  670   continue
      end if
c   only need to do advanced model turbulence B.C.s on finest grid
      if (level .ge. lglobal) then
      if (ivisc(3).ge.4 .or. ivisc(2).ge.4 .or. ivisc(1).ge.4) then
        do l=1,nummem
        do 680 k=ksta,kend1
        do 680 j=jsta,jend1
          ti0(j,k,l,3) = tursav(j,k,idim1,l)
          ti0(j,k,l,4) = tursav(j,k,idim1,l)
  680   continue
        enddo
      end if
      end if
      end if
c
      return 
      end
