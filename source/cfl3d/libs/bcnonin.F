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
      subroutine bcnonin(jdim,kdim,idim,q,qj0,qk0,qi0,sj,sk,si,bcj,bck,
     .                  bci,xtbj,xtbk,xtbi,atbj,atbk,atbi,ista,iend,
     .                  jsta,jend,ksta,kend,nface,
     .                  iuns,nou,bou,nbuf,
     .                  ibufdim,x,y,z,nbl)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Add noninertial rotating component to freestream 
c               boundary conditions 
c               Original coding by Mike Park
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension q(jdim,kdim,idim,5), qi0(jdim,kdim,5,4),
     .          qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4)
      dimension sk(jdim,kdim,idim-1,5),si(jdim,kdim,idim,5),
     .          sj(jdim,kdim,idim-1,5)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
      dimension xtbj(kdim,idim-1,3,2),xtbk(jdim,idim-1,3,2),
     .          xtbi(jdim,kdim,3,2),atbj(kdim,idim-1,3,2),
     .          atbk(jdim,idim-1,3,2),atbi(jdim,kdim,3,2)
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
 
      dimension dx(2),dy(2),dz(2)
c
      common /sklton/ isklton
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /noninertial/ xcentrot,ycentrot,zcentrot,xrotrate,
     .                     yrotrate,zrotrate,noninflag
c
      wx = xrotrate
      wy = yrotrate
      wz = zrotrate
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
      jend1 = jend-1
      kend1 = kend-1
      iend1 = iend-1
c
c            * * * * * * * * * * * * * * * * * * * * * * *
c            * NON-INERT boundary condition bctype=1000  *
c            * * * * * * * * * * * * * * * * * * * * * * *
 
c     To find the noninertial components to the boundary conditions the 
c     the cell centers of the ghost cells must be approximated. This is 
c     done by computing the cell center(cx, cy, cz) of the neighboring 
c     computational cell(NCC) and a vector (dx, dy, dz) between the NCC 
c     and the next interior point. This vector is used to extrapolate 
c     from the NCC to the cell center of the two ghost cells(gx,gy,gz).
 
c     Tue Aug  3 11:38:08 EDT 1999
c     improved the extrapolation for better constancy in a highly 
c     stretched and/or curved mesh. dx(1), dy(1), and dz(1) are the 
c     vectors from the cell center of the border NCC to the border 
c     (nearest) ghost cell. dx(2), dy(2), and dz(2) are vectors 
c     from the border NCC to the second ghost cell.
c   
c     dx(1), dy(1), and dz(1) is the average length and direction 
c     of the sides of the NCC that touch the boundary face
c   
c     dx(2), dy(2), and dz(2) is dx(1), dy(1), and dz(1) plus
c     a vector from the second NCC cell center to the border NCC cell 
c     center
 
c
c******************************************************************************
c      j=1 boundary NON-INERTIAL freestream                         bctype 1000
c******************************************************************************
c
      if (nface.eq.3) then
      if (isklton.eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1001) ista,iend,ksta,kend
      end if
      do 100 i=ista,iend1
      do 100 k=ksta,kend1
 
      j=1

       cx = 0.125 * ( 
     . x(j  , k  , i  ) + x(j  , k  , i+1) +
     . x(j  , k+1, i  ) + x(j  , k+1, i+1) +
     . x(j+1, k  , i  ) + x(j+1, k  , i+1) +
     . x(j+1, k+1, i  ) + x(j+1, k+1, i+1) )
 
      cy = 0.125 * ( 
     . y(j  , k  , i  ) + y(j  , k  , i+1) +
     . y(j  , k+1, i  ) + y(j  , k+1, i+1) +
     . y(j+1, k  , i  ) + y(j+1, k  , i+1) +
     . y(j+1, k+1, i  ) + y(j+1, k+1, i+1) )
 
      cz = 0.125 * ( 
     . z(j  , k  , i  ) + z(j  , k  , i+1) +
     . z(j  , k+1, i  ) + z(j  , k+1, i+1) +
     . z(j+1, k  , i  ) + z(j+1, k  , i+1) +
     . z(j+1, k+1, i  ) + z(j+1, k+1, i+1) )
 
      dx(1) = 0.25 * (
     . x(j  , k  , i  ) - x(j+1, k  , i  ) +
     . x(j  , k  , i+1) - x(j+1, k  , i+1) +
     . x(j  , k+1, i  ) - x(j+1, k+1, i  ) +
     . x(j  , k+1, i+1) - x(j+1, k+1, i+1) )
 
      dy(1) = 0.25 * (
     . y(j  , k  , i  ) - y(j+1, k  , i  ) +
     . y(j  , k  , i+1) - y(j+1, k  , i+1) +
     . y(j  , k+1, i  ) - y(j+1, k+1, i  ) +
     . y(j  , k+1, i+1) - y(j+1, k+1, i+1) )
 
      dz(1) = 0.25 * (
     . z(j  , k  , i  ) - z(j+1, k  , i  ) +
     . z(j  , k  , i+1) - z(j+1, k  , i+1) +
     . z(j  , k+1, i  ) - z(j+1, k+1, i  ) +
     . z(j  , k+1, i+1) - z(j+1, k+1, i+1) )
 
      j=2
 
      dx(2) = dx(1) + cx - 0.125 * ( 
     . x(j  , k  , i  ) + x(j  , k  , i+1) +
     . x(j  , k+1, i  ) + x(j  , k+1, i+1) +
     . x(j+1, k  , i  ) + x(j+1, k  , i+1) +
     . x(j+1, k+1, i  ) + x(j+1, k+1, i+1) )
 
      dy(2) = dy(1) + cy - 0.125 * ( 
     . y(j  , k  , i  ) + y(j  , k  , i+1) +
     . y(j  , k+1, i  ) + y(j  , k+1, i+1) +
     . y(j+1, k  , i  ) + y(j+1, k  , i+1) +
     . y(j+1, k+1, i  ) + y(j+1, k+1, i+1) )
 
      dz(2) = dz(1) + cz - 0.125 * ( 
     . z(j  , k  , i  ) + z(j  , k  , i+1) +
     . z(j  , k+1, i  ) + z(j  , k+1, i+1) +
     . z(j+1, k  , i  ) + z(j+1, k  , i+1) +
     . z(j+1, k+1, i  ) + z(j+1, k+1, i+1) )
 
      cx = cx - xcentrot
      cy = cy - ycentrot
      cz = cz - zcentrot
 
      do 100 l=1,2
 
      gx = cx + dx(l)
      gy = cy + dy(l)
      gz = cz + dz(l)
 
c     Uinf + r x omega is the same as Uinf - omega x r 
 
      qj0(k,i,2,l) = qj0(k,i,2,l) + ( gy * wz - gz * wy )
      qj0(k,i,3,l) = qj0(k,i,3,l) + ( gz * wx - gx * wz )
      qj0(k,i,4,l) = qj0(k,i,4,l) + ( gx * wy - gy * wx )
 
      bcj(k,i,1)   = 0.0
 
 910  format(a,4f7.2,8f11.5)
 911  format(8f11.5)
 
  100 continue
      end if
c
c******************************************************************************
c      j=jdim boundary NON-INERTIAL freestream                      bctype 1000
c******************************************************************************
c
      if (nface.eq.4) then
      if (isklton.eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1002) ista,iend,ksta,kend
      end if
      do 200 i=ista,iend1
      do 200 k=ksta,kend1
 
      j=jdim1
 
      cx = 0.125 * ( 
     . x(j  , k  , i  ) + x(j  , k  , i+1) +
     . x(j  , k+1, i  ) + x(j  , k+1, i+1) +
     . x(j+1, k  , i  ) + x(j+1, k  , i+1) +
     . x(j+1, k+1, i  ) + x(j+1, k+1, i+1) )
 
      cy = 0.125 * ( 
     . y(j  , k  , i  ) + y(j  , k  , i+1) +
     . y(j  , k+1, i  ) + y(j  , k+1, i+1) +
     . y(j+1, k  , i  ) + y(j+1, k  , i+1) +
     . y(j+1, k+1, i  ) + y(j+1, k+1, i+1) )
 
      cz = 0.125 * ( 
     . z(j  , k  , i  ) + z(j  , k  , i+1) +
     . z(j  , k+1, i  ) + z(j  , k+1, i+1) +
     . z(j+1, k  , i  ) + z(j+1, k  , i+1) +
     . z(j+1, k+1, i  ) + z(j+1, k+1, i+1) )
 
      dx(1) = 0.25 * (
     . x(j+1, k  , i  ) - x(j  , k  , i  ) +
     . x(j+1, k  , i+1) - x(j  , k  , i+1) +
     . x(j+1, k+1, i  ) - x(j  , k+1, i  ) +
     . x(j+1, k+1, i+1) - x(j  , k+1, i+1) )
 
      dy(1) = 0.25 * (
     . y(j+1, k  , i  ) - y(j  , k  , i  ) +
     . y(j+1, k  , i+1) - y(j  , k  , i+1) +
     . y(j+1, k+1, i  ) - y(j  , k+1, i  ) +
     . y(j+1, k+1, i+1) - y(j  , k+1, i+1) )
 
      dz(1) = 0.25 * (
     . z(j+1, k  , i  ) - z(j  , k  , i  ) +
     . z(j+1, k  , i+1) - z(j  , k  , i+1) +
     . z(j+1, k+1, i  ) - z(j  , k+1, i  ) +
     . z(j+1, k+1, i+1) - z(j  , k+1, i+1) )
 
      j=jdim1-1
 
      dx(2) = dx(1) + cx - 0.125 * ( 
     . x(j  , k  , i  ) + x(j  , k  , i+1) +
     . x(j  , k+1, i  ) + x(j  , k+1, i+1) +
     . x(j+1, k  , i  ) + x(j+1, k  , i+1) +
     . x(j+1, k+1, i  ) + x(j+1, k+1, i+1) )
 
      dy(2) = dy(1) + cy - 0.125 * ( 
     . y(j  , k  , i  ) + y(j  , k  , i+1) +
     . y(j  , k+1, i  ) + y(j  , k+1, i+1) +
     . y(j+1, k  , i  ) + y(j+1, k  , i+1) +
     . y(j+1, k+1, i  ) + y(j+1, k+1, i+1) )
 
      dz(2) = dz(1) + cz - 0.125 * ( 
     . z(j  , k  , i  ) + z(j  , k  , i+1) +
     . z(j  , k+1, i  ) + z(j  , k+1, i+1) +
     . z(j+1, k  , i  ) + z(j+1, k  , i+1) +
     . z(j+1, k+1, i  ) + z(j+1, k+1, i+1) )
 
      cx = cx - xcentrot
      cy = cy - ycentrot
      cz = cz - zcentrot
 
      do 200 l=3,4
 
      gx = cx + dx(l-2)
      gy = cy + dy(l-2)
      gz = cz + dz(l-2)
 
c     Uinf + r x omega is the same as Uinf - omega x r 
 
      qj0(k,i,2,l) = qj0(k,i,2,l) + ( gy * wz - gz * wy )
      qj0(k,i,3,l) = qj0(k,i,3,l) + ( gz * wx - gx * wz )
      qj0(k,i,4,l) = qj0(k,i,4,l) + ( gx * wy - gy * wx ) 
 
      bcj(k,i,2)   = 0.0
 
  200 continue
      end if
c
c******************************************************************************
c      k=1 boundary NON-INERTIAL freestream                         bctype 1000
c******************************************************************************
c
      if (nface.eq.5) then
      if (isklton.eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1003) ista,iend,ksta,kend
      end if
      do 300 i=ista,iend1
      do 300 j=jsta,jend1
 
      k=1
 
      cx = 0.125 * ( 
     . x(j  , k  , i  ) + x(j  , k  , i+1) +
     . x(j  , k+1, i  ) + x(j  , k+1, i+1) +
     . x(j+1, k  , i  ) + x(j+1, k  , i+1) +
     . x(j+1, k+1, i  ) + x(j+1, k+1, i+1) )
 
      cy = 0.125 * ( 
     . y(j  , k  , i  ) + y(j  , k  , i+1) +
     . y(j  , k+1, i  ) + y(j  , k+1, i+1) +
     . y(j+1, k  , i  ) + y(j+1, k  , i+1) +
     . y(j+1, k+1, i  ) + y(j+1, k+1, i+1) )
 
      cz = 0.125 * ( 
     . z(j  , k  , i  ) + z(j  , k  , i+1) +
     . z(j  , k+1, i  ) + z(j  , k+1, i+1) +
     . z(j+1, k  , i  ) + z(j+1, k  , i+1) +
     . z(j+1, k+1, i  ) + z(j+1, k+1, i+1) )
 
      dx(1) = 0.25 * (
     . x(j  , k  , i  ) - x(j  , k+1, i  ) +
     . x(j  , k  , i+1) - x(j  , k+1, i+1) +
     . x(j+1, k  , i  ) - x(j+1, k+1, i  ) +
     . x(j+1, k  , i+1) - x(j+1, k+1, i+1) )
 
      dy(1) = 0.25 * (
     . y(j  , k  , i  ) - y(j  , k+1, i  ) +
     . y(j  , k  , i+1) - y(j  , k+1, i+1) +
     . y(j+1, k  , i  ) - y(j+1, k+1, i  ) +
     . y(j+1, k  , i+1) - y(j+1, k+1, i+1) )
 
      dz(1) = 0.25 * (
     . z(j  , k  , i  ) - z(j  , k+1, i  ) +
     . z(j  , k  , i+1) - z(j  , k+1, i+1) +
     . z(j+1, k  , i  ) - z(j+1, k+1, i  ) +
     . z(j+1, k  , i+1) - z(j+1, k+1, i+1) )
 
      k=2
 
      dx(2) = dx(1) + cx - 0.125 * ( 
     . x(j  , k  , i  ) + x(j  , k  , i+1) +
     . x(j  , k+1, i  ) + x(j  , k+1, i+1) +
     . x(j+1, k  , i  ) + x(j+1, k  , i+1) +
     . x(j+1, k+1, i  ) + x(j+1, k+1, i+1) )
 
      dy(2) = dy(1) + cy - 0.125 * ( 
     . y(j  , k  , i  ) + y(j  , k  , i+1) +
     . y(j  , k+1, i  ) + y(j  , k+1, i+1) +
     . y(j+1, k  , i  ) + y(j+1, k  , i+1) +
     . y(j+1, k+1, i  ) + y(j+1, k+1, i+1) )
 
      dz(2) = dz(1) + cz - 0.125 * ( 
     . z(j  , k  , i  ) + z(j  , k  , i+1) +
     . z(j  , k+1, i  ) + z(j  , k+1, i+1) +
     . z(j+1, k  , i  ) + z(j+1, k  , i+1) +
     . z(j+1, k+1, i  ) + z(j+1, k+1, i+1) )
 
      cx = cx - xcentrot
      cy = cy - ycentrot
      cz = cz - zcentrot
 
      do 300 l=1,2
 
      gx = cx + dx(l)
      gy = cy + dy(l)
      gz = cz + dz(l)
 
c     Uinf + r x omega is the same as Uinf - omega x r 
 
      qk0(j,i,2,l) = qk0(j,i,2,l) + ( gy * wz - gz * wy )
      qk0(j,i,3,l) = qk0(j,i,3,l) + ( gz * wx - gx * wz )
      qk0(j,i,4,l) = qk0(j,i,4,l) + ( gx * wy - gy * wx )
      bck(j,i,1)   = 0.0
 
  300 continue
      end if
c
c******************************************************************************
c      k=kdim boundary NON-INERTIAL freestream                      bctype 1000
c******************************************************************************
c
      if (nface.eq.6) then
      if (isklton.eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1004) ista,iend,ksta,kend
      end if
      do 400 i=ista,iend1
      do 400 j=jsta,jend1
 
      k=kdim1
 
      cx = 0.125 * ( 
     . x(j  , k  , i  ) + x(j  , k  , i+1) +
     . x(j  , k+1, i  ) + x(j  , k+1, i+1) +
     . x(j+1, k  , i  ) + x(j+1, k  , i+1) +
     . x(j+1, k+1, i  ) + x(j+1, k+1, i+1) )
 
      cy = 0.125 * ( 
     . y(j  , k  , i  ) + y(j  , k  , i+1) +
     . y(j  , k+1, i  ) + y(j  , k+1, i+1) +
     . y(j+1, k  , i  ) + y(j+1, k  , i+1) +
     . y(j+1, k+1, i  ) + y(j+1, k+1, i+1) )
 
      cz = 0.125 * ( 
     . z(j  , k  , i  ) + z(j  , k  , i+1) +
     . z(j  , k+1, i  ) + z(j  , k+1, i+1) +
     . z(j+1, k  , i  ) + z(j+1, k  , i+1) +
     . z(j+1, k+1, i  ) + z(j+1, k+1, i+1) )
 
      dx(1) = 0.25 * (
     . x(j  , k+1, i  ) - x(j  , k  , i  ) +
     . x(j  , k+1, i+1) - x(j  , k  , i+1) +
     . x(j+1, k+1, i  ) - x(j+1, k  , i  ) +
     . x(j+1, k+1, i+1) - x(j+1, k  , i+1) )
 
      dy(1) = 0.25 * (
     . y(j  , k+1, i  ) - y(j  , k  , i  ) +
     . y(j  , k+1, i+1) - y(j  , k  , i+1) +
     . y(j+1, k+1, i  ) - y(j+1, k  , i  ) +
     . y(j+1, k+1, i+1) - y(j+1, k  , i+1) )
 
      dz(1) = 0.25 * (
     . z(j  , k+1, i  ) - z(j  , k  , i  ) +
     . z(j  , k+1, i+1) - z(j  , k  , i+1) +
     . z(j+1, k+1, i  ) - z(j+1, k  , i  ) +
     . z(j+1, k+1, i+1) - z(j+1, k  , i+1) )
 
      k=kdim1-1
 
      dx(2) = dx(1) + cx - 0.125 * ( 
     . x(j  , k  , i  ) + x(j  , k  , i+1) +
     . x(j  , k+1, i  ) + x(j  , k+1, i+1) +
     . x(j+1, k  , i  ) + x(j+1, k  , i+1) +
     . x(j+1, k+1, i  ) + x(j+1, k+1, i+1) )
 
      dy(2) = dy(1) + cy - 0.125 * ( 
     . y(j  , k  , i  ) + y(j  , k  , i+1) +
     . y(j  , k+1, i  ) + y(j  , k+1, i+1) +
     . y(j+1, k  , i  ) + y(j+1, k  , i+1) +
     . y(j+1, k+1, i  ) + y(j+1, k+1, i+1) )
 
      dz(2) = dz(1) + cz - 0.125 * ( 
     . z(j  , k  , i  ) + z(j  , k  , i+1) +
     . z(j  , k+1, i  ) + z(j  , k+1, i+1) +
     . z(j+1, k  , i  ) + z(j+1, k  , i+1) +
     . z(j+1, k+1, i  ) + z(j+1, k+1, i+1) )
 
      cx = cx - xcentrot
      cy = cy - ycentrot
      cz = cz - zcentrot
 
      do 400 l=3,4
 
      gx = cx + dx(l-2)
      gy = cy + dy(l-2)
      gz = cz + dz(l-2)
 
c     Uinf + r x omega is the same as Uinf - omega x r 
 
      qk0(j,i,2,l) = qk0(j,i,2,l) + ( gy * wz - gz * wy )
      qk0(j,i,3,l) = qk0(j,i,3,l) + ( gz * wx - gx * wz )
      qk0(j,i,4,l) = qk0(j,i,4,l) + ( gx * wy - gy * wx )
 
      bck(j,i,2)   = 0.0
 
  400 continue
      end if
c
c******************************************************************************
c      i=1 boundary NON-INERTIAL freestream                         bctype 1000
c******************************************************************************
c
      if (nface.eq.1) then
      if (isklton.eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1005) ista,iend,ksta,kend
      end if
      do 500 k=ksta,kend1
      do 500 j=jsta,jend1
 
      i=1
 
      cx = 0.125 * ( 
     . x(j  , k  , i  ) + x(j  , k  , i+1) +
     . x(j  , k+1, i  ) + x(j  , k+1, i+1) +
     . x(j+1, k  , i  ) + x(j+1, k  , i+1) +
     . x(j+1, k+1, i  ) + x(j+1, k+1, i+1) )
 
      cy = 0.125 * ( 
     . y(j  , k  , i  ) + y(j  , k  , i+1) +
     . y(j  , k+1, i  ) + y(j  , k+1, i+1) +
     . y(j+1, k  , i  ) + y(j+1, k  , i+1) +
     . y(j+1, k+1, i  ) + y(j+1, k+1, i+1) )
 
      cz = 0.125 * ( 
     . z(j  , k  , i  ) + z(j  , k  , i+1) +
     . z(j  , k+1, i  ) + z(j  , k+1, i+1) +
     . z(j+1, k  , i  ) + z(j+1, k  , i+1) +
     . z(j+1, k+1, i  ) + z(j+1, k+1, i+1) )
 
      dx(1) = 0.25 * (
     . x(j  , k  , i  ) - x(j  , k  , i+1) +
     . x(j  , k+1, i  ) - x(j  , k+1, i+1) +
     . x(j+1, k  , i  ) - x(j+1, k  , i+1) +
     . x(j+1, k+1, i  ) - x(j+1, k+1, i+1) )
 
      dy(1) = 0.25 * (
     . y(j  , k  , i  ) - y(j  , k  , i+1) +
     . y(j  , k+1, i  ) - y(j  , k+1, i+1) +
     . y(j+1, k  , i  ) - y(j+1, k  , i+1) +
     . y(j+1, k+1, i  ) - y(j+1, k+1, i+1) )
 
      dz(1) = 0.25 * (
     . z(j  , k  , i  ) - z(j  , k  , i+1) +
     . z(j  , k+1, i  ) - z(j  , k+1, i+1) +
     . z(j+1, k  , i  ) - z(j+1, k  , i+1) +
     . z(j+1, k+1, i  ) - z(j+1, k+1, i+1) )
 
      i=2
 
      dx(2) = dx(1) + cx - 0.125 * ( 
     . x(j  , k  , i  ) + x(j  , k  , i+1) +
     . x(j  , k+1, i  ) + x(j  , k+1, i+1) +
     . x(j+1, k  , i  ) + x(j+1, k  , i+1) +
     . x(j+1, k+1, i  ) + x(j+1, k+1, i+1) )
 
      dy(2) = dy(1) + cy - 0.125 * ( 
     . y(j  , k  , i  ) + y(j  , k  , i+1) +
     . y(j  , k+1, i  ) + y(j  , k+1, i+1) +
     . y(j+1, k  , i  ) + y(j+1, k  , i+1) +
     . y(j+1, k+1, i  ) + y(j+1, k+1, i+1) )
 
      dz(2) = dz(1) + cz - 0.125 * ( 
     . z(j  , k  , i  ) + z(j  , k  , i+1) +
     . z(j  , k+1, i  ) + z(j  , k+1, i+1) +
     . z(j+1, k  , i  ) + z(j+1, k  , i+1) +
     . z(j+1, k+1, i  ) + z(j+1, k+1, i+1) )
 
      cx = cx - xcentrot
      cy = cy - ycentrot
      cz = cz - zcentrot
 
      do 500 l=1,2
 
      gx = cx + dx(l)
      gy = cy + dy(l)
      gz = cz + dz(l)
 
c     Uinf + r x omega is the same as Uinf - omega x r 
 
      qi0(j,k,2,l) = qi0(j,k,2,l) + ( gy * wz - gz * wy )
      qi0(j,k,3,l) = qi0(j,k,3,l) + ( gz * wx - gx * wz )
      qi0(j,k,4,l) = qi0(j,k,4,l) + ( gx * wy - gy * wx )
      bci(j,i,1)   = 0.0
 
  500 continue
      end if
c
c******************************************************************************
c      i=idim boundary NON-INERTIAL freestream                      bctype 1000
c******************************************************************************
c
      if (nface.eq.2) then
      if (isklton.eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1006) ista,iend,ksta,kend
      end if
      do 600 k=ksta,kend1
      do 600 j=jsta,jend1
 
      i=idim1
 
      cx = 0.125 * ( 
     . x(j  , k  , i  ) + x(j  , k  , i+1) +
     . x(j  , k+1, i  ) + x(j  , k+1, i+1) +
     . x(j+1, k  , i  ) + x(j+1, k  , i+1) +
     . x(j+1, k+1, i  ) + x(j+1, k+1, i+1) )
 
      cy = 0.125 * ( 
     . y(j  , k  , i  ) + y(j  , k  , i+1) +
     . y(j  , k+1, i  ) + y(j  , k+1, i+1) +
     . y(j+1, k  , i  ) + y(j+1, k  , i+1) +
     . y(j+1, k+1, i  ) + y(j+1, k+1, i+1) )
 
      cz = 0.125 * ( 
     . z(j  , k  , i  ) + z(j  , k  , i+1) +
     . z(j  , k+1, i  ) + z(j  , k+1, i+1) +
     . z(j+1, k  , i  ) + z(j+1, k  , i+1) +
     . z(j+1, k+1, i  ) + z(j+1, k+1, i+1) )
 
      dx(1) = 0.25 * (
     . x(j  , k  , i+1) - x(j  , k  , i  ) +
     . x(j  , k+1, i+1) - x(j  , k+1, i  ) +
     . x(j+1, k  , i+1) - x(j+1, k  , i  ) +
     . x(j+1, k+1, i+1) - x(j+1, k+1, i  ) )
 
      dy(1) = 0.25 * (
     . y(j  , k  , i+1) - y(j  , k  , i  ) +
     . y(j  , k+1, i+1) - y(j  , k+1, i  ) +
     . y(j+1, k  , i+1) - y(j+1, k  , i  ) +
     . y(j+1, k+1, i+1) - y(j+1, k+1, i  ) )
 
      dz(1) = 0.25 * (
     . z(j  , k  , i+1) - z(j  , k  , i  ) +
     . z(j  , k+1, i+1) - z(j  , k+1, i  ) +
     . z(j+1, k  , i+1) - z(j+1, k  , i  ) +
     . z(j+1, k+1, i+1) - z(j+1, k+1, i  ) )
 
      i=idim1-1
 
      dx(2) = dx(1) + cx - 0.125 * ( 
     . x(j  , k  , i  ) + x(j  , k  , i+1) +
     . x(j  , k+1, i  ) + x(j  , k+1, i+1) +
     . x(j+1, k  , i  ) + x(j+1, k  , i+1) +
     . x(j+1, k+1, i  ) + x(j+1, k+1, i+1) )
 
      dy(2) = dy(1) + cy - 0.125 * ( 
     . y(j  , k  , i  ) + y(j  , k  , i+1) +
     . y(j  , k+1, i  ) + y(j  , k+1, i+1) +
     . y(j+1, k  , i  ) + y(j+1, k  , i+1) +
     . y(j+1, k+1, i  ) + y(j+1, k+1, i+1) )
 
      dz(2) = dz(1) + cz - 0.125 * ( 
     . z(j  , k  , i  ) + z(j  , k  , i+1) +
     . z(j  , k+1, i  ) + z(j  , k+1, i+1) +
     . z(j+1, k  , i  ) + z(j+1, k  , i+1) +
     . z(j+1, k+1, i  ) + z(j+1, k+1, i+1) )
 
      cx = cx - xcentrot
      cy = cy - ycentrot
      cz = cz - zcentrot
 
      do 600 l=3,4
 
      gx = cx + dx(l-2)
      gy = cy + dy(l-2)
      gz = cz + dz(l-2)
 
c     Uinf + r x omega is the same as Uinf - omega x r 
 
      qi0(j,k,2,l) = qi0(j,k,2,l) + ( gy * wz - gz * wy )
      qi0(j,k,3,l) = qi0(j,k,3,l) + ( gz * wx - gx * wz )
      qi0(j,k,4,l) = qi0(j,k,4,l) + ( gx * wy - gy * wx )
      bci(j,i,2)   = 0.0
  600 continue
      end if
c
 1001 format(' ','  j=   1  NONINERTIAL freestream         type 1000',
     .       '  i=',i3,',',i3,'  k=',i3,',',i3)
 1002 format(' ','  j=jdim  NONINERTIAL freestream         type 1000',
     .       '  i=',i3,',',i3,'  k=',i3,',',i3)
 1003 format(' ','  k=   1  NONINERTIAL freestream         type 1000',
     .       '  i=',i3,',',i3,'  j=',i3,',',i3)
 1004 format(' ','  k=kdim  NONINERTIAL freestream         type 1000',
     .       '  i=',i3,',',i3,'  j=',i3,',',i3)
 1005 format(' ','  i=   1  NONINERTIAL freestream         type 1000',
     .       '  j=',i3,',',i3,'  k=',i3,',',i3)
 1006 format(' ','  i=idim  NONINERTIAL freestream         type 1000',
     .       '  j=',i3,',',i3,'  k=',i3,',',i3)
c
      return
      end
