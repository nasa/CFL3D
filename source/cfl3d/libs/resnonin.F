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
      subroutine resnonin(nbl,jdim,kdim,idim,q,x,y,z,sj,sk,si,vol,res,
     .                    nou,bou,nbuf,ibufdim)
c
c     $Id$
c
c***********************************************************************
c     Purpose: Compute the residual contributions to the right-hand-side
c                (due only to performing calculations using relative
c                 velocities in a rotating noninertial reference frame.)
c     Original coding by Mike Park
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension q(jdim*kdim*idim,5)
      dimension x(jdim*kdim*idim),y(jdim*kdim*idim),z(jdim*kdim*idim)
      dimension si(jdim*kdim*idim,5),sj(jdim*kdim*(idim-1),5),
     .          sk(jdim*kdim*(idim-1),5)
      dimension vol(jdim*kdim*(idim-1))
      dimension res(jdim*kdim*(idim-1),5)
 
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .        tur10(7)
      common /noninertial/ xcentrot,ycentrot,zcentrot,xrotrate,
     .                     yrotrate,zrotrate,noninflag
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /sklton/ isklton
 
      if (isklton.gt.0) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),101)
 101        format(34h   adding NONINERTIAL source terms)
      end if
 
c  The residual contributions to the right-hand-side (due only to 
c  performing calculations using velocities in a rotating 
c  noninertial reference frame) are broken into three parts: the 
c  contribution due to centripetal forces f(omega, grid position), 
c  the contribution due to Coriolis effects f(omega, q), and the 
c  acceleration of the reference frame origin f(omega, freestream 
c  velocity). All parts are calculated in this section, but the 
c  centripetal and reference frame acceleration calculations are 
c  independent of the solution or q so it could be moved to the 
c  code initialization and stored in memory to increase speed by 
c  avoiding computing these terms at each iteration.
 
c  rename xrotrate, yrotrate, zrotrate -> wx, wy, wz to save typing
                  
      wx = xrotrate
      wy = yrotrate
      wz = zrotrate
 
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
 
c*********************
c     centripetal part
c*********************
 
c compute Omega x Uinf
c  (constant with respect to solution or q) 
 
      OmegaxUx =  ( wy * qiv(4) - wz * qiv(3) )      
      OmegaxUy =  ( wz * qiv(2) - wx * qiv(4) )
      OmegaxUz =  ( wx * qiv(3) - wy * qiv(2) )
 
      do i = 1,idim1
        do k = 1,kdim1
          do j = 1,jdim1
 
c Find the indices corresponding to the 8 corners of a volume
 
      ind000 = (i-1) * jdim * kdim + (k-1) * jdim + (j-1) + 1
      ind001 = (i-1) * jdim * kdim + (k-1) * jdim + (j  ) + 1
      ind010 = (i-1) * jdim * kdim + (k  ) * jdim + (j-1) + 1
      ind011 = (i-1) * jdim * kdim + (k  ) * jdim + (j  ) + 1
      ind100 = (i  ) * jdim * kdim + (k-1) * jdim + (j-1) + 1
      ind101 = (i  ) * jdim * kdim + (k-1) * jdim + (j  ) + 1
      ind110 = (i  ) * jdim * kdim + (k  ) * jdim + (j-1) + 1
      ind111 = (i  ) * jdim * kdim + (k  ) * jdim + (j  ) + 1
 
c  find the radius vector (rx, ry, rz) from the center of rotation 
c  to the center of the volume 
 
      rx = 0.125 * ( 
     .   x(ind000) + 
     .   x(ind001) + 
     .   x(ind010) + 
     .   x(ind011) + 
     .   x(ind100) + 
     .   x(ind101) + 
     .   x(ind110) + 
     .   x(ind111) ) - xcentrot
 
      ry = 0.125 * ( 
     .   y(ind000) + 
     .   y(ind001) + 
     .   y(ind010) + 
     .   y(ind011) + 
     .   y(ind100) + 
     .   y(ind101) + 
     .   y(ind110) + 
     .   y(ind111) ) - ycentrot
 
      rz = 0.125 * ( 
     .   z(ind000) + 
     .   z(ind001) + 
     .   z(ind010) + 
     .   z(ind011) + 
     .   z(ind100) + 
     .   z(ind101) + 
     .   z(ind110) + 
     .   z(ind111) ) - zcentrot 
 
c  compute centripetal pseudo-acceleration vector 
c  (constant with respect to q) = -w x (w x r)
 
      centaccx = wy * (wx*ry - wy*rx) - wz * (wz*rx - wx*rz) 
      centaccy = wz * (wy*rz - wz*ry) - wx * (wx*ry - wy*rx) 
      centaccz = wx * (wz*rx - wx*rz) - wy * (wy*rz - wz*ry) 
 
c  Coriolis part of pseudo-acceleration
 
      corx = 2. * ( wy * q(ind000,4) - wz * q(ind000,3) )
      cory = 2. * ( wz * q(ind000,2) - wx * q(ind000,4) )
      corz = 2. * ( wx * q(ind000,3) - wy * q(ind000,2) )
 
      totx = -OmegaxUx + centaccx + corx
      toty = -OmegaxUy + centaccy + cory
      totz = -OmegaxUz + centaccz + corz
 
c  increment residual (in conserved vars.)
 
      res(ind000,2) = res(ind000,2) + q(ind000,1) * totx * vol(ind000)
      res(ind000,3) = res(ind000,3) + q(ind000,1) * toty * vol(ind000)
      res(ind000,4) = res(ind000,4) + q(ind000,1) * totz * vol(ind000)
      res(ind000,5) = res(ind000,5) + q(ind000,1) * vol(ind000) * (
     . totx * q(ind000,2) + 
     . toty * q(ind000,3) + 
     . totz * q(ind000,4) )
 
          enddo
        enddo
      enddo

      return 
      end
