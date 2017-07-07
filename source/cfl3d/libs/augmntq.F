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
      subroutine augmntq(q,jdim,kdim,idim,nbl,ldim,qj0,qk0,qi0,qq,
     .                   bcj,bck,bci,nou,bou,nbuf,ibufdim,ighost)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Create an "augmented" q array (qq) of dimensions 
c     jdim+1 x kdim+1 x idim+1, in which the min/max indicies
c     contain cell-face center data and the interior points contain
c     cell-center data at the standard cfl3d cell center locations.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension q(jdim,kdim,idim,ldim),qi0(jdim,kdim,ldim,4),
     .          qj0(kdim,idim-1,ldim,4),qk0(jdim,idim-1,ldim,4)
      dimension qq(jdim+1,kdim+1,idim+1,ldim)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
c
      common /twod/ i2d
c
c     load interior q data into qq array
c
      do l=1,ldim
         do i=1,idim-1
            ii = i+1
            do k=1,kdim-1
               kk = k+ 1
               do j=1,jdim-1
                  jj = j+1
                  qq(jj,kk,ii,l) = q(j,k,i,l)
               end do
            end do
         end do
      end do
c
c     fill in block faces (excluding for now the edges and corners)
c     if ghost-cell data is not available, use duplicate interior
c     data to fill in the block face info.
c
c     note: eddy viscosity and turbulence quantities are always 
c     stored at cell centers, so must override the bcj/k/i arrays
c
      visflg = 1.
      if (ldim .ne. 5) visflg = 0.
c
      if (ighost .ne. 0) then
c
c        calculate and store cell-face center data into qq
c
         do l=1,ldim
            do i=1,idim-1
               ii = i+1
               do k=1,kdim-1
                  kk = k+1
                  aa = 1. + bcj(k,i,1)*visflg
                  bb = 1. - bcj(k,i,1)*visflg
                  cc = 1. + bcj(k,i,2)*visflg
                  dd = 1. - bcj(k,i,2)*visflg
                  qq(1,kk,ii,l)      = (aa*qj0(k,i,l,1) +
     .                                  bb*q(1,k,i,l))*0.5
                  qq(jdim+1,kk,ii,l) = (cc*qj0(k,i,l,3) +
     .                                  dd*q(jdim-1,k,i,l))*0.5     
               end do
            end do
            do i=1,idim-1
               ii = i+1
               do j=1,jdim-1
                  jj = j+1
                  aa = 1. + bck(j,i,1)*visflg
                  bb = 1. - bck(j,i,1)*visflg
                  cc = 1. + bck(j,i,2)*visflg
                  dd = 1. - bck(j,i,2)*visflg
                  qq(jj,1,ii,l)      = (aa*qk0(j,i,l,1) +
     .                                  bb*q(j,1,i,l))*0.5
                  qq(jj,kdim+1,ii,l) = (cc*qk0(j,i,l,3) +
     .                                  dd*q(j,kdim-1,i,l))*0.5
               end do
            end do
            do k=1,kdim-1
               kk = k+1
               do j=1,jdim-1
                  jj = j+1
                  aa = 1. + bci(j,k,1)*visflg
                  bb = 1. - bci(j,k,1)*visflg
                  cc = 1. + bci(j,k,2)*visflg
                  dd = 1. - bci(j,k,2)*visflg
                  qq(jj,kk,1,l)      = (aa*qi0(j,k,l,1) + 
     .                                  bb*q(j,k,1,l))*0.5
                  qq(jj,kk,idim+1,l) = (cc*qi0(j,k,l,3) +
     .                                  dd*q(j,k,idim-1,l))*0.5
               end do
            end do
         end do
c
      else
c
c        load duplicate interior data in block faces (1st order approx)
c
         do l=1,ldim
            do i=1,idim-1
               ii = i+1
               do k=1,kdim-1
                  kk = k+1
                  qq(1,kk,ii,l)      = qq(2,kk,ii,l)
                  qq(jdim+1,kk,ii,l) = qq(jdim,kk,ii,l)
               end do
            end do
            do i=1,idim-1
               ii = i+1
               do j=1,jdim-1
                  jj = j+1
                  qq(jj,1,ii,l)      = qq(jj,2,ii,l)
                  qq(jj,kdim+1,ii,l) = qq(jj,kdim,ii,l)
               end do
            end do
            do k=1,kdim-1
               kk = k+1
               do j=1,jdim-1
                  jj = j+1
                  qq(jj,kk,1,l)      = qq(jj,kk,2,l)
                  qq(jj,kk,idim+1,l) = qq(jj,kk,idim,l)
               end do
            end do
         end do
      end if 
c
c     fill in edge values and corner values via extrapolation/averaging
c
      itwo = 2
      if (i2d.ne.0) itwo = 1
c
      do l=1,ldim
c
c        fill in j-k edge values via extapolation in both j and k
c        directions and averaging
c
         do i=2,idim
            do k=1,kdim+1,kdim
               k1 = 1
               k2 = 2
               if (k.eq.kdim+1) then
                  k1 = -1
                  k2 = -2
               end if
               do j=1,jdim+1,jdim
                  j1 = 1
                  j2 = 2
                  if (j.eq.jdim+1) then
                     j1 = -1
                     j2 = -2
                  end if
                  dqj = qq(j+j1,k,i,l) - qq(j+j2,k,i,l)
                  qj  = qq(j+j1,k,i,l) + dqj
                  dqk = qq(j,k+k1,i,l) - qq(j,k+k2,i,l)
                  qk  = qq(j,k+k1,i,l) + dqk
                  qq(j,k,i,l) = (qj + qk)*0.5
               end do
            end do
         end do
c
c        fill in i-j edge values via extapolation in both i and j
c        directions and averaging
c
         do i=1,idim+1,idim
            i1 = 1
            i2 = itwo
            if (i.eq.idim+1) then
               i1 = -1
               i2 = -itwo
            end if
            do k=2,kdim
               do j=1,jdim+1,jdim
                  j1 = 1
                  j2 = 2
                  if (j.eq.jdim+1)  then
                     j1 = -1
                     j2 = -2
                  end if
                  dqj = qq(j+j1,k,i,l) - qq(j+j2,k,i,l)
                  qj  = qq(j+j1,k,i,l) + dqj
                  dqi = qq(j,k,i+i1,l) - qq(j,k,i+i2,l)
                  qi  = qq(j,k,i+i1,l) + dqi
                  qq(j,k,i,l) = (qj + qi)*0.5
               end do
            end do
         end do
c
c        fill in i-k edge values via extapolation in both i and k
c        directions and averaging
c
         do i=1,idim+1,idim
            i1 = 1
            i2 = itwo
            if (i.eq.idim+1) then
               i1 = -1
               i2 = -itwo
            end if
            do k=1,kdim+1,kdim
               k1 = 1
               k2 = 2
               if (k.eq.kdim+1) then
                  k1 = -1
                  k2 = -2
               end if
               do j=2,jdim
                  dqi = qq(j,k,i+i1,l) - qq(j,k,i+i2,l)
                  qi  = qq(j,k,i+i1,l) + dqi
                  dqk = qq(j,k+k1,i,l) - qq(j,k+k2,i,l)
                  qk  = qq(j,k+k1,i,l) + dqk
                  qq(j,k,i,l) = (qi + qk)*0.5
               end do
            end do
         end do
c
c        fill in corner values via extapolation in all 3 directions
c        and averaging
c
         do k=1,kdim+1,kdim
            k1 = 1
            k2 = 2
            if (k.eq.kdim+1) then
               k1 = -1
               k2 = -2
            end if
            do j=1,jdim+1,jdim
               j1 = 1
               j2 = 2
               if (j.eq.jdim+1) then
                  j1 = -1
                  j2 = -2
               end if
               do i=1,idim+1,idim
                  i1 = 1
                  i2 = itwo
                  if (i.eq.idim+1) then
                     i1 = -1
                     i2 = -itwo
                  end if
                  dqi = qq(j,k,i+i1,l) - qq(j,k,i+i2,l)
                  qi  = qq(j,k,i+i1,l) + dqi
                  dqj = qq(j+j1,k,i,l) - qq(j+j2,k,i,l)
                  qj  = qq(j+j1,k,i,l) + dqj
                  dqk = qq(j,k+k1,i,l) - qq(j,k+k2,i,l)
                  qk  = qq(j,k+k1,i,l) + dqk
                  qq(j,k,i,l) = (qi + qj + qk)/3.
               end do
            end do
         end do
      end do
c
      return
      end
