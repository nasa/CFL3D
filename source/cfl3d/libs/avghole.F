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
      subroutine avghole(q,blank,jdim,kdim,idim,nbl,ldim,int_updt)
c
c     $Id$
c
c***********************************************************************
c      Purpose: Replace the solution at any point with blank = 0 with
c               the average of neighboring points. This affects hole,
c               points, orphan points AND fringe points; however, 
c               fringe points are subsequently updated properly by
c               a call to subroutine xupdt. This allows any orhpan 
c               points to be updated, even though no interpolation
c               stencils are available for them. 
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension q(jdim,kdim,idim,ldim),blank(jdim,kdim,idim)
c
c
      do l=1,ldim
         do i=1,idim-1
            do k=1,kdim-1
               do j=1,jdim-1
                  if (blank(j,k,i) .eq. 0.) then
                     jp = min(j+1,jdim-1)
                     jm = max(j-1,1)
                     kp = min(k+1,kdim-1)
                     km = max(k-1,1)
                     ip = min(i+1,idim-1)
                     im = max(i-1,1)
                     q(j,k,i,l) = (q(jp, k,ip,l) + q(jm, k,ip,l)
     .                          +  q( j,kp,ip,l) + q( j,km,ip,l)
     .                          +  q(jp, k,im,l) + q(jm, k,im,l)
     .                          +  q( j,kp,im,l) + q( j,km,im,l))*0.125
                  end if
               end do
            end do
         end do
      end do
      return
      end
