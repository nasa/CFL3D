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
	subroutine get_bvals(jdim1,kdim1,idim1,jmax1,kmax1,
     .                       ldim,q,qbou,mtype,i1)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  To pick out values at the boundary from global values
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
        dimension q(jdim1, kdim1, idim1, ldim)
        dimension qbou(jmax1, kmax1, ldim)
c
        if (mtype. eq.1) then
c
           do j = 1,jmax1
              do k = 1,kmax1
                 do l = 1,ldim
                    qbou(j,k,l) = q(j,k,i1,l)
                 end do
              end do
           end do
c
        else if (mtype. eq. 2) then
c
           do j = 1,jmax1
              do k = 1,kmax1
                 do l = 1,ldim
                    qbou(j,k,l) = q(i1,j,k,l)
                 end do
              end do
           end do             
c
        else if (mtype. eq. 3) then
c
           do j = 1,jmax1
              do k = 1,kmax1
                 do l = 1,ldim
                    qbou(j,k,l) = q(j,i1,k,l)
                 end do
              end do
           end do
c
        end if
c
        return
        end
