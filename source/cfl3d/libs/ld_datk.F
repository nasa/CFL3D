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
      subroutine ld_datk(data,jdim,kdim,idim,datak,ldim,mp,np,
     .                   ista,iend,jsta,jend,ksta,kend,nou,bou,nbuf,
     .                   ibufdim,myid)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Load data from full 3d array to smaller work array for
c               efficient message passing
c***********************************************************************
c
c     Description of variables:
c       data              = full 3D array of data array
c       jdim,kdim,idim    = dimensions of 3D data array
c       datak             = np-plane subset of data in k direction
c       mp                = array containing plane indicies to store
c       np                = number of planes of data needed (up to 4)
c       ista,iend,etc.    = starting and ending values for the
c                           boundary segment of data array to be loaded
c                           into the subset arrays
c       ldim              = number of variables in data array; e.g.
c                           ldim = 5 for q data, ldim = nummem for turb.
c                           variables, ldim = 1 for eddy viscosity,
c                           ldim = 3 for grid data
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension data(jdim,kdim,idim,ldim),
     .          datak(jdim,np,idim,ldim)
      dimension mp(4)
c
      if (np.gt.4) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),*)'stopping...increase dimension of mp ',
     .                'array in ld_datk and its calling routine'
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
c
c     load data values into dataX array:
c
      do i=ista,iend
         do j=jsta,jend
           do m=1,np
              do l=1,ldim
                 datak(j,m,i,l)=data(j,mp(m),i,l)
              end do
            end do
         end do
      end do
c
      return
      end
