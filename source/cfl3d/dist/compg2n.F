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
      subroutine compg2n (nblock,ngrid,ncgg,nblg,idimg,jdimg,kdimg,
     .                    nblcg,nodes,iwrk,myid,myhost,mblk2nd,mycomm,
     .                    maxgr,maxbl,ierrflg,ibufdim,nbuf,bou,nou)
c
c     $Id$
c
c***********************************************************************
c     Purpose: To assign blocks to nodes based on size considerations
c     and multigrid structure
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension ncgg(maxgr),nblg(maxgr)
      dimension idimg(maxbl),jdimg(maxbl),kdimg(maxbl)
      dimension iwrk(maxbl,3),nblcg(maxbl),mblk2nd(maxbl)
c
c     simple case: 1 node
c
      if (nodes.eq.1) then
         do ibloc = 1,nblock
            mblk2nd(ibloc) = 1
         end do
         return
      end if
c
c     zero out block assignments
c
      do ibloc = 1,nblock
         mblk2nd(ibloc)  = 0
         iwrk(ibloc,3) = 0
      end do
c
c     find parent grid sizes and blocks associated with coarser grids
c     ...sizes are in terms of number of cells
c
      do igrid = 1,ngrid
         ibloc = nblg(igrid)
         iwrk(igrid,1) = igrid
         iwrk(igrid,2) =
     .   (idimg(ibloc)-1)*(jdimg(ibloc)-1)*(kdimg(ibloc)-1)
      end do
c
c     sort grids from largest to smallest
c
      do igrid = 1, ngrid-1
         do jgrid = igrid+1, ngrid 
            if (iwrk(igrid,2) .lt. iwrk(jgrid,2)) then
               itmp          = iwrk(igrid,1)
               iwrk(igrid,1) = iwrk(jgrid,1)
               iwrk(jgrid,1) = itmp
               itmp          = iwrk(igrid,2)
               iwrk(igrid,2) = iwrk(jgrid,2)
               iwrk(jgrid,2) = itmp
            end if
         end do
      end do
c
c     option 1: allow embeded blocks to reside on any processor in
c               order to load balance (some cases may not work 
c               correctly, and is just about as slow as option 2)
c
c     option 2: restrict embeded blocks to the same processor as
c               the block in which they are embeded (many cases
c               won't load balance very well). 
c
c     option 2 is recommended at this point, 1 is only for testing
c
      ioption = 2
c
      if (ioption.eq.1) then
c
c        assign global-level grids in order (largest to smallest) to
c        the nodes with the fewest points
c
         do jgrid = 1,ngrid
            igrid = iwrk(jgrid,1)
            ibloc = nblg(igrid)
            if (mblk2nd(ibloc).eq.0) then
               imin  = 1
               do inode=2,nodes
                  if (iwrk(inode,3).lt.iwrk(imin,3)) imin = inode
               end do
               mblk2nd(ibloc) = imin
               iwrk(imin,3) = iwrk(imin,3) + iwrk(jgrid,2)
            end if
         end do
c
      else
c
c       assign global-level grids in order (largest to smallest) to
c       the nodes with the fewest points, skipping embedded blocks
c
         do jgrid = 1,ngrid
            igrid = iwrk(jgrid,1)
            ibloc = nblg(igrid)
            if (nblcg(ibloc).eq.ibloc) then
c              block is not embeded
               if (mblk2nd(ibloc).eq.0) then
                  imin  = 1
                  do inode=2,nodes
                     if (iwrk(inode,3).lt.iwrk(imin,3)) imin = inode
                  end do
                  mblk2nd(ibloc) = imin
                  iwrk(imin,3) = iwrk(imin,3) + iwrk(jgrid,2)
               end if
            end if
         end do
c
c        assign each embedded block to the same node as the block 
c        into  which it is embeded
c
         do jgrid = 1,ngrid
            igrid = iwrk(jgrid,1)
            ibloc = nblg(igrid)
            if (nblcg(ibloc).ne.ibloc) then
c               block is embeded
                mblk2nd(ibloc) = mblk2nd(nblcg(ibloc))
            end if
         end do
c
      end if
c
c     assign each coarser-level grid to the same node as its parent
c
      nbl = 0
      do igrid = 1, ngrid
         nbl = nbl + 1
         ibloc = nblg(igrid)
         if (nbl .ne. ibloc) then
            write(11,*) 'block numbers mismatch, in compg2n'
            call termn8(myid,ierrflg,ibufdim,nbuf,bou,nou)
         end if
         nodel = mblk2nd(ibloc)
         ncg   = ncgg(igrid)
         if (ncg.gt.0) then
            do n = 1, ncg
               nbl = nbl + 1
               mblk2nd(nbl) = nodel
            end do
         end if
      end do
c
      return
      end
