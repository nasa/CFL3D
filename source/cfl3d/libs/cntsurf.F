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
      subroutine cntsurf(nsurf,maxbl,maxgr,maxseg,ngrid,nblg,nbci0,
     .                   nbcj0,nbck0,nbcidim,nbcjdim,nbckdim,ibcinfo,
     .                   jbcinfo,kbcinfo,ibctyp)
c
c     $Id$
c
c***********************************************************************
c     Purpose: To count the number of surface points (nsurf) with
c     boundary condition type ibctyp 
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension nbci0(maxbl),nbcidim(maxbl),nbcj0(maxbl),nbcjdim(maxbl),
     .          nbck0(maxbl),nbckdim(maxbl),ibcinfo(maxbl,maxseg,7,2),
     .          jbcinfo(maxbl,maxseg,7,2),kbcinfo(maxbl,maxseg,7,2)
      dimension nblg(maxgr)
c
      nsurf=0
      do igrid = 1,ngrid
         nbl  = nblg(igrid)
c
c        i=constant surfaces
c
         do m = 1,2
         if (m.eq.1) then
            ns    = nbci0(nbl)
            nface = 1
         else
            ns    = nbcidim(nbl)
            nface = 2
         end if
            do iseg = 1,ns
               nbctype  =  ibcinfo(nbl,iseg,1,m)
               n1beg    =  ibcinfo(nbl,iseg,2,m)
               n1end    =  ibcinfo(nbl,iseg,3,m)
               n2beg    =  ibcinfo(nbl,iseg,4,m)
               n2end    =  ibcinfo(nbl,iseg,5,m)
               if(abs(nbctype).eq.ibctyp) then
                  nsurf=nsurf+(n2end-n2beg+1)*(n1end-n1beg+1)
               end if
            end do
         end do

c        j=constant surfaces
c
         do m = 1,2
         if (m.eq.1) then
            ns    = nbcj0(nbl)
            nface = 3
         else
            ns    = nbcjdim(nbl)
            nface = 4
         end if
            do iseg = 1,ns
               nbctype  =  jbcinfo(nbl,iseg,1,m)
               n1beg    =  jbcinfo(nbl,iseg,4,m)
               n1end    =  jbcinfo(nbl,iseg,5,m)
               n2beg    =  jbcinfo(nbl,iseg,2,m)
               n2end    =  jbcinfo(nbl,iseg,3,m)
               if(abs(nbctype).eq.ibctyp) then
                  nsurf=nsurf+(n2end-n2beg+1)*(n1end-n1beg+1)
               end if
            end do
         end do
c
c        k=constant surfaces
c
         do m = 1,2
         if (m.eq.1) then
            ns    = nbck0(nbl)
            nface = 5
         else
            ns    = nbckdim(nbl)
            nface = 6
         end if
            do iseg = 1,ns
               nbctype  =  kbcinfo(nbl,iseg,1,m)
               n1beg    =  kbcinfo(nbl,iseg,2,m)
               n1end    =  kbcinfo(nbl,iseg,3,m)
               n2beg    =  kbcinfo(nbl,iseg,4,m)
               n2end    =  kbcinfo(nbl,iseg,5,m)
               if(abs(nbctype).eq.ibctyp) then
                  nsurf=nsurf+(n2end-n2beg+1)*(n1end-n1beg+1)
               end if
            end do
         end do

      end do
c
      return
      end
