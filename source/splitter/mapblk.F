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
      subroutine mapblk(nbl,id,imin,imax,jmin,jmax,kmin,kmax,mbloc,
     .                  msegt,mtot,imap,idbloc,ivisb,itrb,val,xdum,
     .                  iold,nxtseg,intrfc,ipatch,nsubbl,idobl,
     .                  nseg,idno,ijk,idseg,idnext)
c
c     $Id$
c
      dimension nsubbl(mbloc),idobl(mbloc),nseg(mbloc),idno(mbloc),
     .          ijk(6,mbloc),idseg(mbloc),idnext(mbloc)
      dimension imap(msegt,mtot),idbloc(mtot),
     .          ivisb(msegt,mbloc),itrb(7,mbloc),
     .          val(mtot),xdum(msegt,mtot),
     .          iold(4,mtot),nxtseg(mtot),intrfc(mtot),
     .          ipatch(mtot)
      dimension imin(*),imax(*),jmin(*),jmax(*),kmin(*),kmax(*)
      dimension id(*)
c
      common /oldbl/ noldbl
c
      nbl = 0
      do 100 n=1,noldbl
         ibl = idobl(n)
         do 101 n2=1,nsubbl(n)
            nbl = nbl + 1
            id(nbl) = n
            imin(nbl) = ijk(1,ibl)
            imax(nbl) = ijk(2,ibl)
            jmin(nbl) = ijk(3,ibl)
            jmax(nbl) = ijk(4,ibl)
            kmin(nbl) = ijk(5,ibl)
            kmax(nbl) = ijk(6,ibl)
            ibl = idnext(ibl)
 101     continue
 100  continue
c
      return
      end
