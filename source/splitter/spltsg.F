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
      subroutine spltsg(iseg,ix1,ix2,indx,is,it,mbloc,msegt,mtot,
     .                  imap,idbloc,ivisb,itrb,val,xdum,iold,
     .                  nxtseg,intrfc,ipatch,nsubbl,idobl,
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
c
      common /segment/ nsgtop
c
      call mkseg(iseg,it,ns,mbloc,msegt,mtot,
     .           imap,idbloc,ivisb,itrb,val,xdum,iold,
     .           nxtseg,intrfc,ipatch,nsubbl,idobl,
     .           nseg,idno,ijk,idseg,idnext)
c
      imap(ix2,iseg) = indx
      imap(ix1,ns) = indx
      if (imap(13,iseg).ne.0) then
         if (imap(ix1+11,iseg).ge.indx) then
            imap(13,iseg) = 0
         else if (imap(ix2+11,iseg).gt.indx) then
            imap(ix2+11,iseg) = indx
         endif
         if (imap(ix2+11,ns).le.indx) then
            imap(13,ns) = 0
         else if (imap(ix1+11,ns).lt.indx) then
            imap(ix1+11,ns) = indx
         endif
      endif
      iold(ix2-2,iseg) = iold(ix1-2,iseg) + (indx-imap(ix1,iseg))
      iold(ix1-2,ns) = iold(ix2-2,iseg)
      jseg = intrfc(iseg)
      if (jseg.gt.0) then
         if (imap(8,iseg).gt.0) then
            jx1 = ix1
            jx2 = ix2
            is1 = ix1 + 6
            is2 = ix2 + 6
         else if (ix1.eq.3) then
            jx1 = 5
            jx2 = 6
            is1 = 11
            is2 = 12
         else
            jx1 = 3
            jx2 = 4
            is1 = 9
            is2 = 10
         endif
         ia = (imap(is2,iseg)-imap(is1,iseg)) / 
     .        (imap(ix2,ns)-imap(ix1,iseg))
         ib = imap(is1,iseg) - ia*imap(ix1,iseg)
         jndx = ia*indx + ib
         if (iseg.eq.jseg) then
            if (jndx.gt.imap(is1,iseg)) then
               imap(is2,iseg) = jndx
               imap(is1,ns) = jndx
               intrfc(ns) = ns
            else if (jndx.eq.indx) then
               imap(is2,iseg) = jndx
               imap(is1,ns) = jndx
               intrfc(ns) = iseg
               intrfc(iseg) = ns
            else if (jndx.gt.indx) then
               call splt2(ns,ix1,ix2,jndx,it,it,mbloc,msegt,mtot,
     .                    imap,idbloc,ivisb,itrb,val,xdum,iold,
     .                    nxtseg,intrfc,ipatch,nsubbl,idobl,
     .                    nseg,idno,ijk,idseg,idnext)
               imap(ix1,ns) = 1
               imap(ix2,ns) = jndx - indx + 1
               imap(is1,ns) = imap(ix2,ns)
               imap(is2,ns) = imap(ix1,ns)
               intrfc(ns) = ns
               imap(7,ns) = it
               imap(8,ns) = imap(2,ns)
               ns = nsgtop
               intrfc(ns) = iseg
               intrfc(iseg) = ns
            else if (jndx.lt.indx) then
               call splt2(iseg,ix1,ix2,jndx,is,is,mbloc,msegt,mtot,
     .                    imap,idbloc,ivisb,itrb,val,xdum,iold,
     .                    nxtseg,intrfc,ipatch,nsubbl,idobl,
     .                    nseg,idno,ijk,idseg,idnext)
               nsold = ns
               ns = nsgtop
               imap(ix1,ns) = jndx
               imap(ix2,ns) = indx
               imap(is1,ns) = imap(ix2,ns)
               imap(is2,ns) = imap(ix1,ns)
               intrfc(ns) = ns
               imap(7,ns) = is
               imap(8,ns) = imap(2,ns)
               ns = nsold
               intrfc(ns) = iseg
               intrfc(iseg) = ns
            endif
            call renmbr(iseg,ix1,ix2,1,mbloc,msegt,mtot,
     .                  imap,idbloc,ivisb,itrb,val,xdum,iold,
     .                  nxtseg,intrfc,ipatch,nsubbl,idobl,
     .                  nseg,idno,ijk,idseg,idnext)
         else
            imap(is2,iseg) = jndx
            imap(is1,ns) = jndx
            ibl = idbloc(jseg)
            call splt2(jseg,jx1,jx2,jndx,ibl,ibl,mbloc,msegt,mtot,
     .                 imap,idbloc,ivisb,itrb,val,xdum,iold,
     .                 nxtseg,intrfc,ipatch,nsubbl,idobl,
     .                 nseg,idno,ijk,idseg,idnext)
            if (jndx.gt.imap(is1,iseg)) then
               intrfc(ns) = idseg(ibl)
               intrfc(idseg(ibl)) = ns
               intrfc(jseg) = iseg
            else
               intrfc(ns) = jseg
               intrfc(idseg(ibl)) = iseg
               intrfc(jseg) = ns
               intrfc(iseg) = idseg(ibl)
            endif
            call renmbr(iseg,ix1,ix2,1,mbloc,msegt,mtot,
     .                  imap,idbloc,ivisb,itrb,val,xdum,iold,
     .                  nxtseg,intrfc,ipatch,nsubbl,idobl,
     .                  nseg,idno,ijk,idseg,idnext)
         endif
      endif
      call renmbr(ns,ix1,ix2,indx,mbloc,msegt,mtot,
     .            imap,idbloc,ivisb,itrb,val,xdum,iold,
     .            nxtseg,intrfc,ipatch,nsubbl,idobl,
     .            nseg,idno,ijk,idseg,idnext)
c
      return
      end
