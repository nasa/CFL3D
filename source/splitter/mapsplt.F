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
      subroutine mapsplt(nbl,ndir,indx,mbloc,msegt,mtot,
     .                   imap,idbloc,ivisb,itrb,val,xdum,iold,
     .                   nxtseg,intrfc,ipatch,nsubbl,idobl,
     .                   nseg,idno,ijk,idseg,idnext,xdmold,idmold)
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
      dimension xdmold(msegt),idmold(msegt)
c
      common /block/ nbltop
c
      if (ndir.eq.1) then
         ix1 = 1
         ix2 = 2
         imn = 1
         imx = 2
         jmn = 3
         jmx = 4
      else if (ndir.eq.2) then
         ix1 = 3
         ix2 = 4
         imn = 3
         imx = 4
         jmn = 5
         jmx = 6
      else
         ix1 = 5
         ix2 = 6
         imn = 5
         imx = 6
         jmn = 1
         jmx = 2
      endif
c
      ibl = idobl(nbl)
      nsub = nsubbl(nbl)
      do 100 n=1,nsub
         if ((ijk(ix1,ibl).lt.indx).and.(ijk(ix2,ibl).gt.indx)) then
            nsubbl(nbl) = nsubbl(nbl) + 1
            nbltop = nbltop + 1
            nbnew = nbltop
            idnext(nbnew) = idobl(nbl)
            idobl(nbl) = nbnew
            nseg(nbnew) = 0
            idseg(nbnew) = 0
            do 200 i=1,6
               ijk(i,nbnew) = ijk(i,ibl)
 200        continue
            do 201 i=1,msegt
               ivisb(i,nbnew) = ivisb(i,ibl)
 201        continue
            do 202 i=1,7
               itrb(i,nbnew) = itrb(i,ibl)
 202        continue
            ijk(ix1,nbnew) = indx
            ijk(ix2,ibl) = indx
            local = (indx+1) - ijk(ix1,ibl)
            loc2  = (ijk(ix2,nbnew)+1) - indx
            icnt = 1
            ipatch(1) = 2
            iseg = idseg(ibl)
 900        if (iseg.eq.0) goto 999
               if (imap(1,iseg).eq.-1) then
                  call ipadd(imap(7,iseg),imap(8,iseg),mbloc,msegt,mtot,
     .                 imap,idbloc,ivisb,itrb,val,xdum,iold,
     .                 nxtseg,intrfc,ipatch,nsubbl,idobl,
     .                 nseg,idno,ijk,idseg,idnext)
               endif
               nft = imap(2,iseg)
               if (nft.eq.imn) then
                  if (imap(13,iseg).ne.0) then
                     if (imap(19,iseg).gt.local) imap(19,iseg) = local
                     if (imap(20,iseg).gt.local) imap(20,iseg) = local
                  endif
                  iseg = nxtseg(iseg)
               else if (nft.eq.imx) then
                  if (imap(13,iseg).ne.0) then
                     if (imap(19,iseg).gt.loc2) imap(19,iseg) = loc2
                     if (imap(20,iseg).gt.loc2) imap(20,iseg) = loc2
                  endif
                  itmp = nxtseg(iseg)
                  call movseg(iseg,ibl,nbnew,mbloc,msegt,mtot,
     .                        imap,idbloc,ivisb,itrb,val,xdum,iold,
     .                        nxtseg,intrfc,ipatch,nsubbl,idobl,
     .                        nseg,idno,ijk,idseg,idnext)
                  iseg = itmp
               else if ( (nft.eq.jmn).or.(nft.eq.jmx) ) then
                  if (imap(6,iseg).le.local) then
                     iseg = nxtseg(iseg)
                  else if (imap(5,iseg).ge.local) then
                     itmp = nxtseg(iseg)
                     call movseg(iseg,ibl,nbnew,mbloc,msegt,mtot,
     .                           imap,idbloc,ivisb,itrb,val,xdum,iold,
     .                           nxtseg,intrfc,ipatch,nsubbl,idobl,
     .                           nseg,idno,ijk,idseg,idnext)
                     call renmbr(iseg,5,6,local,mbloc,msegt,mtot,
     .                           imap,idbloc,ivisb,itrb,val,xdum,iold,
     .                           nxtseg,intrfc,ipatch,nsubbl,idobl,
     .                           nseg,idno,ijk,idseg,idnext)
                     iseg = itmp
                  else
                     call spltsg(iseg,5,6,local,ibl,nbnew,mbloc,
     .                           msegt,mtot,imap,idbloc,
     .                           ivisb,itrb,val,xdum,iold,nxtseg,
     .                           intrfc,ipatch,nsubbl,idobl,
     .                           nseg,idno,ijk,idseg,idnext)

                     iseg = idseg(ibl)
                  endif
               else
                  if (imap(4,iseg).le.local) then
                     iseg = nxtseg(iseg)
                  else if (imap(3,iseg).ge.local) then
                     itmp = nxtseg(iseg)
                     call movseg(iseg,ibl,nbnew,mbloc,msegt,mtot,
     .                           imap,idbloc,ivisb,itrb,val,xdum,iold,
     .                           nxtseg,intrfc,ipatch,nsubbl,idobl,
     .                           nseg,idno,ijk,idseg,idnext)
                     call renmbr(iseg,3,4,local,mbloc,msegt,mtot,
     .                           imap,idbloc,ivisb,itrb,val,xdum,iold,
     .                           nxtseg,intrfc,ipatch,nsubbl,idobl,
     .                           nseg,idno,ijk,idseg,idnext)
                     iseg = itmp
                  else
                     call spltsg(iseg,3,4,local,ibl,nbnew,mbloc,
     .                           msegt,mtot,imap,idbloc,
     .                           ivisb,itrb,val,xdum,iold,nxtseg,
     .                           intrfc,ipatch,nsubbl,idobl,
     .                           nseg,idno,ijk,idseg,idnext)
                     iseg = idseg(ibl)
                  endif
               endif
               icnt = icnt + 1
               goto 900
 999        continue
            call mkintr(imn,ibl,nbnew,mbloc,msegt,mtot,
     .                  imap,idbloc,ivisb,itrb,val,xdum,iold,
     .                  nxtseg,intrfc,ipatch,nsubbl,idobl,
     .                  nseg,idno,ijk,idseg,idnext,xdmold,idmold)
            call ipfix(ibl,ndir,local,nbnew,mbloc,msegt,mtot,
     .                 imap,idbloc,ivisb,itrb,val,xdum,iold,
     .                 nxtseg,intrfc,ipatch,nsubbl,idobl,
     .                 nseg,idno,ijk,idseg,idnext)
         endif
         ibl = idnext(ibl)
 100  continue
      return
      end
