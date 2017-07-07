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
      subroutine ipfix(ibl,ndir,indx,inewbl,mbloc,msegt,mtot,
     .                 imap,idbloc,ivisb,itrb,val,xdum,iold,
     .                 nxtseg,intrfc,ipatch,nsubbl,idobl,
     .                 nseg,idno,ijk,idseg,idnext)
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
c	write(6,*) 'Fix: ',ibl,ndir,indx,inewbl
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

      itop = ipatch(1)
      i = 2
 1    if (i.ge.itop) goto 999
      j = i + 1
      jbl = ipatch(i)
      jfc = ipatch(j)
c      write(6,*) 'Fixing block ',jbl,', face ',jfc
      is = idseg(jbl)
 2    if (is.eq.0) goto 998
      
      if (imap(1,is).ne.-1)  goto 997
      if (imap(2,is).ne.jfc) goto 997
      if (imap(7,is).ne.ibl) goto 997
      nft = abs(imap(8,is))
      if (imap(8,is).gt.0) then
	 if (nft.eq.jmn .or. nft.eq.jmx) then
	    jx1 = 11
	    jx2 = 12
         else
	    jx1 = 9
	    jx2 = 10
	 endif
      else
	 if (nft.eq.jmn .or. nft.eq.jmx) then
	    jx1 = 9
	    jx2 = 10
         else
	    jx1 = 11
	    jx2 = 12
	 endif
      endif
      if (imap(jx1,is).gt.imap(jx2,is)) then
         jtmp = jx1
         jx1 = jx2
         jx2 = jtmp
      endif
c      write(6,*) 'jx1,jx2: ',jx1,jx2
      if (nft.eq.imn) goto 997
      if (nft.eq.imx) then
         imap(7,is) = inewbl
      else if (imap(jx1,is).ge.indx) then
         imap(7,is) = inewbl
         imap(jx1,is) = (imap(jx1,is)+1) - indx
         imap(jx2,is) = (imap(jx2,is)+1) - indx
      else if (imap(jx2,is).gt.indx) then
         call mkseg(is,jbl,ns,mbloc,msegt,mtot,
     .              imap,idbloc,ivisb,itrb,val,xdum,iold,
     .              nxtseg,intrfc,ipatch,nsubbl,idobl,
     .              nseg,idno,ijk,idseg,idnext)
         imap(jx2,is) = indx
         imap(7,ns) = inewbl
         imap(jx1,ns) = 1
         imap(jx2,ns) = (imap(jx2,ns)+1) - indx
c	 write(6,*) 'new seg: ',imap(jx1,ns),imap(jx2,ns)
      endif
 997  continue
      is = nxtseg(is)
      goto 2
      
 998  continue
      i = j + 1
      goto 1
      
 999  continue

      return
      end
