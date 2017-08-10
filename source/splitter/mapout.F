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
      subroutine mapout(iwrit,nbloc,ni,nj,nk,mbloc,msegt,mtot,
     .                  imap,idbloc,ivisb,itrb,val,xdum,iold,
     .                  nxtseg,intrfc,ipatch,nsubbl,idobl,nseg,idno,
     .                  ijk,idseg,idnext)
c
c     $Id$
c
c-----------------------------------------------------------------------
c
c     object : to read in and set up the topological mapping 
c              funcions, 'imap'
c
c     initial coding : by m.d.sanetrik (August, 1990)
c
      dimension nsubbl(mbloc),idobl(mbloc),nseg(mbloc),idno(mbloc),
     .          ijk(6,mbloc),idseg(mbloc),idnext(mbloc)
      dimension imap(msegt,mtot),idbloc(mtot),
     .          ivisb(msegt,mbloc),itrb(7,mbloc),
     .          val(mtot),xdum(msegt,mtot),
     .          iold(4,mtot),nxtseg(mtot),intrfc(mtot),
     .          ipatch(mtot),ni(mbloc),nj(mbloc),nk(mbloc)
c
      common /oldbl/ noldbl
      common /titles/ title1,title2,title3
c
      character*80 title1,title2,title3
c
      nbl = 0
      do 100 n=1,noldbl
         ibl = idobl(n)
         do 101 n2=1,nsubbl(n)
            nbl = nbl + 1
            idno(ibl) = nbl
            ibl = idnext(ibl)
 101     continue
 100  continue
c
      do 110 n=1,noldbl
         ibl = idobl(n)
         do 111 n2=1,nsubbl(n)
            ics = idseg(ibl)
            do 112 iseg=1,nseg(ibl)
               ityp = imap(1,ics)
	       if (0.eq.ityp .or. 1.eq.ityp) then
c		  jbl = imap(7,ics)
		  imap(7,intrfc(ics)) = idno(ibl)
	       endif
	       if ((-1).eq.ityp) then
		  jbl = imap(7,ics)
	          imax = (ijk(2,jbl)+1) - ijk(1,jbl)
	          jmax = (ijk(4,jbl)+1) - ijk(3,jbl)
	          kmax = (ijk(6,jbl)+1) - ijk(5,jbl)
		  imap(7,ics) = idno(jbl)
		  jface = imap(8,ics)
		  if (1.eq.jface .or. 2.eq.jface) then
		     mx1 = jmax
		     mx2 = kmax
		  else if ((-1).eq.jface .or. (-2).eq.jface) then
		     mx1 = kmax
		     mx2 = jmax
		  else if (3.eq.jface .or. 4.eq.jface) then
		     mx1 = kmax
		     mx2 = imax
		  else if ((-3).eq.jface .or. (-4).eq.jface) then
		     mx1 = imax
		     mx2 = kmax
		  else if (5.eq.jface .or. 6.eq.jface) then
		     mx1 = imax
		     mx2 = jmax
		  else if ((-5).eq.jface .or. (-6).eq.jface) then
		     mx1 = jmax
		     mx2 = imax
		  endif
		  if(imap(9,ics).gt.mx1)  imap(9,ics) = mx1
		  if(imap(10,ics).gt.mx1) imap(10,ics) = mx1
		  if(imap(11,ics).gt.mx2) imap(11,ics) = mx2
		  if(imap(12,ics).gt.mx2) imap(12,ics) = mx2
	       endif
               ics = nxtseg(ics)
 112        continue
            ibl = idnext(ibl)
 111     continue
 110  continue
c
      write(iwrit,1000)
      write(iwrit,'(1x,i5)') nbl
 1000 format(1x,'nbloc')
 1001 format(1x,i5)
 1002 format(2i5,i6,11i5)
 1003 format(15x,8i5,f5.1)
 1004 format(1x,'  nbt  nst   typ  nft 1tmn 1tmx 2tmn 2tmx  nbs',
     .          '  nfs 1smn 1smx 2smn 2smx')
 1005 format(1x,'                ntrb 1beg 1end 2beg 2end fmx1',
     .          ' fmx2 nvis')
 2000 format(/,/,' new block  old block  i0  i1  j0  j1  k0  k1')
 2001 format(2x,i4,7x,i4,4x,6i4)
      write(6,2000)
      do 200 nn=1,noldbl
         ibloc = idobl(nn)
         do 201 n2 = 1,nsubbl(nn)
           write(6,2001) idno(ibloc),nn,(ijk(l,ibloc),l=1,6)
         write(iwrit,'(1x)')
         write(iwrit,'(a79)') title1
         write(iwrit,1002) nseg(ibloc),ivisb(1,ibloc),
     .        ivisb(2,ibloc),ivisb(3,ibloc),(itrb(n,ibloc),n=1,7)
         write(iwrit,'("begin non-standard map data")')
         write(iwrit,'(12i6)') (ivisb(n,ibloc),n=4,15)
         write(iwrit,'(12i6)') (ivisb(n,ibloc),n=16,27)
         write(iwrit,'(12i6)') (ivisb(n,ibloc),n=28,34)
         write(iwrit,'("end non-standard map data")')
         write(iwrit,'(a79)') title2
         write(iwrit,'(a79)') title3
         loc = idseg(ibloc)
         do 202 iseg  = 1,nseg(ibloc)
c
c            if (imap(1,loc).eq.-1) then
c               imap(7,loc) = idno(imap(7,loc))
c            endif
            write(iwrit,'(1x)')
            write(iwrit,1002) idno(ibloc),iseg,(imap(n,loc),n= 1,12)
            write(iwrit,1003) (imap(n,loc),n=13,20),
     .                        val(loc)
            write(iwrit,'(15x,5f5.1)') 0.0,0.0,0.0,1.0,0.0
            write(iwrit,'("begin non-standard map data")')
            write(iwrit,'(15x,5f10.4)') (xdum(n,loc),n=1,5)
            write(iwrit,'(15x,5f10.4)') (xdum(n,loc),n=6,10)
            write(iwrit,'(15x,5f10.4)') (xdum(n,loc),n=11,15)
            write(iwrit,'(15x,5f10.4)') (xdum(n,loc),n=16,20)
            write(iwrit,'(15x,5f10.4)') (xdum(n,loc),n=21,25)
            write(iwrit,'(15x,5f10.4)') (xdum(n,loc),n=26,26)
            write(iwrit,'(15x,8i5)') (imap(n,loc),n=21,23)
            write(iwrit,'("end non-standard map data")')
c
            loc = nxtseg(loc)
 202     continue
         ibloc = idnext(ibloc)
 201     continue
c
 200  continue
      write(iwrit,'(1x)')

      return
      end
