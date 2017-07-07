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
      subroutine mkintr(imn,i1,i2,mbloc,msegt,mtot,
     .                  imap,idbloc,ivisb,itrb,val,xdum,iold,
     .                  nxtseg,intrfc,ipatch,nsubbl,idobl,
     .                  nseg,idno,ijk,idseg,idnext,xdmold,idmold)
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
      common /segment/ nsgtop
c
      imx = 0
      jmx = 0
      iomn = 1000000
      iomx = 0
      jomn = 1000000
      jomx = 0
      iseg = idseg(i1)
 100  if (iseg.eq.0) goto 199
         if (imap(2,iseg).eq.imn) then
            do 101 i=1,msegt
               xdmold(i) = xdum(i,iseg)
 101        continue
            do i=21,msegt
               idmold(i) = imap(i,iseg)
            end do
            if (imap(4,iseg).gt.imx) imx = imap(4,iseg)
            if (imap(6,iseg).gt.jmx) jmx = imap(6,iseg)
            if (iold(1,iseg).lt.iomn) iomn = iold(1,iseg)
            if (iold(2,iseg).gt.iomx) iomx = iold(2,iseg)
            if (iold(3,iseg).lt.jomn) jomn = iold(3,iseg)
            if (iold(4,iseg).gt.jomx) jomx = iold(4,iseg)
         endif
         iseg = nxtseg(iseg)
         goto 100
 199  continue
      nsgtop = nsgtop + 1
      iseg = nsgtop
      nxtseg(iseg) = idseg(i1)
      idseg(i1) = iseg
      nseg(i1) = nseg(i1) + 1
      idbloc(iseg) = i1
      nsgtop = nsgtop + 1
      jseg = nsgtop
      nxtseg(jseg) = idseg(i2)
      idseg(i2) = jseg
      nseg(i2) = nseg(i2) + 1
      idbloc(jseg) = i2
      intrfc(iseg) = jseg
      intrfc(jseg) = iseg
      iold(1,iseg) = iomn
      iold(2,iseg) = iomx
      iold(3,iseg) = jomn
      iold(4,iseg) = jomx
      iold(1,jseg) = iomn
      iold(2,jseg) = iomx
      iold(3,jseg) = jomn
      iold(4,jseg) = jomx
      imap(1,iseg) = 1
      imap(2,iseg) = imn+1
      imap(3,iseg) = 1
      imap(4,iseg) = imx
      imap(5,iseg) = 1
      imap(6,iseg) = jmx
      imap(7,iseg) = i2
      imap(8,iseg) = imn
      imap(9,iseg) = 1
      imap(10,iseg) = imx
      imap(11,iseg) = 1
      imap(12,iseg) = jmx
      imap(1,jseg) = 1
      imap(2,jseg) = imn
      imap(3,jseg) = 1
      imap(4,jseg) = imx
      imap(5,jseg) = 1
      imap(6,jseg) = jmx
      imap(7,jseg) = i1
      imap(8,jseg) = imn+1
      imap(9,jseg) = 1
      imap(10,jseg) = imx
      imap(11,jseg) = 1
      imap(12,jseg) = jmx
      do 200 i=13,20
         imap(i,iseg) = 0
         imap(i,jseg) = 0
 200  continue
      do 201 i=1,msegt
         xdum(i,iseg) = xdmold(i)
         xdum(i,jseg) = xdmold(i)
 201  continue
         imap(21,iseg) = 0  ! no bc2000 series data on new interface
         imap(21,jseg) = 0
         imap(22,iseg) = 1  ! don't turn off force on new interface
         imap(22,jseg) = 1
         imap(23,iseg) = idmold(23) ! preserve iovrlp on new interface
         imap(23,jseg) = idmold(23)
      val(iseg) = 0.0
      val(jseg) = 0.0
c
      return
      end
