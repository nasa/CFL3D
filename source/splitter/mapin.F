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
      subroutine mapin(imapfil,iwrit,nbloc,ni,nj,nk,mbloc,msegt,mtot,
     .                 imap,idbloc,ivisb,itrb,val,xdum,iold,
     .                 nxtseg,intrfc,ipatch,nsubbl,idobl,nseg,idno,
     .                 ijk,idseg,idnext,nou,bou,nbuf,ibufdim)
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
      common /block/ nbltop
      common /segment/ nsgtop
      common /titles/ title1,title2,title3
c
      character*80 title1,title2,title3
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
c
c*************************************************************************
c--------------------  mapping function description ---------------------
c
c      imap  : mapping function containing topological information
c      msegt : maximum no. of types of operations/boundary conditons
c              (currently set to 20)
c      msegn : maximum no. of segments permitted on all faces of a block
c              (currently set to 20)
c
c      imap(1 , ) : specifies boundary/operation type
c      imap(2 , ) : specifies face number (1-6)
c      imap(3 , ) : specifies beginning of direction 1
c      imap(4 , ) : specifies end       of direction 1
c      imap(5 , ) : specifies beginning of direction 2
c      imap(6 , ) : specifies end       of direction 2
c
c      imap(8,  ) : if the boundary type is symmetry
c                   takes the value 1,2 or 3 to indicate symmetry 
c                   about x, y or z = constant planes, respectively
c 
c  if the boundary type is a cut/interface
c  additional information described below is required 
c
c      imap(7 , ) : specifies block no. of source segment
c     |imap(8 , )|: specifies face no. of source segment
c      imap(8 , ) > 0 corresponds to direction 1 of source segment
c                      matching with direction 1 of target segment
c                      and same  for direction 2
c      imap(8 , ) < 0 corresponds to direction 2 of source segment
c                      matching with direction 1 of target segment
c                      and vice-versa
c
c      imap(9 , ) : specifies beginning of direction 1 of source segment
c      imap(10, ) : specifies end       of direction 1 of source segment
c      imap(11, ) : specifies beginning of direction 2 of source segment
c      imap(12, ) : specifies end       of direction 2 of source segment
c
c      turbulence/transitional related information
c    
c      imap(13, ) : specifies if there is turbulent flow on this segment
c                   1  corresponds to turbulent flow
c                   0  corresponds to laminar   flow
c      imap(14, ) : begining index in direction 1 for turbulent flow
c      imap(15, ) : ending   index in direction 1 for turbulent flow
c      imap(16, ) : begining index in direction 2 for turbulent flow
c      imap(17, ) : ending   index in direction 2 for turbulent flow
c      imap(18, ) : begining index in normal direction for fmax 
c      imap(19, ) : ending   index in normal direction for fmax 
c      imap(20, ) : ending   index in normal direction for turb. flow 
c
c      the following items are added to facilitate the splitting of
c      the cfl3d input files, and are not part of the standard tlns3d
c      imap array
c
c      the following additions apply to face segments
c      imap(21, ) : number of additional data for 2000 series bc's
c      imap(22, ) : force flag for this segment
c      xmap( 1, ) : 1st additional data entry for 2000 series bc's
c      xmap( 2, ) : 2nd additional data entry for 2000 series bc's
c      xmap( 3, ) : 3rd additional data entry for 2000 series bc's
c      xmap( 4, ) : 4th additional data entry for 2000 series bc's
c      xmap( 5, ) : 5th additional data entry for 2000 series bc's
c
c      the following additions apply to blocks; ivisb(1-3, ) pertain 
c      to the turb model in the i,j,k directions, as before.
c      ivisb( 4, ) : ncg
c      ivisb( 5, ) : iem
c      ivisb( 6, ) : iadvance
c      ivisb( 7, ) : iforce
c      ivisb( 8, ) : ilamlo
c      ivisb( 9, ) : ilamhi
c      ivisb(10, ) : jlamlo
c      ivisb(11, ) : jlamhi
c      ivisb(12, ) : klamlo
c      ivisb(13, ) : klamhi
c      ivisb(14, ) : inewg
c      ivisb(15, ) : igridc
c      ivisb(16, ) : is
c      ivisb(17, ) : js
c      ivisb(18, ) : ks
c      ivisb(19, ) : ie
c      ivisb(20, ) : je
c      ivisb(21, ) : ke
c      ivisb(22, ) : idiag(i)
c      ivisb(23, ) : idiag(j)
c      ivisb(24, ) : idiag(k)
c      ivisb(25, ) : iflim(i)
c      ivisb(26, ) : iflim(j)
c      ivisb(27, ) : iflim(k)
c      ivisb(28, ) : ifds(i)
c      ivisb(29, ) : ifds(j)
c      ivisb(30, ) : ifds(k)
c      ivisb(31, ) : rkap0(i) (integer used to represent a real value)
c      ivisb(32, ) : rkap0(j)                    "
c      ivisb(33, ) : rkap0(k)                    "
c      ivisb(34, ) : iovrlp
c
c--------------------  boundary/operation type description ---------------------
c
c      nbctype    = imap(1 , )
c                 = 0  corresponds to an inner cut that maps to same block
c                 = 1  corresponds to a cut that maps to another block
c                 = 2  corresponds to a slip (inviscid) wall
c                 = 3  corresponds to a noslip (viscous) wall
c                 = 4  symmetry condition (imap(8) tells about which plane)
c                 = 5  downstream/extrapolation condition
c                 = 6  far-field condition (Riemann invariants)
c                 = 7  singular line
c                 = 8  fan-upstream engine boundary
c                 = 9  fan-downstream engine boundary
c
c***************************************************************************
c
c     ivisb(1)    = 1  viscous fluxes evaluated along i-direction
c     ivisb(2)    = 1  viscous fluxes evaluated along j-direction
c     ivisb(3)    = 1  viscous fluxes evaluated along k-direction
c
c***************************************************************************
c

c      print *,'is =',is
      read(imapfil,'(1x)')
      read(imapfil,*) idum
      if(idum.ne.nbloc) then
         write(iwrit,'(2x,"no. of blocks on map file ",i3,
     .         " do not match the no. of blocks on grid file",
     .          i3)') idum,nbloc
         call termn8(0,-6,ibufdim,nbuf,bou,nou)
      endif
c
c
c----------  read in imap array and (don't) shift due to ghost cells  --------
c
      nbltop = nbloc
      noldbl = nbloc
      nsgtop = 0
      do 100 ibloc = 1,nbloc
        nsubbl(ibloc) = 1
        idobl(ibloc) = ibloc
        idno(ibloc) = ibloc
        ijk(1,ibloc) = 1
        ijk(2,ibloc) = ni(ibloc)
        ijk(3,ibloc) = 1
        ijk(4,ibloc) = nj(ibloc)
        ijk(5,ibloc) = 1
        ijk(6,ibloc) = nk(ibloc)
        idseg(ibloc) = nsgtop + 1
        idnext(ibloc) = 0
        read(imapfil,'(1x)')
        read(imapfil,'(a79)') title1
        read(imapfil,*) nseg(ibloc),ivisb(1,ibloc),ivisb(2,ibloc),
     .                  ivisb(3,ibloc),(itrb(n,ibloc),n=1,7)
        read(imapfil,*)
        read(imapfil,*) (ivisb(n,ibloc),n=4,15)
        read(imapfil,*) (ivisb(n,ibloc),n=16,27)
        read(imapfil,*) (ivisb(n,ibloc),n=28,34)
        read(imapfil,*)
        read(imapfil,'(a79)') title2
        read(imapfil,'(a79)') title3
        do 101 iseg  = 1,nseg(ibloc)
          nsgtop = nsgtop + 1
          idbloc(nsgtop) = ibloc
          intrfc(nsgtop) = 0
c
          read(imapfil,'(1x)')
c         write(6,*) ibloc,iseg
          read(imapfil,*) idum,idum,(imap(n,nsgtop),n= 1,12)
          read(imapfil,*) (imap(n,nsgtop),n=13,20)
c         read(imapfil,*) (imap(n,nsgtop),n=13,20),val(nsgtop)
          read(imapfil,*) 
          read(imapfil,*)
          read(imapfil,*) (xdum(n,nsgtop),n=1,5)
          read(imapfil,*) (xdum(n,nsgtop),n=6,10)
          read(imapfil,*) (xdum(n,nsgtop),n=11,15)
          read(imapfil,*) (xdum(n,nsgtop),n=16,20)
          read(imapfil,*) (xdum(n,nsgtop),n=21,25)
          read(imapfil,*) (xdum(n,nsgtop),n=26,26)
          read(imapfil,*) (imap(n,nsgtop),n=21,23)
          read(imapfil,*)
c
          nxtseg(nsgtop) = nsgtop + 1
          iold(1,nsgtop) = imap(3,nsgtop)
          iold(2,nsgtop) = imap(4,nsgtop)
          iold(3,nsgtop) = imap(5,nsgtop)
          iold(4,nsgtop) = imap(6,nsgtop)
 101    continue
        nxtseg(nsgtop) = 0
c
c
 100  continue
c
      do 200 ibloc=1,nbloc
        loc = idseg(ibloc)
        do 201 iseg = 1,nseg(ibloc)
          if((imap(1,loc).eq.0).or.
     .       (imap(1,loc).eq.1)) then
            nbs  = imap(7,loc)
            nfs  = imap(8,loc)
            mns1 = min0(imap( 9,loc),imap(10,loc))
            mxs1 = max0(imap( 9,loc),imap(10,loc))
            mns2 = min0(imap(11,loc),imap(12,loc))
            mxs2 = max0(imap(11,loc),imap(12,loc))
            nfs  = abs(nfs)
            match = 0
            nloc = idseg(nbs)
            do 310 ns=1,nseg(nbs)
               if((imap(1,nloc).eq.0.or.imap(1,nloc).eq.1)
     .         .and.imap(2,nloc).eq.nfs) then
                 if((imap(3,nloc).eq.mns1 .and. 
     .               imap(4,nloc).eq.mxs1).and.
     .              (imap(5,nloc).eq.mns2 .and.
     .               imap(6,nloc).eq.mxs2))     then
                    match = nloc
                 end if
               end if
               nloc = nxtseg(nloc)
  310       continue
            if(match .eq. 0) then
              write(6,'(2x,"segment ",i3," of block ",i3,
     .           " does not match target info")') iseg,ibloc
              write(6,*) nbs,nfs,mns1,mxs1,mns2,mxs2
              call termn8(0,-6,ibufdim,nbuf,bou,nou)
            endif
            intrfc(loc) = match
          endif
          loc = nxtseg(loc)
 201    continue
 200  continue
c
      return
      end
