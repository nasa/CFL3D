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
      subroutine getibk(blank,jdim,kdim,idim,nbl,intpts,nblpts,ibpnts,
     .                  iipnts,iitot,maxbl,iibg,kkbg,jjbg,ibcg,lig,
     .                  lbg,dxintg,dyintg,dzintg,iiig,jjig,kkig,
     .                  ibpntsg,iipntsg,myid,ibufdim,nbuf,bou,nou)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Read the output from MAGGIE (not the grids).
c***********************************************************************
c
#ifdef CMPLX
      implicit complex(a-h,o-z)
c
      integer stats
c
      real, dimension(:,:,:), allocatable :: blankr
      real, dimension(:), allocatable :: dxintgr
      real, dimension(:), allocatable :: dyintgr
      real, dimension(:), allocatable :: dzintgr
c
#endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension blank(jdim,kdim,idim)
      dimension intpts(4),iibg(iitot),kkbg(iitot),jjbg(iitot),
     .          ibcg(iitot),lig(maxbl),lbg(maxbl),ibpntsg(maxbl,4),
     .          iipntsg(maxbl)
      dimension dxintg(iitot),dyintg(iitot),dzintg(iitot),
     .          iiig(iitot),jjig(iitot),kkig(iitot)


#ifdef CMPLX
c
c     allocate real array in which to read MAGGIE data
c
      memuse = 0
      allocate( blankr(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'blankr',memuse,stats)
      allocate( dxintgr(iitot), stat=stats )
      call umalloc(iitot,0,'dxintgr',memuse,stats)
      allocate( dyintgr(iitot), stat=stats )
      call umalloc(iitot,0,'dyintgr',memuse,stats)
      allocate( dzintgr(iitot), stat=stats )
      call umalloc(iitot,0,'dzintgr',memuse,stats)
#endif
c
      idim1 = idim-1
      jdim1 = jdim-1
      kdim1 = kdim-1
c
c     initialize all iblank values to 1
c
      do 10 i=1,idim
      do 10 j=1,jdim
      do 10 k=1,kdim
      blank(j,k,i) = 1.0
 10   continue
c
      read(21) jchk,kchk,lchk
c
      if (jchk.ne.jdim1 .or. kchk.ne.kdim1 .or. lchk.ne.idim1) then
         write(11,*) '  mismatch in indices.....stopping in',
     .   ' getibk'
         write(11,*) '  jdim,kdim,idim=',jdim,kdim,idim
         write(11,*) '  jmax,kmax,lmax=',jchk,kchk,lchk
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
c
      read(21) ibpnts,intpts,iipnts,idum,idum
      ibpntsg(nbl,1) = intpts(1)
      ibpntsg(nbl,2) = intpts(2)
      ibpntsg(nbl,3) = intpts(3)
      ibpntsg(nbl,4) = intpts(4)
      iipntsg(nbl)   = iipnts 
c
      lsta = lig(nbl)
      lend = lsta+iipnts-1
      if (lend.gt.iitot) then
         write(11,33)
         write(11,34)lend,iitot
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
   33 format(' stopping in getibk - insufficient memory allocation')
   34 format(' lend,iitot = ',2i20)
c
#ifdef CMPLX
      read(21)(jjig(l),kkig(l),iiig(l),dyintgr(l),dzintgr(l),
     .         dxintgr(l),l=lsta,lend)
      do l=lsta,lend
         dyintg(l) = dyintgr(l)
         dzintg(l) = dzintgr(l)
         dxintg(l) = dxintgr(l)
      end do
#else
      read(21)(jjig(l),kkig(l),iiig(l),dyintg(l),dzintg(l),
     .         dxintg(l),l=lsta,lend)
#endif
c
      lsta = lbg(nbl)
      lend = lsta+ibpnts-1
      if (lend.gt.iitot) then
         write(11,33)
         write(11,34)lend,iitot
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
c
      read(21)(jjbg(l),kkbg(l),iibg(l),ibcg(l),l=lsta,lend)
c
#ifdef CMPLX
      read(21)(((blankr(j,k,i),j=1,jdim1),k=1,kdim1),i=1,idim1)
      do i=1,idim1
         do k=1,kdim1
            do j=1,jdim1
               blank(j,k,i) = blankr(j,k,i)
            end do
         end do
      end do
#else
      read(21)(((blank(j,k,i),j=1,jdim1),k=1,kdim1),i=1,idim1)
#endif
c
      do 8877 i=1,idim1
      isum = 0
      do 8876 k=1,kdim1
      do 8876 j=1,jdim1
      if (int(blank(j,k,i)).ne.1) isum = isum+1
 8876 continue
      if (isum.eq.0) then
c        write(15,*) ' for i=',i,' all iblank values are 1'
      else
c        write(15,*) ' iblank array...i=',i
         do 8866 k=1,kdim1
c        write(15,8833) (int(blank(j,k,i)),j=1,jdim1)
 8866    continue
      end if
 8877 continue
 8833 format(1x,71i1)
c
      nblpts  = 0
      do 9100 i=1,idim1
      do 9100 k=1,kdim1
      do 9100 j=1,jdim1
      if (blank(j,k,i).ne.1.) then
         nblpts = nblpts + 1
      end if
 9100 continue
c
      if (nbl.lt.maxbl) then
         lig(nbl+1) = lig(nbl)+iipntsg(nbl)
         lbg(nbl+1) = lbg(nbl)+ibpnts
      end if
c
#ifdef CMPLX
c
c     deallocate real array in which MAGGIE data was read
c
      deallocate(dxintgr)
      deallocate(dyintgr)
      deallocate(dzintgr)
      deallocate(blankr)
#endif
c
      return 
      end
