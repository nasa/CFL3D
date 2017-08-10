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
      subroutine rpatch(maxbl,maxxe,intmax,nsub1,windex,ninter,iindex,
     .                  nblkpt,dthetxx,dthetyy,dthetzz,jdimg,kdimg,
     .                  idimg,nou,bou,nbuf,ibufdim,myid)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Read in generalized-coordinate interpolation data
c     from a file
c***********************************************************************
c
#ifdef CMPLX
      implicit complex(a-h,o-z)
c
      integer stats
c
      real, dimension(:,:), allocatable :: windexr
#endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension windex(maxxe,2),iindex(intmax,6*nsub1+9),nblkpt(maxxe),
     .          dthetxx(intmax,nsub1),dthetyy(intmax,nsub1), 
     .          dthetzz(intmax,nsub1),jdimg(maxbl),kdimg(maxbl),
     .          idimg(maxbl)
c

#ifdef CMPLX
c     allocate real array in which to read patch data
c
      memuse = 0
      allocate( windexr(maxxe,2), stat=stats )
      call umalloc(maxxe*2,0,'windexr',memuse,stats)
#endif
c
c     read patched-grid connectivity file
c
      read(22) ninter
      if (ninter.gt.intmax) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),9220)ninter
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
 9220 format(52h too many blocks to interpolate; intmax should be at,
     .        7h least ,i4)

      do 1500 n=1,abs(ninter)
      read(22) iindex(n,1)
      nfb = iindex(n,1)
      read(22) (iindex(n,1+ll),ll=1,nfb)
      read(22) iindex(n,nfb+2)
      read(22) (iindex(n,nfb+2+ll),ll=1,nfb)
      read(22) iindex(n,2*nfb+3)
      read(22) iindex(n,2*nfb+4)
      read(22) iindex(n,2*nfb+5)
      read(22) iindex(n,2*nfb+6)
      read(22) iindex(n,2*nfb+7)
      read(22) iindex(n,2*nfb+8)
      read(22) iindex(n,2*nfb+9)
      nbl = iindex(n,nfb+2)
      lst  =iindex(n,2*nfb+5)
      len = lst + iindex(n,2*nfb+4) - 1
      if (n.eq.ninter .and. len.gt.maxxe) then 
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1)len
    1    format(/,1x,42hstopping...parameter maxxe is too small...,
     .               16hmust be at least,i6)
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
      read(22) (nblkpt(nnn),nnn=lst,len)
#ifdef CMPLX
      read(22) ((windexr(nnn,ll),nnn=lst,len),ll=1,2)
      do ll=1,2
         do nnn=lst,len
            windex(nnn,ll) = windexr(nnn,ll)
         end do
      end do
#else
      read(22) ((windex(nnn,ll),nnn=lst,len),ll=1,2)
#endif
 1500 continue
      ninter = abs(ninter)
c
c     if angular displacements are in the patch file, read them
c     (obtained with ioflag=2 in ronnie input); if not (ioflag=0 or 1)
c     set angular patch displacements to zero
c
      itest = 0
      do n=1,ninter
         nfb = iindex(n,1)
         do ll=1,nfb
#ifdef CMPLX
            read(22,end=999) dthetxr,dthetyr,dthetzr
            dthetxx(n,ll) = dthetxr
            dthetyy(n,ll) = dthetyr
            dthetzz(n,ll) = dthetzr
#else
            read(22,end=999) dthetxx(n,ll),dthetyy(n,ll),
     .                       dthetzz(n,ll)
#endif
         end do
      end do
      itest = 1
c
  999 continue
c
      if (itest .eq. 0) then 
         do n=1,ninter
            nfb = iindex(n,1)
            do ll=1,nfb
               dthetxx(n,ll) = 0.
               dthetyy(n,ll) = 0.
               dthetzz(n,ll) = 0.
            end do
         end do
      end if
c
c     set search range on from side to full index range
c     (added to be compatible with dynamic patch storage
c     - not needed in flow code for static patching)
c
      do 1700 n=1,ninter
      nfb = iindex(n,1)
      do 1700 ll=1,nfb
      nbl = iindex(n,ll+1)
      itop = (iindex(n,ll+nfb+2)-iindex(n,ll+nfb+2)/100*100)/10
      if (itop .eq. 1) then
         jd1 = jdimg(nbl)
         kd1 = kdimg(nbl)
      end if
      if (itop .eq. 2) then
         jd1 = kdimg(nbl)
         kd1 = idimg(nbl)
      end if
      if (itop .eq. 3) then
         jd1 = jdimg(nbl)
         kd1 = idimg(nbl)
      end if
      iindex(n,2*nfb+9+ll) = 1 
      iindex(n,3*nfb+9+ll) = jd1
      iindex(n,4*nfb+9+ll) = 1
      iindex(n,5*nfb+9+ll) = kd1
 1700 continue
#ifdef CMPLX
c
c     deallocate real array in which patch data was read
c
      deallocate(windexr)
#endif
c
      return
      end
