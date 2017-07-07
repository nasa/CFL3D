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
      subroutine readdat(bcdata,mdim,ndim,filname,bcdat,nou,bou,nbuf,
     .                   ibufdim,myid,mblk2nd,maxbl)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Read in auxiliary data arrays for 2000 series boundary
c     conditions
c***********************************************************************
c
#ifdef CMPLX
      implicit complex(a-h,o-z)
c
      integer stats
c
      real, dimension(:,:,:,:), allocatable :: bcdatar
#endif
c
      character*80 filname
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension bcdata(mdim,ndim,2,12),bcdat(12),mblk2nd(maxbl)
c

#ifdef CMPLX
c     allocate real array in which to read bc data
c
      memuse = 0
      allocate( bcdatar(mdim,ndim,2,12), stat=stats )
      call umalloc(mdim*ndim*2*12,0,'bcdatar',memuse,stats)
#endif
c
c     initialize data array for check in bc routine
c
      do 5 l=1,12
      do 5 ip=1,2
      do 5 m=1,mdim
      do 5 n=1,ndim
      bcdata(m,n,ip,l) = -1.e15
    5 continue
c
      if (filname.ne.'null') then
         iunit = 26
         open(unit=iunit,file=filname,form='formatted',
     .   status='old')
         write(11,*)
         write(11,'('' reading bc data file '',a80)') filname
         read(iunit,*)
         read(iunit,*) mdum,ndum,np
         if (mdum.ne.mdim .or. ndum.ne.ndim) then
            write(11,'('' stopping in readdat...data file'',
     .      '' has dimensions '',2i6)') mdum,ndum
            write(11,'('' dimensions set from input file'',
     .      '' are '',2i6)') mdim,ndim
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         if (np.lt.1 .or. np.gt.2) then
            write(11,'('' stopping in readdat...data file'',
     .      '' must have 1 or 2 planes of data - there are '',i6)') np
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         read(iunit,*) nvar
         if (nvar.gt.12 .or. nvar.lt.1) then
            write(11,'('' stopping in readdat...no more'',
     .      '' than 12 variables are allowed - there are '',i6)') nvar
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if 
#ifdef CMPLX
         read(iunit,*) ((((bcdatar(m,n,ip,l),m=1,mdim),n=1,ndim),
     .                                      ip=1,np),l=1,nvar)
         do l=1,nvar
            do ip=1,np
               do n=1,ndim
                  do m=1,mdim
                     bcdata(m,n,ip,l) = bcdatar(m,n,ip,l)
                  end do
               end do
            end do
         end do
#else
         read(iunit,*) ((((bcdata(m,n,ip,l),m=1,mdim),n=1,ndim),
     .                                     ip=1,np),l=1,nvar)
#endif
c
c        add 2nd plane of data if not read in (duplicate 1st plane)
c
         if (np.lt.2) then
            do 10 l=1,nvar
            do 10 m=1,mdim
            do 10 n=1,ndim
            bcdata(m,n,2,l) = bcdata(m,n,1,l)
   10       continue
         end if
c
         rewind(iunit)
         close(iunit)
c
      else
         do 30 l=1,12
         do 30 ip=1,2
         do 30 m=1,mdim
         do 30 n=1,ndim
         bcdata(m,n,ip,l) = bcdat(l)
   30 continue
      end if
#ifdef CMPLX
c
c     deallocate real array in which bc data was read
c
      deallocate(bcdatar)
#endif
c
      return
      end
