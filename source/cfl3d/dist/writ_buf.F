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
#   ifdef FASTIO
      subroutine writ_buffast(nbl,iunit,nou,bou,nbuf,ibufdim,myhost,
     .                    myid,mycomm,mblk2nd,maxbl,idbugp)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Passes the internal buffer that resids on the processor
c     that owns block nbl to the host, and writes the contents of the
c     buffer to the specified unit number, iunit.
c     "Fast I/O" version written by Bob Bergeron, NASA Ames, 2005
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
      dimension istat(1)
c
      character*120 bou(ibufdim,nbuf)
      character*120 boutemp
      dimension nou(nbuf)
      dimension mblk2nd(maxbl)
      common/rjbdbgi/lunfio0
c
c
      iou = 1
      if (iunit.eq.11) iou = 1
      if (iunit.eq. 9) iou = 2
      if (iunit.eq.14) iou = 3
      if (iunit.eq.25) iou = 4
c
c     set baseline tag values
c
        myidchk = mblk2nd(nbl)
        if(myid.eq.myidchk.or.myid.eq.myhost) then
          if(myid.eq.myhost) then
             myidb=0
          else
             myidb=1
          endif
          if(myid.eq.myidchk) then
             write(lunfio0+myid,'(" ppflag ",i3,i6)')idbugp,nou(iou)
            if(nou(iou).gt.0) then
crjb...once per entry worker writes flag for reassembly code
c             write(lunfio0+myid,'(" ppflag ",i3)')idbugp
              do kou=1,nou(iou)
                 write(lunfio0+myid,'(a)')bou(kou,iou)
             enddo
            endif
            do kou=1,nou(iou)
             bou(kou,iou) = ' '
            end do
            nou(iou) = 0
crjb..myid.eq.myidchk
          else
c           if (nou(iou) .ge. ibufdim) then
c               write(iunit,'(''WARNING: internal buffer length'',
c    .          '' exceeded -  make parameter ibufdim > '',i6)') nou(iou)
c               write(iunit,'(''continuing, but you will be missing'',
c    .          '' some output data '')')
c           end if
crjb..myid.eq.myhost
crjb...flag for reassembly code
            write(11,'(" ppflag ",i3)')idbugp
            do kou=1,nou(iou)
             bou(kou,iou) = ' '
            end do
            nou(iou) = 0
          endif
        endif
      return
      end
#   else
      subroutine writ_buf(nbl,iunit,nou,bou,nbuf,ibufdim,myhost,myid,
     .                    mycomm,mblk2nd,maxbl)
c
c***********************************************************************
c     Purpose:  Passes the internal buffer that resids on the processor 
c     that owns block nbl to the host, and writes the contents of the
c     buffer to the specified unit number, iunit.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
#if defined DIST_MPI
#     include "mpif.h"
      dimension istat(MPI_STATUS_SIZE)
#endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension mblk2nd(maxbl)
c
      iou = 1
      if (iunit.eq.11) iou = 1
      if (iunit.eq. 9) iou = 2
      if (iunit.eq.14) iou = 3
      if (iunit.eq.25) iou = 4
c
#if defined DIST_MPI
c     set baseline tag values
c
      ioffset  = maxbl
      itag_nou = 1
      itag_bou = itag_nou + ioffset
c
      myidchk = mblk2nd(nbl)
      mytag = itag_nou + nbl
      if (myid.eq.myidchk) then
         call MPI_Send (nou, nbuf, MPI_INTEGER, myhost,
     .                  mytag, mycomm, ierr)
      else if (myid.eq.myhost) then
         call MPI_Recv (nou, nbuf, MPI_INTEGER, myidchk,
     .                  mytag, mycomm, istat, ierr)
      end if
c
#endif
      if (nou(iou).gt.0) then
#if defined DIST_MPI
         myidchk = mblk2nd(nbl)
         mytag = itag_bou + nbl
         nvals = 120*ibufdim*nbuf
         if (myid.eq.myidchk) then
            call MPI_Send (bou, nvals, MPI_CHARACTER, myhost,
     .                     mytag, mycomm, ierr)
         else if (myid.eq.myhost) then
            call MPI_Recv (bou, nvals, MPI_CHARACTER,
     .                     myidchk, mytag, mycomm, istat, ierr)
         end if
c
         if (myid.eq.myhost) then
#endif
         do kou=1,nou(iou)
            call outbuf(bou(kou,iou),iunit)
         end do
c
         if (nou(iou) .ge. ibufdim) then
            write(iunit,'(''WARNING: internal buffer length'',
     .      '' exceeded -  make parameter ibufdim > '',i6)') nou
            write(iunit,'(''continuing, but you will be missing'',
     .      '' some output data '')')
         end if
#if defined DIST_MPI
         end if
c
         if (myid.eq.myidchk .or. myid.eq.myhost) then
#endif
         do kou=1,nou(iou)
            bou(kou,iou) = ' '
         end do
         nou(iou) = 0
#if defined DIST_MPI
         end if
#endif
      end if
c
      return
      end
#   endif
