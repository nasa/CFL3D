      subroutine writ_buf(nbl,iunit,nou,bou,nbuf,ibufdim,myhost,myid,
     .                    mycomm,mblk2nd,maxbl)
c
c     $Id$
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
