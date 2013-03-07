      subroutine umalloc_r(n_words,intflag,text,memuse,status)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Keep track of memory allocated via allocate in F90
c
c     NOTE: in this special-purpose version, allocates memory
c     for REAL (not complex real) arrays, as well as integer arrays.
c
c     text    - text string containing allocatable array name 
c     intflag - flag to set variable type:
c              -1...real*4 variable
c               0...real variable
c               1...integer variable
c     memuse  - running total of memory that has been allocated
c     status  - flag indicating status of call to allocate: 0 for
c               successful allocation, otherwise the allocation was 
c               unsuccessful; gracefully abort
c***********************************************************************
c
#if defined DIST_MPI
#     include "mpif.h"
#endif
c
      parameter(len_i=4,len_r4=4)
c
      integer status
c
      character(*) text
c
      common /mydist2/ nnodes,myhost,myid,mycomm
c
#ifdef CRAY
      len_r = 8
#else
#  ifdef DBLE_PRECSN
      len_r = 8
#  else 
      len_r = 4
#  endif
#endif
c
      if (intflag.eq.1) then
        len = len_i
      else if (intflag.eq.0) then
        len = len_r
      else if (intflag.eq.-1) then
        len = len_r4
      end if
c
      n_bytes = n_words * len
c
c     check status of allocate call
      if (status.ne.0) then
         write(6,'(''stopping...failed trying to allocate memory '',
     .         ''for array '',a20)') text
         if (n_bytes.lt.1e10) then
            write(6,12) n_words, n_bytes
         else if (n_bytes.lt.1e15) then
           write(6,13) n_words, n_bytes
         end if
        write(6,14) memuse
c        call my_flush(6)
#ifdef DIST_MPI
         call MPI_ABORT(MPI_COMM_WORLD, myid, mpierror)
#else
         stop
#endif
      else
c        write(6,'(''allocating:'',a20,i15,'' bytes'')') text,n_bytes
c        call my_flush(6)
      end if
c
  12  format('requested ', i10, ' words, (',i11, ' bytes)')
  13  format('requested ', i15, ' words, (',i16, ' bytes)')
  14  format('prior to this call, total allocation was',
     .        i16, ' bytes')
c
      memuse  = memuse + n_bytes
c
      return
      end
