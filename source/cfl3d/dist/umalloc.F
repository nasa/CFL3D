      subroutine umalloc(ipoint,n_words,intflag,icall,memuse)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  universal malloc...allocates memory via fortran 77
c     pointer extension.
c
c     icall   - integer counter to aid in debugging
c     intflag - flag to set variable type:
c              -1...real*4 variable
c               0...real variable
c               1...integer variable
c     memuse  - running total of memory that has been allocated
c
c     original version provied courtesy Steven J. Massey, AS&M Inc.
c     modified by r. t. biedron, January, 1999.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
#if defined DIST_MPI
#     include "mpif.h"
#endif
c
      parameter(len_i=4,len_r4=4)
c
      integer errcode,abort
c
      dimension dum(1)
c
      pointer (ipoint, dum)
c
      common /mydist2/ nnodes,myhost,myid,mycomm
c
#ifdef CRAY_TIME
c     cray_time implies cray
      len_r = 8
#else
#  ifdef DBLE_PRECSN
      len_r = 8
#  else 
      len_r = 4
#  endif
#endif
c
#ifdef CMPLX
      len_r = 2*len_r
#endif
c
      icall = icall+1
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
#ifdef CRAY_TIME
c     cray_time implies cray
      abort = 0
      errcode = 0
      nwords = n_words
#ifdef CMPLX
      if (intflag.eq.0) then
         nwords = nwords*2
      end if
#endif
      call hpalloc( ipoint, nwords, errcode, abort)
#else
#  ifdef IBM
      ipoint = malloc(%val(n_bytes))
#  else
      ipoint = malloc(n_bytes)
#  endif
#endif
c
      if (ipoint.eq.0) then
        write(6,'(a,i4)') 
     .  'stopping...malloc failed to allocate the memory on call ',
     .  icall
        if (n_bytes.lt.1e10) then 
          write(6,12) n_words, n_bytes
        else if (n_bytes.lt.1e15) then 
          write(6,13) n_words, n_bytes
        else 
          write(6,*) 'requested more than 10^(15) bytes'
        end if
        write(6,14) memuse
#ifdef CRAY_TIME then
#else
#   ifdef IBM then
        call flush_(6)
#   else
        call flush(6)
#   endif
#endif
#ifdef DIST_MPI then
        call MPI_ABORT(MPI_COMM_WORLD, myid, mpierror)
#else
        stop
#endif
      else
c        write(6,*) 'call',icall,' umalloc: allocated ',
c    .   n_bytes,' bytes'
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
