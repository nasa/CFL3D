      subroutine umalloc_r(ipoint,n_words,intflag,icall,memuse)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  universal malloc...allocates memory via fortran 77
c     pointer extension.
c
c     NOTE: in this special-purpose version, allocates memory
c     for REAL (not complex real) arrays, as well as integer arrays.
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
      parameter(len_i=4,len_r4=4)
c
      integer errcode,abort
c
      dimension dum(1)
c
      pointer (ipoint, dum)
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
        stop
      else
c        write(6,*) 'call',icall,' umalloc: allocated ',
c    .   n_bytes,' bytes'
      end if
c
  12  format('requested ', i10, ' words, (',i10, ' bytes)')
  13  format('requested ', i15, ' words, (',i15, ' bytes)')
  14  format('prior to this call, total allocation was',
     .        i15, ' bytes')
c
      memuse  = memuse + n_bytes
c
      return
      end
