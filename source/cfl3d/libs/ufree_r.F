      subroutine ufree_r(ipoint)
c***********************************************************************
c     Purpose:  universal free...deallocates memory via fortran 77
c     pointer extension.
c
c     NOTE: in this special-purpose version, frees REAL (not complex 
c     real) allocation, as well as integer allocation.
c
c     original version provied courtesy Steven J. Massey, AS&M Inc.
c***********************************************************************
c
      integer errcode,abort
c
      pointer (ipoint,dum)
c
#ifdef CRAY_TIME
c     cray_time implies cray
      abort   = 0
      errcode = 0
      call hpdeallc(ipoint,errcode,abort )
#else
#   ifdef IBM
      call free(%val(ipoint))
#   else
      call free(ipoint)
#   endif
#endif
      return
      end
