      subroutine ufree_c(ipoint)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  universal free...deallocates memory via fortran 77
c     pointer extension.
c
c     NOTE: in this special-purpose version, frees COMPLEX REAL (not
c     plin real) allocation, as well as integer allocation.
c
c     original version provied courtesy Steven J. Massey, AS&M Inc.
c***********************************************************************
c
      implicit complex(a-h,o-z)
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
