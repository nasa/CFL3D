      subroutine my_flush(iunit)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Generic routine to flush an output buffer
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      integer iunit, istat
c
#   ifdef CRAY
c     no flush on the crays....
#   endif
c
#   ifdef IBM
      call flush_(iunit)
#   endif
c
#   ifdef SGI
      call flush(iunit,istat)
#   endif
c
#   ifdef SUN
      call flush(iunit)
#   endif
c
#   ifdef HP
      call flush(iunit)
#   endif
c
#   ifdef ALPHA
      call flush(iunit)
#   endif
c
#   ifdef INTEL
c     no flush on the intel...
#   endif
c
#   ifdef PG
      call flush(iunit)
#   endif
c
#   ifdef LAHEY
      call flush(iunit)
#   endif
c
#   ifdef SUN
      call flush(iunit)
#   endif
c
#   ifdef GENERIC
      call flush(iunit)
#   endif

c
      return
      end

