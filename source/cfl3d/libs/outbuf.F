      subroutine outbuf(str,iunit)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Write the output buffer (array).
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*1 str(120)
      do i=120,1,-1
         if(str(i).ne.' ') go to 10
      end do
   10 continue
      write(iunit,1) (str(k),k=1,i)
    1 format(120a1)
#ifdef CRAY_TIME
#else
#   ifdef IBM
      call flush_(iunit)
#   else
      call flush(iunit)
#   endif
#endif
      return
      end
