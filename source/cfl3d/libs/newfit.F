      subroutine newfit(ifits,ifit,icount)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Determine new polynomial fit for cells with stubborn
c     convergence.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      common /hist/ histry(4)
c
      if(icount.eq.1) then
        histry(1) = ifits
        iff = ifits
        do 10 n=2,4
        if(iff .eq. 4)iff = 0
        iff = iff+1
        histry(n) = iff
10      continue
      end if
c
      ifit = histry(icount+1)
      return
      end
