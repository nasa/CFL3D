#   ifdef CMPLX
      complex function q8smax(n,a)
#   else
      real function q8smax(n,a)
#   endif
c***********************************************************************
c     Purpose:  Find the maximum value in an array.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension a(n)
      q8smax = a(1)
      do 10 i=2,n
   10 if (real(q8smax).lt.real(a(i))) q8smax = a(i)
      return
      end
