#   ifdef CMPLX
      complex function q8sdot(n,a,idum1,b)
#   else
      real function q8sdot(n,a,idum1,b)
#   endif
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Compute the dot product between two vectors.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension a(n),b(n)
      q8sdot = 0.
      do 10 i=1,n
   10 q8sdot = q8sdot+a(i)*b(i)
      return
      end
