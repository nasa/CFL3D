      subroutine q8vrev(n,a,ndum,b)
c***********************************************************************
c     Purpose:  Reverse the elements in an array.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension a(n),b(ndum)
      np   = n+1
      do 10 i=1,n
   10 b(i) = a(np-i)
      return
      end
