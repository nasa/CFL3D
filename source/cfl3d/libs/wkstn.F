#if defined WKSTN_OFF
      integer function isrcheq_wkstn(n,x,incx,target)
#else
      integer function isrcheq(n,x,incx,target)
#endif
c***********************************************************************
c     Purpose:  To find the first occurance in the array x that is equal
c     to target; a replacement for the standard Cray function
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension x(*)
      j=1
      isrcheq=0
      if(n.le.0) return
      if(incx.lt.0) j=1-(n-1)*incx
      do 100 i=1,n
      if(real(x(j)).eq.real(target)) go to 200
      j=j+incx
  100 continue
  200 isrcheq=i
#if defined WKSTN_OFF
c     must assign a value to the function name, even if a dummy name!
      isrcheq_wkstn=-999
#endif
      return
      end
