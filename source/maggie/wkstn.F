#if defined WKSTN_OFF
      function ismin_wkstn(n,x,incx)
#else
      function ismin(n,x,incx)
#endif
c
c     $Id$
c
c***********************************************************************
c     Purpose:  To return the index corresponding to the minimum value 
c     in an arrary; a replacement for the standard Cray function
c***********************************************************************
c
      dimension x(1+(n-1)*incx)
c
      ismin = 1
      xmin = x(1)
c
      do 100 nn=2,n,incx
      if(x(nn).lt.xmin) then
        ismin = nn
        xmin  = x(nn)
      end if
100   continue
c
#if defined WKSTN_OFF
c     must assign a value to the function name, even if a dummy name!
      ismin_wkstn=-999
#endif
c
      return
      end 
#if defined WKSTN_OFF
      function ismax_wkstn(n,x,incx)
#else
      function ismax(n,x,incx)
#endif
c***********************************************************************
c     Purpose:  To return the index corresponding to the maximum value
c     in an arrary; a replacement for the standard Cray function
c***********************************************************************
c
c
      dimension x(1+(n-1)*incx)
c
      ismax = 1
      xmax  = x(1)
c
      do 100 nn=2,n,incx
      if(x(nn).gt.xmax) then
        ismax= nn
        xmax  = x(nn)
      end if
100   continue
c
#if defined WKSTN_OFF
c     must assign a value to the function name, even if a dummy name!
      ismax_wkstn=-999
#endif
c
      return
      end
