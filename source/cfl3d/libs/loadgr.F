      subroutine loadgr(w,mgwk,lx,ly,lz,jindex,x,y,z,mdim,ndim,
     .                  idimg,jdimg,kdimg)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Load the proper grid from 1-d storage array to 2-d
c     work array.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension w(mgwk),x(mdim,ndim),y(mdim,ndim),z(mdim,ndim)
c      
c     patch surface is an i=constant surface
c
      if (jindex/10.eq.1) then
         if (jindex.eq.11) then
            i = 1
         else
            i = idimg
         end if
         do 10 j=1,jdimg
         do 10 k=1,kdimg
         l      = (i-1)*jdimg*kdimg+(k-1)*jdimg+(j-1)
         x(j,k) = w(lx+l)
         y(j,k) = w(ly+l)
         z(j,k) = w(lz+l)
   10    continue      
      end if
c      
c     patch surface is a j=constant surface
c
      if (jindex/10.eq.2) then
         if (jindex.eq.21) then
            j = 1
         else
            j = jdimg
         end if
         do 20 i=1,idimg
         do 20 k=1,kdimg
         l      = (i-1)*jdimg*kdimg+(k-1)*jdimg+(j-1)
         x(k,i) = w(lx+l)
         y(k,i) = w(ly+l)
         z(k,i) = w(lz+l)
   20    continue      
      end if      
c      
c     patch surface is a k=constant surface
c
      if (jindex/10.eq.3) then
         if (jindex.eq.31) then
            k = 1
         else
            k = kdimg
         end if
         do 30 i=1,idimg
         do 30 j=1,jdimg
         l      = (i-1)*jdimg*kdimg+(k-1)*jdimg+(j-1)
         x(j,i) = w(lx+l)
         y(j,i) = w(ly+l)
         z(j,i) = w(lz+l)
   30    continue      
      end if
      return
      end
