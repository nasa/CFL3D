      subroutine mreal(xyz,mdim1,ndim1,mdim2,ndim2,x,y,z)
c
c     $Id$
c
c***********************************************************************
c      Purpose: Unload grid data from transfer array xyz to x,y,z
c      work arrays
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension xyz(mdim1,ndim1,3)
      dimension x(mdim2,ndim2),y(mdim2,ndim2),z(mdim2,ndim2)
c
c     initialize the individual x,y,z arrays, since they are typically
c     larger that the corresponding sections in xyz (i.e. mdim2*ndim2 
c     vs. mdim1*ndim1)
c
      do m=1,mdim2
         do n=1,ndim2
            x(m,n) = 0.
            y(m,n) = 0.
            z(m,n) = 0.
         end do
      end do
c
      do m=1,mdim1
         do n=1,ndim1
            x(m,n) = xyz(m,n,1)
            y(m,n) = xyz(m,n,2)
            z(m,n) = xyz(m,n,3)
         end do
      end do
c
      return
      end
